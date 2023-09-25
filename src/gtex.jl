
function loadgtexgct(metafilename="gtex_meta.tsv.gz", tpmfilename="gtex_znf808.txt.gz", projdir=getprojectdir())
    metafile = joinpath(projdir, "data", metafilename)
    tpmfile = joinpath(projdir, "data", tpmfilename)
    
    meta = CSV.read(metafile, DataFrame)

    io = open(tpmfile) |> GzipDecompressorStream 
    
    verline = readline(io)
    
    rows, cols = parse.(Int, split(readline(io)))
    header = split(readline(io))
    
    M = zeros(rows, cols)
    
    genes = String[]
    trans = String[]
    
    for i = 1:rows
        line = readline(io)
        fields = split(line)
        push!(trans, fields[1])
        push!(genes, fields[2])
        
        for j = 1:cols
           M[i, j] = parse(Float64, fields[j+2])
        end
        
    end
    
    close(io)
    
    df = [DataFrame([trans, genes], header[1:2]) DataFrame(M, header[3:end])]
    sdf = stack(df, names(df)[3:end], variable_name=:SampleID, value_name=:TPM)

    
    tpm = innerjoin(meta, sdf, on=:SampleID)  ### isoform level
    genetpm = combine(groupby(tpm, [:SampleID, :Tissue, :gene_id]), :TPM => sum => :TPM) ## gene level

    tissues = sort(unique(genetpm.Tissue))
    
    tissueorder = Dict(t => i for (i, t) in enumerate(tissues));
    mediantissues = sort(combine(groupby(genetpm, :Tissue), :TPM => median => :Median), :Median, rev=true).Tissue
    mediantissueorder = Dict(t => i for (i, t) in enumerate(mediantissues))
    genetpm.TissueOrder = getindex.(Ref(tissueorder), genetpm.Tissue);
    genetpm.MedianOrder = getindex.(Ref(mediantissueorder), genetpm.Tissue);

    tpm.TissueOrder = getindex.(Ref(tissueorder), tpm.Tissue);
    tpm.MedianOrder = getindex.(Ref(mediantissueorder), tpm.Tissue);



    (tpm=tpm, genetpm=genetpm)
end

function gtexboxplot(genetpm; rev=false, orderfield=:MedianOrder)
    plot(grid=false)
    tissues = sort(unique(genetpm[!, [:Tissue, orderfield]]), orderfield, rev=rev).Tissue

    cc = cgrad(:viridis, length(tissues), categorical=true)
    order = genetpm[!, orderfield]
    order = ifelse(rev, maximum(order) .- order .+ 1, order)

    boxplot!(order, genetpm.TPM*20, group=genetpm[!, orderfield], c=hcat(cc...), marker=(stroke(0), 0.5, 2), size=(1000,600),  lab="")
    plot!(xticks=(1:length(tissues), tissues), xrotation=45, bottom_margin=50mm, right_margin=10mm, fontfamily="helvetica", ylabel="TPM", title="GTEx ZNF808 Expression", left_margin=5mm, xlabel="Tissue")
end


function load_gtex(; gtex_liver_quantile_file="gtex_liver_nonliver_quantiles.tsv.gz", gtex_median_file="gtex_median.tsv.gz", tdo2_file="gtex_v8_tdo_tpm.tsv.gz", projdir=getprojectdir())
    gtex_quants = CSV.read(joinpath(projdir, "data", gtex_liver_quantile_file), DataFrame)
    median = CSV.read(joinpath(projdir, "data", gtex_median_file), DataFrame)
    tdo2   = CSV.read(joinpath(projdir, "data", tdo2_file), DataFrame)

    liver_exclusive = @subset(gtex_quants, :LiverLower .> :MaxNonLiverU)

    gtex = (; liver_exclusive, median, tdo2)

    gtex
end



function gtex_enrichments(gtex, deresults)

    gtex_enrich = DataFrame(Stage=String[], Direction=String[], count=Int[],exp=Float64[], pvalue=Float64[], OR=Float64[])
    for dir in [1, -1]
        for st in deresults.stages
            df = @subset(deresults.deg, :Stage .== st)
            df[!, :GTEXLiver] = df.GeneName .∈ Ref(Set(gtex.liver_exclusive.GeneName))
            count = sum((df.dsig .== dir) .& (df.GTEXLiver))
            exp = sum(df.dsig .== dir)*sum(df.GTEXLiver)/size(df, 1)
            p, or = ProximityEnrichment.hypertest(df.dsig .== dir, df.GTEXLiver)
            push!(gtex_enrich, (st, ifelse(dir == 1, "A", "R"), count, exp, p, or))
            
        end
    end
    gtex_enrich
end