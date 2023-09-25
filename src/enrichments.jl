"""
    enrichments(deresult, laser_dorsal_hep, rois)

Main enrichment pipeline.
"""
function enrichments(deresult, laser_dorsal_hep, enrichrgenes)
    ent = enrichrtable(deresult.deg);

    ## recalc enrichments for selected genesets to ensure proper background is reflected
    enset = @subset(ent, :GeneSet .∈ Ref(Set(keys(enrichrgenes))))
    enc = recalcenrichr(deresult.deg, enset, enrichrgenes)
    len = laserenrichment(deresult.deg, laser_dorsal_hep);
    #proxdf = proxde(deresult.deg, rois, xp = range(0, 8, length=50));
    merge(deresult, (ent=ent, enfisher=enc, len=len))#, proxdf=proxdf))
end

function loadgeneset(file)
    geneset = Dict{String, Vector{String}}()
    io = open(file) |> GzipDecompressorStream
    for line in eachline(io)
        fields = split(line, '\t', keepempty=false)

        if occursin(r",[0-9]", fields[2])
            geneset[first(fields)] = replace.(fields[2:end], r",[0-9\.]*$" => "")
        else
            geneset[first(fields)] = fields[2:end]
        end
    end
    close(io)
    geneset
end


function loadenrichrgenesets(enrichr_folder="enrichr", projdir=getprojectdir())
    dir = joinpath(projdir, "data", enrichr_folder)
    files = glob("*.txt.gz", dir)
 
    samples = replace.(basename.(files), ".txt.gz" => "")
    Dict(s => loadgeneset(f) for (s, f) in zip(samples, files))
end

enrichr_genesets() = [ "KEGG_2019_Human", "BioPlanet_2019", "WikiPathways_2019_Human", "GO_Molecular_Function_2018", "GO_Biological_Process_2018", "GO_Cellular_Component_2018", "Human_Gene_Atlas"]


function enrichrtable(deg;  genesets=enrichr_genesets(), pt = 0.05)

    stages = unique(deg.Stage)
    
    ens = DataFrame[]
    n = length(stages)*2*length(genesets)

    p = Progress(n, desc="Enrichr API   : ")
    for st in stages, dir in [1, -1]
        ind = (deg.Stage .== st) .& (deg.dsig .== dir)
        genes = deg.GeneName[ind]
        config = load_genelist(genes)
        for gs in genesets
            et = enrichment(config, gs)
            et[!, :Stage] .= st
            et[!, :Direction] .= ifelse(dir == 1, "Up", "Down")
            push!(ens, et)
            next!(p)
        end
    end
    ent = reduce(vcat, ens)
    ent[!, :NumGenes] =  length.(ent.Genes)


    ent
end

function recalcenrichr(deg, ent, enrichrgenes) 

    tg = unique(ent[!, [:GeneSet, :Term]])
    df = recalcenrichr(tg.GeneSet, tg.Term, deg.GeneName, deg.Stage, deg.dsig, unique(deg.Stage), [-1, 1], enrichrgenes) 
    df.Direction = ifelse.(df.Direction .== 1, "Up", "Down")
    
    ens = names(ent)
    edf = setdiff(names(df), ens)
    leftjoin(ent, df, on=[:GeneSet, :Term, :Stage, :Direction])[!, [ens ; edf]]
end


function recalcenrichr(terms, genesets, genes, stages, dsig, ustages, uds, enrichrgenes)
   
    df = DataFrame(GeneSet=String[], Term=String[], Stage=String[], Direction=Int[], FG=Int[], FisherPvalue=Float64[], OR=Float64[])
    sti = [stages .== st for st in ustages]
    @showprogress "Recalc Enrichr:" for (t, g) in zip(terms, genesets)
        gs = Set(enrichrgenes[t][g])
        gind = uppercase.(genes) .∈ Ref(gs)
        for (us, si) in zip(ustages, sti), d in uds
            gi = gind[si]
            di = dsig[si].== d
            fg = sum(gi .& di)
            p, or = ProximityEnrichment.hypertest(gi, di)    
            push!(df, (t, g, us, d, fg, p, or))
        end
    end
    
    DataFrames.transform!(groupby(df, [:GeneSet, :Stage]), :FisherPvalue => (x -> MultipleTesting.adjust(Float64.(x), MultipleTesting.BenjaminiHochberg())) => :FDR)    
    df[!, :FDRAll] = MultipleTesting.adjust(df.FisherPvalue, MultipleTesting.BenjaminiHochberg())
    df
end



termind(ent, term) = termind(ent, term, "")
termind(ent, term::Tuple{V, T}) where {V, T} = termind(ent, term[1], term[2])

function termind(ent, term, gs)

    ind = occursin.(term, ent.Term)
    if gs != ""
        ind .&= ent.GeneSet .== gs
    end


    if (length(unique(ent.Term[ind])) > 1) || (length(unique(ent.GeneSet[ind])) > 1)
        error("Multiple terms: $gs $term\n$(unique(ent.Term[ind]))\n$(unique(ent.GeneSet[ind]))")
    end
    ind
end

function genesetind(ent, gs, term)
    ind = (ent.GeneSet .== gs) .& occursin.(term, ent.Term)
    if length(unique(ent.Term[ind])) > 1
        error("Multiple terms: $gs $term\n$(unique(ent.Term[ind]))")
    end
    ind
end




function saveenrichr(enfisher, stagelabels; projdir=getprojectdir())
    enf =  sort(enfisher, [order(:GeneSet, by=x->ifelse(x == "Human_Gene_Atlas", "A", x)), :FDR])
    enf.Genes = join.(enf.Genes, ", ")
    enf.Stage = getindex.(Ref(stagelabels), enf.Stage)
    enf.Direction = ifelse.(enf.Direction .== "Up", "Activated", "Repressed")
    filepath = joinpath(projdir, "results", "znf808_enrichr_genesets.csv.gz")

    mkpath(dirname(filepath))
    CSV.write(filepath, enf[!, Not([:FDRAll])], compress=true)
    enf
end

#### Comparison to LASER RNA-seq


function loadlaserde(file="laser_dorsal_panc_vs_hep_coords.csv.gz" ; projdir=getprojectdir())
    filepath = joinpath(projdir, "data", file)

    dp_hc = CSV.read(filepath, DataFrame, skipto=5, header=4)
    dp_hc = dp_hc[!, 1:7]
    rename!(dp_hc, [:GeneID, :log2FoldChange, :pvalue, :FDR, :genetype, :genestatus, :GeneName])
    dropmissing!(dp_hc)
    sort!(dp_hc, :FDR)
    dp_hc.LSD = ifelse.(dp_hc.FDR .< 0.05, ifelse.(dp_hc.log2FoldChange .> 0, 1, -1), 0)
    dp_hc
end

function laserenrichment(deg, lasers, codes = Dict(-1 => "D", 1 => "U"))

    rdf = DataFrame(Stage=String[], Label=String[], count=Int[], pvalue=Float64[], OR=Float64[])
    for df in groupby(deg, :Stage)
        jt = innerjoin(df, lasers[!, [:GeneName, :LSD]], on=:GeneName)
    
        lb = String[]
        cc = Int[]
        pv = Float64[]
        or = Float64[]
        
        for dd in [-1, 1], ld in [-1, 1]
            indA = jt.dsig .== dd
            indB = jt.LSD  .== ld
            pv, or = ProximityEnrichment.hypertest(indA, indB)

            push!(rdf, (df.Stage[1], string(codes[dd], codes[ld]), sum(indA .& indB), pv, or))
        end
    end
    rdf
    rdf.FDR = MultipleTesting.adjust(rdf.pvalue, MultipleTesting.BenjaminiHochberg());
    rdf
end



############## Proximity enrichments

function midmember(x)
    si = sortperm(x)
    x[si[cld(length(x), 2)]]
end

function loadtss(file=""; projdir=getprojectdir())
    if isempty(file)
        file = joinpath(projdir, "data", "gencode.v36lift37.tss.tsv.gz")
        @assert isfile(file)
    end
    tss = CSV.read(file, DataFrame)

    ### need a strategy for selecting different isoform TSS, ideally would do this based on isoform expression data
    ### however, here will select the median 

    genetss = combine(groupby(tss, [:Gene, :GeneID, :GeneName]), nrow => :count, :chrom => first => :chrom, :strand => first => :strand, :TSS => midmember => :TSS)
    sort!(genetss,[:chrom, :TSS])
    genetss
end

function annotatetss(deg, genetss)
    genefield = ifelse(:Gene ∉ propertynames(deg), :GeneName, :Gene)
    jt = leftjoin(deg, genetss[!, [genefield, :chrom, :TSS, :strand]], on=genefield)
    @assert !any(ismissing, jt.TSS)
    dropmissing!(jt)
    jt
end



tableintervals(df) = [Interval(c, s, e) for (c, s, e) in zip(df.chrom, df.start, df.stop)]

# function load_rois(; projdir=getprojectdir())
#     mer11 = CSV.read(joinpath(projdir, "data", "hg19.repeats.mer11.bed"), DataFrame)
    
#     #mer11clusters = CSV.read(joinpath(projdir, "results", "merclusters_all.bed"), DataFrame, header=[:chrom, :start, :stop, :Cluster, :ZNF808])
    

#     rors = [mer11, hnf4pro, mer11clusters]
#     labels = ["MER11", "Hepatic HNF4A", "MER11CLAll"]

#     for c in sort(unique(mer11clusters.Cluster))
#         push!(rors, @subset(mer11clusters, :Cluster .== c))
#         push!(labels, string("MER11CL", c))

#     end
    
#     for r in rors
#         sort!(r, [:chrom, :start, :stop])
#     end
#     intervals = tableintervals.(rors)
#     Dict(l => i for (l, i) in zip(labels, intervals))
# end

function mer11_cluster_rois(mer11, label="MER11", field=:Cluster, cln=2)
    cs = unique(mer11[!, field])

    labels = [label ; string.(label, "_CL", lpad.(cs, cln, "0"))]
    ms = sort(mer11, [:chrom, :start, :stop])
    inds = [[trues(size(mer11, 1))] ; [ms[!, field] .== c for c in cs]]

    Dict(l => tableintervals(ms[i, :]) for (l, i) in zip(labels, inds))
end

function proxde(deg, rois; xp = range(0, 8, length=50))

    proxdf = DataFrame(ROI=String[], Stage=String[], Direction=String[], x=Float64[], pvalue=Float64[], or=Float64[], count=Int[])
    n = length(rois)*length(unique(deg.Stage))
    p = Progress(n, "Proximity Enrichments: ")

    
    D0genes = @subset(deg, :Stage .== "D0", :dsig .!= 0).Gene |> Set
    DEgenes = @subset(deg, :Stage .== "DE", :dsig .!= 0).Gene |> Set
    
    stage_ex_genes = Dict("D0" => DEgenes, "DE" => D0genes) ## exclude D0 genes from DE and vice versa
    
    

    for (label, roi) in rois
        for df in groupby(deg, :Stage)
            sdf = sort(df, [:chrom, :TSS])
            geneintervals = @with sdf Interval.(:chrom, :TSS, :TSS, first.(:strand))
    


            up_ind = sdf.dsig .==  1
            dn_ind = sdf.dsig .== -1

            pu = proxenrich(xp, geneintervals, roi, up_ind)
            pd = proxenrich(xp, geneintervals, roi, dn_ind)
            
            
            udf = DataFrame(ROI=label, Stage=df.Stage[1], Direction="Up",   x=pu.xg, pvalue=pu.pvalue, or=pu.or, count=pu.count)
            ddf = DataFrame(ROI=label, Stage=df.Stage[1], Direction="Down", x=pd.xg, pvalue=pd.pvalue, or=pd.or, count=pd.count)

            append!(proxdf, udf)
            append!(proxdf, ddf)

            next!(p)
        end


        
    end
    proxdf
end