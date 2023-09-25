

function loadrnaseq(metafile, countfile)
    meta   = CSV.read(metafile, DataFrame)
    counts = CSV.read(countfile, DataFrame)
    meta, counts
end

"""
    filterdesign(design, counts; min_thresh, groupfield=:Group, samplefield=:Sample)

Ensures a min of min_thresh counts in all samples of a group, or at least `ns` samples if there are more than that.
"""
function filterdesign(design, counts; min_thresh=10, groupfield=:Group, samplefield=:Sample, ns=3, verbose=false)

    gdf = combine(groupby(design, groupfield), nrow => :Count, samplefield => Ref => :Samples)
    M = mapreduce(ind -> sum(Matrix(counts[!, ind]) .≥ min_thresh, dims=2) .≥ min(length(ind), ns), hcat, gdf.Samples)

    selind = vec(sum(M, dims=2)) .> 0
    verbose && println("Filtered ", sum(selind), " of ", length(selind), " : ", sum(selind)/length(selind))
    
    fields = [first(names(counts)) ; design[!, samplefield]]
    design, counts[selind, fields]

end



"""
    calcmeancounts(design, counts)

Calculates mean counts
"""
function calcmeancounts(design, counts)
    designgroup = combine(groupby(design, [:Stage, :GenoType, :Group]), nrow => :count, :Sample => Ref => :Samples)
    sort!(designgroup, [order(:GenoType, rev=true), :Stage])
    [DataFrame(Gene=counts.Gene) DataFrame(mapreduce(s -> mean(Matrix(counts[!, s]), dims=2), hcat, designgroup.Samples), designgroup.Group)]
end


function decounttable(counts, deg, lt=log2(0))
    deind = deg.sig .& (abs.(deg.log2FoldChange) .> lt)
    degenes = unique(deg.Gene[deind])
    @subset(counts, :Gene .∈ Ref(Set(degenes)))
end



function saveresults(res; projdir=getprojectdir(), resultsdir="results")

    files = Dict(:deg => "znf808_deseq.tsv",
                 :ncounts => "znf808_normcounts.tsv",
                 :design  => "znf808_tested_meta.tsv",
                 :ent     => "znf808_enrichr.tsv")
    respath = joinpath(projdir, resultsdir)
    mkpath(respath)

    for (s, f) in files
        fp = CSV.write(joinpath(respath, f), res[s], delim='\t')
        println("Writing: ", fp)
    end
end

function loadresults(;projdir=getprojectdir(), resultsdir="results")

    files = Dict(:deg => "znf808_deseq.tsv",
                 :ncounts => "znf808_normcounts.tsv",
                 :design  => "znf808_tested_meta.tsv",
                 :ent     => "znf808_enrichr.tsv")
    respath = joinpath(projdir, resultsdir)
    data = Dict(s => CSV.read(joinpath(respath, f), DataFrame) for (s, f) in files)
    sdf = (stages=unique(data[:design].Stage),)
    merge(sdf, (; data...))
end




function ipsfc(deg, design, ncounts)
    
    siggenes = @subset(deg, :sig).Gene |> unique |> Set
        
    ips_run = @subset(design, :Run .== "R3", :Stage .∈ Ref(["D0", "DE", "S2", "S3", "S4"]))
    
    RC = @subset(ncounts, :Gene .∈ Ref(siggenes))[!, ["Gene" ; ips_run.Sample]]
    
    ipsfc = @chain RC begin
        stack(ips_run.Sample, variable_name=:Sample, value_name=:count)
        innerjoin(ips_run[!, [:Sample, :Stage, :GenoType]], on=:Sample)
        unstack([:Gene, :Stage], :GenoType, :count)
    
    end
    
    @assert all(!ismissing, ipsfc.KO)
    @assert all(!ismissing, ipsfc.WT)
    ipsfc.ips_l2fc = @with ipsfc log2.(:KO .+ 1) .- log2.(:WT .+ 1)
    
    deg_ips = innerjoin(deg, ipsfc[!, [:Gene, :Stage, :ips_l2fc]], on=[:Gene, :Stage])
    @subset(deg_ips, :sig)
end


function ips_comparison(deresults, design_ips_full, ncounts)
    deg_ips = ipsfc(deresults.deg, design_ips_full, ncounts);

    plot()
    hline!([0], c=:black, ls=:dash, lab="")
    @with stack(@subset(deg_ips, :Stage .!= "S4"), [:log2FoldChange, :ips_l2fc], [:Gene, :Stage, :dsig]) groupedboxplot!(:Stage, :value,
        group=(ifelse.(:dsig .== 1, "Activated", "Repressed"), ifelse.(:variable .== "log2FoldChange", "ZNF808 KO", "iPSC")), outliers=false, c=^([:orange :darkorange2 :steelblue1 :steelblue]))
    pb = plot!(leg=:bottomleft, ylims=(-8, 6), ylabel="log2 Fold Change", fontfamily="helvetica")

    phs = [@with @subset(deg_ips, :Stage .== st) scattline(:log2FoldChange, :ips_l2fc, label=st, framestyle=^(:box), ylims=(-11, 11), aspect_ratio=^(:equal)) for st in ["D0", "DE", "S2", "S3"]]
    [plot!(p, yformatter=y -> "", left_margin=-5mm) for p in phs[[2, 4]]]
    [plot!(p, xformatter=x -> "", bottom_margin=-4mm) for p in phs[[1, 2]]]
    sw(x) = x/sum(x)
    lt = @layout [a{0.55w} [b c ; d e]]
    
    plot(pb, phs..., layout=lt, size=(700*1.1, 300*1.1), fontfamily="helvetica", grid=false)

end
