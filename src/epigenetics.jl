
#### Loading ChIP data

function loadchipmeta(; projdir=getprojectdir(), file="chipseq_meta.tsv")
    filepath = joinpath(projdir, "data", file)
    meta = CSV.read(filepath, DataFrame)
    
    meta = @subset(meta, :Chip .!= "IN")
    sort!(meta, [order(:Run), order(:Chip, rev=true), order(:GenoType, rev=true), :Stage])

    meta.Index = 1:size(meta, 1)
    meta
end

function loadepigenetics(; projdir=getprojectdir(), file="chipseq_meta.tsv")
    meta = loadchipmeta(projdir=projdir, file=file)
    if "SignalFile" ∈ names(meta)
        FM = load_frag_matrix.(meta.SignalFile)
    else
        FM = nothing
    end
    meta, FM
end

function loadznfsmeta(; projdir=getprojectdir(), file="znfsmeta.tsv")
    filepath = joinpath(projdir, "data", file)
    znfs = CSV.read(filepath, DataFrame)
    znfs
end

function loadmer11(; projdir=getprojectdir(), file="MER11-ABCD.bed.gz")
    filepath = joinpath(projdir, "data", file)
    mer11 = CSV.read(filepath, DataFrame, header=[:chrom, :start, :stop, :name, :score, :strand, :fam])
    sort!(mer11, [:chrom, :start, :stop])
    mer11
end

function loadmer11znfs(; projdir=getprojectdir())
    files = glob("GSM*.gz", joinpath(projdir, "data", "mer11tfs"))
    peakdict = Dict(split(basename(f), "_")[2] => CSV.read.(f, DataFrame, header=[:chrom, :start, :stop, :name, :score, :strand]) for f in files);
    ### load ZNF578 from chip-atlas
    peakdict["ZNF578"] = CSV.read(joinpath(projdir, "data", "mer11tfs", "SRX4745235.05.bed.gz"), DataFrame, header=[:chrom, :start, :stop, :name, :score, :strand, :FC, :nlog10p, :nlog10q, :summit])
    peakdict
end



#### Recluster peaks from Michael Imbeault
function reordercluster(KMK, nr = 1:KMK.k)
    
    KSI = KMK.KSI[nr]
    KSSI = sortperm(KSI)
    A = KSSI[KMK.KM.assignments]
    SI = sortperm(A)
    (KSI = KSI, A=A, SI=SI, k=KMK.k, KM=KMK.KM, X=KMK.X)
end

zs(X) = (X .- mean(X, dims=2))./std(X, dims=2)

function clusterpeaks(deresults, mer11, file="ZNF808-K27ac-clusters-jan-22_broad_v4.csv.gz"; projdir=getprojectdir(), ks=[5])

    mcl = CSV.read(joinpath(projdir, "data", file), DataFrame)
    ("Column1" ∈ names(mcl)) && rename!(mcl, :Column1 => :Loc)
    ("end" ∈ names(mcl)) && rename!(mcl, :end => :stop)
    ("chr" ∈ names(mcl)) && rename!(mcl, :chr => :chrom)
    
    sort!(mcl, [:chrom, :start])
    
    mers = sort(unique(mer11.name))
    [markintersection!(mcl, @subset(mer11, :name .== m), m) for m in mers]

    ## separate into K27/K9 WT/KO for normalisation prior to clustering
    MHK27 = mcl[!, r"H1[DS][0-9E]K27"];
    MHK9  = mcl[!, r"H1[DS][0-9E]K9"];
    MKK27 = mcl[!, r"s4C[DS][0-9E]K27"];
    MKK9  = mcl[!, r"s4C[DS][0-9E]K9"];
    
    MK27 = [Matrix(MHK27) Matrix(MKK27)]
    MK9  = [Matrix(MHK9) Matrix(MKK9)]

    MMS = [MK9 MK27]
    ZMK27 = zs(MK27)
    ZMK9 = zs(MK9)
    ZMS = [ZMK9 ZMK27];


    ### calculate silhouettes to determine K
    ## reveals that k = 5 offers a good balance of capturing behaviours and generalisation.
    sks = 2:20
    SZMS = [silhouettes(kmeansorder(ZMS', k).KM, pairwise(Euclidean(), ZMS')) for k = sks]


    silhouetteplot = plot()
    plot!(sks, median.(SZMS), marker=:auto, lab="All", xticks=sks, xlabel="Cluster Size", ylabel="Median Silhouette Score")
    


    ### setup dict to store clusterings, in this instance only default is ks = [5] and only k = 5 is calculated
    kd = Dict(string("All_CL", k) => kmeansorder(ZMS', k) for k = ks)
    kd["All_CL5"] = reordercluster(kd["All_CL5"], [1, 3, 2, 4, 5]) ## reordering to ensure similar K27 behaviours appear consecutively
    
    for (k, KMK) in kd
       mcl[!, k] = KMK.A 
    end
    mcl.Cluster5Label = replace(mcl.All_CL5, 1 => 1.0, 2 => 1.5, 3 => 2, 4 => 3, 5 => 4)
    mcl.Cluster5Map = Int.(floor.(mcl.Cluster5Label))
    mcl.Cluster5Order = kd["All_CL5"].SI

    #### proximity enrichments
    ## first MI's region labels
    # mclrois = Dict(first(df.cluster) => Interval.(df.chrom, df.start, df.stop) for df in groupby(mcl, :cluster))
    mclrois = Dict{String, Vector{Interval{Nothing}}}()
    mclrois["D0A"] = @with @subset(mcl, occursin.("D0", :cluster)) Interval.(:chrom, :start, :stop);
    mclrois["DEA"] = @with @subset(mcl, occursin.("DE", :cluster)) Interval.(:chrom, :start, :stop);
    mclrois["S2A"] = @with @subset(mcl, occursin.("S2", :cluster)) Interval.(:chrom, :start, :stop);
    mclrois["S2A"] = @with @subset(mcl, occursin.("S2", :cluster)) Interval.(:chrom, :start, :stop);
    mclrois["S3A"] = @with @subset(mcl, occursin.("S3", :cluster)) Interval.(:chrom, :start, :stop);


    ## regions for each kmeans cluster
    for k in ks
        for i = 1:k
            ind_all = mcl[!, string("All_CL", k)] .== i
            mclrois["All_CL$(k)_$(i)"] = @with mcl[ind_all, :] Interval.(:chrom, :start, :stop)
        end
    end

    ### combinations of cluster 5
    mclrois["All_GCL5_12"] = @with @subset(mcl, :All_CL5 .∈ Ref([1, 2])) Interval.(:chrom, :start, :stop)
    mclrois["All_GCL5_123"] = @with @subset(mcl, :All_CL5 .∈ Ref([1, 2, 3])) Interval.(:chrom, :start, :stop)
    mclrois["All_GCL5_34"]  = @with @subset(mcl, :All_CL5 .∈ Ref([3, 4])) Interval.(:chrom, :start, :stop)

    ### All MER elements that lose k9 ad gain k27
    mclrois["MER_ALL"] = @with mcl Interval.(:chrom, :start, :stop)

    ### calculate prox enrichments
    proxmcl = proxde(deresults.deg, mclrois, xp=range(0, 8, length=75));
    ### summary table
    proxsum = @chain proxmcl begin
        @subset(:Direction .== "Up")
        groupby([:ROI, :Stage])
        combine(:pvalue => minimum => :pvalue)
        unstack(:ROI, :Stage, :pvalue)
    end


    ### annotate chipatlas and MER11 znfs
    # [markintersection!(mcl, p, c)  for (c, p) in chipatlaspeaks]
    # [markintersection!(mcl, p, c)  for (c, p) in znfpeaks];


    ### mark on MER11 elements to test them as a baseline
    # [markintersection!(mer11, p, c)  for (c, p) in chipatlaspeaks]
    # [markintersection!(mer11, p, c)  for (c, p) in znfpeaks];

    # annotatecol!(mer11, mcl, :All_CL5, :All_CL5);

    ### Association between annotations and clusters
    #labels = mapreduce(collect ∘ keys, vcat, [znfpeaks, chipatlaspeaks])
    clusterassociations = DataFrame(DataSet=String[], Cluster=Int[], Factor=String[], OddsRatio=Float64[], FisherPvalue=Float64[])

    # @showprogress for (cl, KMK) in kd, k = 1:KMK.k, l in labels
    #     p, or = ProximityEnrichment.hypertest(mcl[!, cl] .== k, mcl[!, l])
    #     push!(clusterassociations, (cl, k, l, or, p))
                
    #     p, or = ProximityEnrichment.hypertest(mer11[!, cl] .== k, mer11[!, l])
    #     push!(clusterassociations, (string("MER111_", cl), k, l, or, p))
    # end

    # for (cl, KMK) in kd, k = 1:KMK.k, n in sort(unique(mer11.name))
    #     p, or = ProximityEnrichment.hypertest(mcl[!, cl] .== k, mer11[!, :name] .== n)
    #     push!(clusterassociations, (string("MERFAM_", cl), k, n, or, p))
    # end
    


    (mcl=mcl, MMS=MMS, ZMK27=ZMK27, ZMK9=ZMK9, ZMS=ZMS, proxmcl=proxmcl, proxsum=proxsum, kd=kd, silhouetteplot=silhouetteplot, clusterassociations=clusterassociations)

end



kernelsmooth(x::Vector{T}, σ, pad=NA()) where {T} = imfilter(x, KernelFactors.gaussian(σ, Int(ceil(4σ)) + 1), pad)
smoothpile(p, dx, σ=2*dx) = average_heatmap(kernelsmooth(p, σ), dx)
#### plot epigentic data 

function plot_gene_snapshot(chrom, loc, peakmeta, FM, mer11, genemodels; dx = 150, merex=500, label = "epi_all", transcripts=String[], ns=15, plottype=:grid, yld=Dict{String, Float64}())


    if isempty(transcripts) ## if no transcripts defined, select longest transcript per gene
        genecoords = @chain genemodels begin
            @subset(:chrom .== chrom, (:start .∈ Ref(loc)) .| (:stop .∈ Ref(loc)))
            groupby(:GeneName)
            combine(df -> df[argmin(df.Rank), :])
        end
    else
        genecoords = @subset(genemodels, :TranscriptID .∈ Ref(Set(transcripts)))
        
        if isempty(genecoords)
            
            genecoords = @subset(genemodels, :GeneName .∈ Ref(Set(transcripts)))
            genecoords = combine(groupby(genecoords, :GeneName), df -> df[argmin(df.Rank), :])
        end
        
    end

    PS = [smoothpile(P, dx) for P in piles(chrom, loc, peakmeta.Index, label, peakmeta.Label, FM)]
    peakindex = Dict(peakmeta.Index[i] => i for i = 1:size(peakmeta, 1))
    peakmax = maximum.(PS)

    peakgroup = combine(groupby(peakmeta, [:Chip, :Run, :GenoType]), nrow => :Count, :Index => Ref => :Inds)
    
    maxplots = maximum(peakgroup.Count)
    DataFrames.transform!(groupby(peakgroup, [:Run, :Chip]), :Inds => (x -> maximum(maximum.(PS[getindex.(Ref(peakindex), vcat(x..., ))]))) => :PeakMax)

    # showwide(peakgroup)

    mersel = @subset(mer11, :chrom .== chrom, (:start .∈ Ref(loc)) .| (:stop .∈ Ref(loc)))

    stages = unique(peakmeta.Stage)
   
    ### colour by stages
    ##cc = cgrad(:viridis, length(stages) + 2, categorical=true)[2:end-1]
    ##cdict = Dict(s => c for (s, c) in zip(stages, cc))

    cdict = Dict(("WT", "K9") => :steelblue1, ("WT", "K27") => :tomato, ("KO", "K9") => :steelblue, ("KO", "K27") => :orange, ("iPSC", "K9") => :darkblue, ("iPSC", "K27") => :darkorange2)


    lxp = first(loc):dx:last(loc)

    aphs = Vector{Plots.Plot}[]

    pblank = plot(framestyle=:none)

    pgenemodel = plot(ylims=(0, 2), yticks=false, yaxis=false, xlabel=string(chrom, ":", first(loc), "-", last(loc)), xguidefont=font(8, "helvetica"), xlims=(first(loc), last(loc)))
    yo = 1
    for (t, g) in zip(genecoords.TranscriptID, genecoords.GeneName)
        
        plotgenemodel!(t, genecoords, label=g, pos=:right, bottom_margin=5mm, ns=ns, y=yo)
        # yo += 1
    end
    plot!(ylims=(0, yo+1))
    annotate!([((s + e)/2, 1.5, text("MER11", font("helvetica", 10))) for (s, e) in zip(mersel.start, mersel.stop)])


    for row in eachrow(peakgroup)
        phs = [plot(bottom_margin=-1mm, framestyle=:box) for i = 1:length(row.Inds)]
        for (p, ri) in zip(phs, row.Inds)
            i = peakindex[ri]
            pm = max(1.1*row.PeakMax, 0.5)
            pm = get(yld, row.Chip, pm)
            
            plot!(p, ylims=(0, pm))
            for (s, e) in zip(mersel.start, mersel.stop)
                se = s - merex
                ee = e + merex
                plot!(p, rectangle(ee - se, pm, se, 0), c=:black, fillalpha=0.15, line=stroke(0))
            end
            #plot!(p, lxp, PS[i], lab="", fill=0, c=cdict[peakmeta.Stage[i]], yticks=ifelse(row.PeakMax > 2.1, [0, 2], [0, 1]))
            plot!(p, lxp, PS[i], lab="", fill=0, c=cdict[(row.GenoType, row.Chip)], yticks=ifelse(pm >= 2.1, [0, 2], [0, 1]))
            annotate!(p, [(first(loc) + 1000, 1, text(string(ifelse(peakmeta.Chip[i] .== "K9", "H3K9me3", "H3K27ac"), " ", peakmeta.GenoType[i]), font("helvetica", 9, :left)))])
        end

        for k in (length(phs) + 1):maxplots
            push!(phs, pblank)
        end

        (plottype == :grid) && push!(phs, pgenemodel)

        push!(aphs, phs)
    end

    pphs = mapreduce(p -> vcat(p..., ), vcat, zip(aphs...))
    if plottype == :grid
        for p in aphs[1]
            plot!(p, right_margin=-1mm)
        end

        for phs in aphs[2:end]
            for p in phs
                plot!(p, left_margin=-1mm, right_margin=-1mm)
            end
        end

        
        plot(pphs..., layout=(maxplots + 1, size(peakgroup, 1)), size=(1200, 600), xlims=(first(loc), last(loc)), leg=false, xticks=false, top_margin=-2mm, grid=false)
    else
        push!(pphs, pgenemodel)
        plot(pphs..., layout=(length(pphs), 1), size=(1200, 600), xlims=(first(loc), last(loc)), leg=false, xticks=false, top_margin=-2mm, grid=false)
    end

end


function annotatemerdist!(deg, rois, label)
    gi = @with deg Interval.(:chrom, :TSS, :TSS)
    c = ProximityEnrichment.closest_tss_peak(gi, rois) 
    deg[!, label] = ifelse.(c .== -1, missing, c)    
      
    
    deg
end
