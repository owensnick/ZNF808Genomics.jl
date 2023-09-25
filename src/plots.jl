
### savedisplay save svg and display current figure
function savedisplay(label, fmts=[".svg", ".png"], projdir=getprojectdir(), verbose=false)
    # file = joinpath(projdir, "figures", string(label, fmt))
    files = joinpath.(projdir, "figures", string.(label, fmts))
    verbose && println("Writing ", first(files))
    mkpath(dirname(first(files)))
    map(savefig, files)
    
    p = plot!()
    display(p)
end


### humanreadable axis ticks
function humanread(x, s)
    
    c = floor(log10(abs(x)))
    
    if c < 3
        return s
    else
        n = replace(Humanize.datasize(x), r"\.0| |B$" => "")
        n
    end
end


function plottotalde(deg, stagelabels ; kwargs...)
    stages = unique(deg.Stage)
    
    dtf = combine(groupby(deg, :Stage), nrow => :count, :sig => sum => :sig, :dsig => (x -> sum(x .== 1)) => :Up, :dsig => (x -> sum(x .== -1)) => :Down)
    groupedbar(dtf.Stage, [dtf.Up dtf.Down]/1000, bar_position=:stack, lab=["Activated" "Repressed"], c=[:orange :steelblue], fontfamily="helvetica", framestyle=:box, grid=false, line=stroke(:white), ylabel="Total Genes (K)", xlabel="Stage", title="Dysregulated Genes (DGs) in ZNF808 KO" ; kwargs...)
    annotate!([(i - 0.5, (dtf.Down[i] + dtf.Up[i])/1000, text(dtf.Down[i],   font(10, :steelblue,  :bottom, "helvetica"))) for i = 1:length(stages)])
    annotate!([(i - 0.5, (dtf.Down[i] + dtf.Up[i])/1000 + 0.0, text(string(dtf.Up[i], "\n"), font(10, :orange, :bottom, "helvetica"))) for i = 1:length(stages)])
    ylims!(0, ylims()[2] + 0.5)
    plot!(xticks=((1:length(stages)) .- 0.5, getindex.(Ref(stagelabels), stages)))
end


function plotpca(design, counts; group=design.Stage, pcc = (1, 2), thr=10, kwargs...)

    M = Matrix(counts[!, design.Sample])
    ind = vec(maximum(M, dims=2)) .> thr 
    LX = log10.(M[ind, :] .+ 1)
    @show extrema(LX)
    pcaD = fit(PCA, LX)
    TD = MultivariateStats.transform(pcaD, LX)'
    explained = principalvars(pcaD)/tvar(pcaD)
    pcastats = DataFrame([1:length(explained) explained cumsum(explained)], [:PC, :VarianceExplained, :CumlativeVarianceExplained])
#     first(pcastats, 3) |> showwide
    
    scatter(TD[:, pcc[1]], TD[:, pcc[2]], group=group, marker=(:auto), framestyle=:box, grid=false, leg=:outerright; kwargs...)
    # annotate!([(TD[i, pcc[1]], TD[i, pcc[2]], text(design.Sample[i], font(:left, 7))) for i = 1:size(LX, 2)])
    xlabel!(string("PC", pcc[1], " : ", round(explained[pcc[1]], digits=2)))
    ylabel!(string("PC", pcc[2], " : ", round(explained[pcc[2]], digits=2)))

end


### clustering
function clusterlines!(KMK)
    c = KMK.KM.counts[KMK.KSI]
    cs = cumsum(c) 
    hline!(cs, lab="", c=:white)
    yticks!(cs .- c/2, string.(1:KMK.k))
end


function bubbleplot(C, S=C; mso=5, msm=2, sp=0.5, xoff=0, kwargs...)
    @assert size(C) == size(S)
    XY = [(j, i) for i = 1:size(C, 1), j=1:size(C, 2)][:]    
    scatter(first.(XY) .+ xoff, last.(XY), ms=msm.*sqrt.(S[:]) .+ mso, zcolor=C[:], yflip=true, grid=false, xlims=(1 - sp, size(C, 2) + sp) .+ xoff, ylims=(1 - sp, size(C, 1)+sp), marker=stroke(0); kwargs...)
end


################################# proximity enrichment


proxplot(proxdf, roi; data=:pvalue, kwargs...) = proxplot(proxdf, [roi], data=data; vl=[], kwargs...)
function proxplot(proxdf, rois::Vector{T} ; data=:pvalue, cs = [:steelblue :orange], group=:Direction, labs=[""], titles=String[], yts=0:5:25, vl=[], dx=1, kwargs...) where T
    # @show rois
    if data == :pvalue
        field = :pvalue
        trf = x -> -log10(x)
        label = "-log10 p Fisher"
    elseif data == :count
        field = :count
        trf = x -> x
        label = "Count"
    elseif data == :or
        field = :or
        trf = x -> x
        label = "Odds Ratio"
    else
        error("Plotting data = $(data) unknown")
    end

    stages = unique(proxdf.Stage)
    if isempty(titles)
        titles = [string(join(rois, ", "), ":", s) for s in stages]
    end
    phs = [@with @subset(proxdf, :ROI .∈ Ref(Set(rois)), :Stage .== s) plot(:x[1:dx:end], trf.(cols(field))[1:dx:end], group=cols(group)[1:dx:end], fill=0, fillalpha=0.2, title=t, c=cs, lab=labs) for (s, t) in zip(stages, titles)]
    
    if !isempty(vl)
        for p in phs
            vline!(p, vl, c=:black, ls=:dash, lab="")
        end
    end

    for p in phs[2:end]
        plot!(p, yticks=false, left_margin=-2mm)
    end
    plot!(phs[1], yticks=yts, ylabel=label, left_margin=5mm)
    for p in phs[1:end-1]
        plot!(p, leg=false)
    end
    ylm = max(maximum(ylims(p)[2] for p in phs), 15)

    plot(phs..., layout=(1, length(phs)), link=:y, ylims=(0, ylm), size=(1000, 200), xlims=(2.8, 8), grid=false, xticks=(3:2:7, ["1kb", "100kb", "10Mb"]), top_margin=5mm, titlefont=font(11, "helvetica"), fontfamily="helvetica" ; kwargs...)

    
end


function epiproxproption(deresults)
    dwd = @subset(unstack(combine(groupby(deresults.deg, [:Stage, :dsig, :pMER11]), nrow => :count), [:Stage, :dsig], :pMER11, :count), abs.(:dsig) .> 0)
    rename!(dwd, [:Stage, :dsig, :Out, :In])
    dwd.PR = dwd.In./(dwd.Out .+ dwd.In)
    
    plot(dwd.Stage, dwd.PR, group=dwd.dsig, marker=(:circle, stroke(:white)), c=[:steelblue :orange], lab=["Repressed" "Activated"], ylabel="Proportion")
    plot!(xticks=(((1:length(deresults.stages)) .- 0.5), getindex.(Ref(stagelabels), deresults.stages)))
    hline!([mean(unique(deresults.deg[!, [:Gene, :pMER11]]).pMER11)], c=:black, lab="Geneome wide background", ls=:dash)
    
    plot!(leg=:topright, grid=false, xlabel="Stage", title="Proportion DGs ≤ 1MB to\nMER11 H3K9me3 → H3K27ac", fontfamily="helvetica", titlefont=font("helvetica", 12))
end



function offsets(x)
    if length(x) == 2
        return [-.2, .2]
    elseif length(x) == 3
        return [-.25, 0, 0.25]
    elseif length(x) == 4
        return range(-.25, .25, length=4)
    else
        return 0.1.*randn(length(x))
    end
end

sw(x) = x/sum(x)

function plotgenesimple(gene, deresults, stagelabels; meta = deresults.design, ncounts=deresults.ncounts, kwargs...)
    meta.GenoTypeIPS = ifelse.(meta.Clone .== "Patient", "iPSC", meta.GenoType)

    plotgene(gene, meta, ncounts, size=(600, 300), groupf=[:GenoTypeIPS], leg=:outertopright, annotlabel=false, plotpair=false, ps=[:steelblue, :orange, :darkgrey]; kwargs...)# trans = x -> log10.(x .+ 1))
    
    gn = replace(string(gene), r"[r\"\|\\\$\^]" => "")
    plot!(title=gn, ylabel="Normalised Reads", grid=false)
    
    stages = unique(meta.Stage)
    plot!(xticks=(1:length(stages), getindex.(Ref(stagelabels), stages)))
    vline!((1:(length(stages) -1)) .+ 0.5, c=:grey, lab="")
    
    p = plot!(framestyle=:box, xlims=(0.5, length(stages) + .5), fontfamily="helvetica")
end


function plotgene(gene, meta, tpm; kwargs...)
    ind = findall(occursin.(gene, tpm.Gene))
    if isempty(ind)
        error("$ind not found")
    elseif length(ind) > 1
        error("Multiple matches for $ind: $(tpm.Gene[ind])")
    else
        plotgene(first(ind), meta, tpm; kwargs...)
    end
end


function plotgene(i::Int, meta, tpm; fmt=:short, label="", annot=:sc, trans=identity, groupf=[:Clone], plotpair=true, annotlabel=true, ps = PlotThemes.wong_palette, markips=false, kwargs...)
    y = [trans(tpm[i, l]) for l in meta.Sample]
    ynt = [tpm[i, l] for l in meta.Sample]
    if fmt == :long
        si = sortperm(meta, [:Stage, :GenoType, :Clone])
        scatter(meta.Index, y[si], xticks=(meta.Index, replace.(meta.Sample[si], Ref("_" => "\n"))), c=ifelse.(meta.GenoType[si] .== "WT", :blue, :red); title=tpm.Gene[i], kwargs...)
    else
            
        
    if annot == :sc ## only for genotype
        gmeta = copy(meta)
        gmeta[!, :y] = y;
        gmeta = DataFrames.transform(groupby(gmeta, [:Stage ; groupf]), :y => (x -> offsets(x)/4) => :ROff)

        gg = combine(groupby(gmeta, [:Stage ; groupf]), nrow => :count, :y => Ref => :y, :Sample => Ref => :Sample, :ExpNum => Ref => :ExpNum, :ROff => Ref => :ROff, :Clone  => Ref => :Clones)

      
        ust = unique(gg.Stage)
        ugt = unique(gg[!, first(groupf)])
        xp = Dict(s => i for (i, s) in enumerate(ust))
            
            pp = Dict(gt => p for (gt, p) in zip(ugt, ps[1:length(ugt)]))
            gto = Dict(gt => o for (gt, o) in zip(ugt, offsets(ugt)))
            
            labels = Set{String}()
            p = plot(; kwargs...)        
            
            fgl = first(groupf)
            if plotpair
                # display(gmeta)
                for row in eachrow(combine(groupby(gmeta, [:Stage, :ExpNum]), :Clone => Ref => :Clone, :y => Ref => :y, :Sample => Ref => :Sample, :ROff => Ref => :ROff, first(groupf) => Ref => :GroupF))
                    # showwide(row)
                    # display(gto)
                    # display(ugt)
                    x = xp[row.Stage] .+ getindex.(Ref(gto), row.GroupF) .+ row.ROff
                    indh1 = row.Clone .== "H1"
                    fi = findall(.!indh1)
                    for f in fi
                        plot!([x[indh1] ; x[f]], [row.y[indh1] ; row.y[f]], c=ifelse(row.ExpNum .== "E11", :red, :grey), lab="")
                    end
                end
            end
            for row in eachrow(gg)
                #roff = offsets(row.y)/4
                x = ones(length(row.y)).*xp[row.Stage] .+ gto[row[fgl]] .+ row.ROff # roff
                
                if markips && ("Patient" ∈ row.clones)
                    # ipsind = row.Clones .== "Patient"
                    # scatter!(x[.!ipsind], row.y[.!ipsind], c=pp[row[fgl]], lab=ifelse(row[fgl] ∈ labels, "", row[fgl]), marker=ifelse.(row.Clones .== "Patient", :utriangle, :circle))
                else
                    scatter!(x, row.y, c=pp[row[fgl]], lab=ifelse(row[fgl] ∈ labels, "", row[fgl]), marker=ifelse.(row.Clones .== "Patient", :utriangle, :circle))
                end
                # scatter!(x, row.y, group=row.ExpNum, lab=ifelse(row.Clone ∈ labels, "", row.ExpNum |> permutedims), c = [1 2 3 4])
                
                
                if annotlabel
                    lr = ifelse(row.Clone == "H1", :right, :left)
                    z = ifelse(row.Clone == "H1", -.2, .2)
                    annotate!([(x[i] + z, row.y[i], text(row.Sample[i]*row.ExpNum[i], font(lr, 7))) for i = 1:length(row.y)])
                end
                push!(labels, row[fgl])
            end
            
                  
            plot!(xticks=(1:length(ust), ust), xlims=(0.1, length(ust) .+ .5), ylims=(0, ylims()[2]), title=string(label, "\n", tpm.Gene[i]))
        
        elseif  annot == :plots
            phs = Plots.Plot[]
            gs = last(split(tpm.Gene[i], "|"))
            for smeta in groupby(meta, :Stage)
                yg = [tpm[i, l] for l in smeta.Sample]
                p = groupeddotplot(smeta.GenoType, yg, group=(map(f -> smeta[!, f], groupf)...,), marker=:auto, title=string(first(smeta.Stage), " : ", gs), leg=false, xlims=(0, 2.5))
                yl = ylims()
                plot!(p, ylims=(0, yl[2]))
                push!(phs, p)
            end
            smeta = @subset(meta, :Stage .== "D0")
            yg = [tpm[i, l] for l in smeta.Sample]
            p = groupeddotplot(smeta.GenoType, yg, group=(map(f -> smeta[!, f], groupf)...,), marker=:auto, title="", framestyle=:none, xlims=(-20, -19), leg=:outerleft)
            push!(phs, p)
            
            
            
            plot(phs..., layout=grid(1, length(phs), widths=sw([ones(length(phs)-1) ; 2])), titlefont=font(10), ytickfont=font(6), left_margin=-1mm, right_margin=-1mm; kwargs...)
            
        else
            groupeddotplot(meta.Stage, y, group=(meta.GenoType, meta.Run, meta.Clone), marker=:auto, title=string(label, "\n", tpm.Gene[i]); kwargs...)
    
        end
    end
end


function geneheatmap_scatter(pattern, deresults; meta=deresults.design, ncounts=deresults.ncounts, k = 0, kwargs...)

    sel = @subset(deresults.meancounts, occursin.(pattern, :Gene))
    sort!(sel, :Gene)
    # deg = @subset(deresults.deg, occursin.(pattern, :Gene))
    # fcdeg = coalesce.(unstack(deg, :Gene, :Stage, :log2FoldChange), 0)
    # sgdeg = coalesce.(unstack(deg, :Gene, :Stage, :sig), 0)
    
    # sort!(fcdeg, :Gene)
    # sort!(sgdeg, :Gene)
    # @assert sel.Gene == fcdeg.Gene
    # @assert sel.Gene == sgdeg.Gene
    
    stages = deresults.stages
    
    
    nwt = string.("WT_", stages)
    MW = Matrix(sel[!, nwt])
    m = maximum(MW, dims=2)
    ZW = MW./m

    if k > 1
        KMK = kmeansorder(ZW', k)
        si = reverse(KMK.SI)
    else

        si = sortperm(DataFrame(ZW, :auto))
    end
    
    gns = last.(split.(sel.Gene, "|"))
    
    stl = getindex.(Ref(stagelabels), stages)
    
    pw = heatmap(ZW[si, :], xticks=(1:size(ZW, 2), stl), yticks=(1:size(ZW, 1), gns[si]), title="WT", clims=(0, 1))

    phs = [plotgenesimple(Regex(string("\\|", mtf, "\$")), deresults, stagelabels, meta = meta, ncounts=ncounts, trans=x->x,  leg=:topright, ylabel="") for mtf in gns[reverse(si)]]

    for (tf, p) in zip(gns[reverse(si)], phs)
        plot!(p, ylabel="", left_margin=-2mm, title="", bottom_margin=-1mm)
        yt = first(yticks(p))
        ytn = first(yt)[[1, end]]
        yts = humanread.(ytn[[1, end]], last(yt)[[1, end]])
        yticks!(p, ytn, yts)
        annotate!((0.75, 0.1*(ytn[2] - ytn[1]), text(tf, font("helvetica", 10, :left, :bottom))))
    end

    nr = cld(length(phs), 3)

    for i = (length(phs)+1):(nr*3)
        push!(phs, plot(framestyle=:none))
    end
    lt = @layout [a{0.35w} grid(nr, 3)]
    plot(pw, phs..., layout=lt, size=(1000, 400); kwargs...)

end



#### human embryo rnaseq

## bar chart version
function plot_human_embryo_rnaseq(hem, pattern=r"ZNF808|ZNF440|ZNF525|ZNF468|ZNF433")
    hem_znf = @chain hem begin
        @subset(occursin.(pattern, :gene_name))
        stack(names(hem, r"_[0-9]$"), :gene_name, variable_name=:tissue, value_name=:count)
        DataFrames.transform(:tissue => (x -> replace.(replace.(x, r"_[0-9]*$" => ""), "_" => "\n")) => :Group)
        groupby([:gene_name, :Group])
        combine(:count => mean => :count)
    end
    
    hem_mean = combine(groupby(hem_znf, :Group), :count => mean => :meancount)
    hem_mean.Order = sortperm(sortperm(hem_mean.meancount))
    sort!(hem_mean, :Order)
    hem_znf = innerjoin(hem_znf, hem_mean, on=:Group)
    
    n = length(unique(hem_znf.gene_name))

    cc = cgrad(:viridis, max(n, 2), categorical=true)[1:end]
    @with hem_znf groupedbar(:Order, :count, group=:gene_name, size=(1000, 300), margin=5mm, leg=^(:topleft), ylabel="Normalised Read Counts", title="Human Embryo RNA-seq", xlabel="Tissue", c=permutedims(cc))
    xticks!(hem_mean.Order, hem_mean.Group)
end

### boxplot version
function plot_human_embryo_rnaseq_boxplot(human_embryo_rnaseq, pattern=r"ZNF808|ZNF440|ZNF525|ZNF468|ZNF433|ZNF578")

        ### table of MER11 binding ZNFs
    hem_znf = @chain human_embryo_rnaseq begin
        @subset(occursin.(pattern, :gene_name))
        stack(names(human_embryo_rnaseq, r"_[0-9]$"), :gene_name, variable_name=:tissue, value_name=:count)
        DataFrames.transform(:tissue => (x -> replace.(replace.(x, r"_[0-9]*$" => ""), "_" => " ")) => :Group)
    end

    ### order tissues by expression of ZNF808
    hem_znf808 = @subset(hem_znf, :gene_name .== "ZNF808")
    gdf = combine(groupby(hem_znf808, :Group), :count => mean => :count)
    sort!(gdf, :count)
    gdf.GI = 1:size(gdf, 1)
    hem_znf808 = leftjoin(hem_znf808, gdf[!, [:Group, :GI]], on=:Group)

    ### plot ZNF808
    cc = cgrad(:viridis, categorical=true, size(gdf, 1))[1:end] |> permutedims
    pz = plot()
    @with hem_znf808 boxplot!(:GI, :count, group=:GI, title=:gene_name[1], right_margin=-2mm, c=cc, outliers=false)
    @with hem_znf808 dotplot!(:GI, :count, group=:GI, title=:gene_name[1], right_margin=-2mm, c=cc, marker=(3, ^(:black)))
    plot!(xticks=(gdf.GI, gdf.Group), leg=false, xrotation=45, ylims=(0, 450), grid=false, ylabel="Normalised Read Counts", left_margin=5mm)

    ### Plot other ZNFs
    phs = Plots.Plot[]
    for df in groupby(@subset(hem_znf, :gene_name .!= "ZNF808"), :gene_name)
        
        df = leftjoin(df, gdf[!, [:Group, :GI]], on=:Group)
        p = plot()
        cc = cgrad(:viridis, categorical=true, size(gdf, 1))[1:end] |> permutedims
        @with df boxplot!(:GI, :count, group=:GI, title=:gene_name[1], right_margin=-2mm, c=cc, outliers=false)
        @with df dotplot!(:GI, :count, group=:GI, title=:gene_name[1], right_margin=-2mm, c=cc, marker=(3, ^(:black)))
        plot!(xticks=(gdf.GI, gdf.Group),  grid=false,xtickfont=font("helvetica", 7))
        push!(phs, p)
    end
    ### tidy axis
    for p in phs[[2, 3, 5]]
    plot!(p, yformatter=y -> "", left_margin=-2mm) 
    end

    lt = @layout [a{0.5w} grid(2, 3)]

plot(pz, phs..., plot(yticks=false, axis=false, xticks=false), layout=lt, size=(1000, 500), link=:y, leg=false, xrotation=45, bottom_margin=15mm, fontfamily="helvetica",)
end




function horzenrichmentplot(enfisher, len; geneset = "Human_Gene_Atlas", upterms = ["Fetalliver", "Liver"], downterms=["Pancreas", "PancreaticIslet"], thr=25)
    
    ent_up = @subset(enfisher, :Direction .== "Up", :GeneSet .== geneset, :Term .∈ Ref(Set(upterms)))[!, [:Stage, :Term, :OR, :FDR]]
    ent_dn = @subset(enfisher, :Direction .== "Down", :GeneSet .== geneset, :Term .∈ Ref(Set(downterms)))[!, [:Stage, :Term, :OR, :FDR]]
    ent = [ent_up ; ent_dn]
    rename!(ent, [:Stage, :Label, :OR, :FDR])
    len = @subset(deresults.len, :Label .∈ Ref(Set(["UD", "DU"])))[!, [:Stage, :Label, :OR, :FDR]]
    ort = unstack([len ; ent], :Stage, :Label, :OR)
    fdrt = unstack([len ; ent], :Stage, :Label, :FDR)
    
    dirvec = hcat(ones(length(upterms))', -ones(length(downterms))', 1, -1)

    FDR = -log10.(coalesce.(Matrix(fdrt[!, [upterms ; downterms ; ["UD", "DU"]]]), 1.0))
    
    FDR = min.(FDR, thr).*dirvec
    
    OR  = coalesce.(Matrix(ort[!, [upterms ; downterms ; ["UD", "DU"]]]), 1.0)
    OR = max.(0.0, log2.(OR))
    
    bubbleplot(FDR', abs.(FDR'), mso=1, msm=4,xoff=-0.5, c=cgrad([:steelblue, :white, :orange]), colorbar=false, clims=(-thr, thr), xticks=((1:length(fdrt.Stage)) .- 0.5, getindex.(Ref(stagelabels),fdrt.Stage)))

    ### annotate 
    annotate!([(i .- 0.5, j, text(string(round(abs(FDR[i, j]), digits=1)), font("helvetica", 7, :center, :middle))) for i = 1:size(FDR, 1), j = 1:size(FDR, 2)][:])


    plot!(size=(350, 300), leg=false, yticks=(1:(length(upterms) + length(downterms) + 2), [upterms ; replace.(downterms, "PancreaticIslet" => "Pancreatic\nIslet") ; ["Hepatic\nCoords", "Dorsal\nPancreas"]]), fontfamily="helvetica")
end

function panc_hep_enrichment_plot(deresults, laser_dorsal_hep)

    pancgenes  = unique(@subset(innerjoin(deresults.deg, laser_dorsal_hep[!, [:GeneName, :LSD]], on=:GeneName), :LSD .== 1,  :dsig .!= 0, :Stage .∈ Ref(Set(["S2", "S3", "S4"]))).GeneName)
    livergenes = unique(@subset(innerjoin(deresults.deg, laser_dorsal_hep[!, [:GeneName, :LSD]], on=:GeneName), :LSD .== -1, :dsig .!= 0, :Stage .∈ Ref(Set(["S2", "S3", "S4"]))).GeneName);

    pfc = plot(fmt=:png)
    @with @subset(deresults.deg, :GeneName .∈ Ref(Set(pancgenes)))  violin!(:Stage, :log2FoldChange, lab="Dorsal Pancreas", side=^(:left), c=^(:steelblue), fillalpha=0.2, bandwidth=0.09)
    @with @subset(deresults.deg, :GeneName .∈ Ref(Set(livergenes))) violin!(:Stage, :log2FoldChange, lab="Hepatic Cords", side=^(:right), c=^(:orange), fillalpha=0.2, bandwidth=0.09)
    @with @subset(deresults.deg, :GeneName .∈ Ref(Set(pancgenes)))  dotplot!(:Stage, :log2FoldChange, c=^(:steelblue), marker=(2, 0.45, stroke(^(:white))), side=^(:left), lab="")
    @with @subset(deresults.deg, :GeneName .∈ Ref(Set(livergenes))) dotplot!(:Stage, :log2FoldChange, c=^(:orange), marker=(2, 0.45, stroke(^(:white))), side=^(:right), lab="")
    plot!(ylims=(-2.5, 2.5), ylabel="log2 Fold Change WT/KO", leg=:bottomright, size=(400, 315), xlabel="Stage", fontfamily="helvetica")

    plot!(xticks=((1:length(deresults.stages)) .- 0.5, getindex.(Ref(stagelabels), deresults.stages)), grid=false)
    # p = plot!(title="Relative change of Dorsal pancreas and Hepatic cord genes over differentiation")
    hline!([0], lab="", c=:black)
    


    ph = horzenrichmentplot(deresults.enfisher, deresults.len, upterms=[], downterms=[],thr=25)
    plot!(yaxis=false)

    plot(ph, pfc, layout=grid(2, 1, heights=[0.3, 0.7]), size=(300, 350), link=:x, xlims=(0, 5))
    
end

function scattline(x, y; label="", kwargs...)

    p = histogram2d(x, y, bins=200, ; lab="", smooth=true, colorbar=false, kwargs...) 
    m = lm(@formula(Y ~ X), DataFrame(X=x, Y=y))
    ct = coeftable(m)
    pv = ct.cols[ct.pvalcol][2]
    c = coef(m)
    xp = [-10, 10]
    yp = predict(m, DataFrame(X=xp))
    plot!(xp, yp, lab="", c=:black)
    annotate!((-10, 10, text(string("p ", StatsBase.PValue(pv)), font(8, "helvetica", :top, :left))))
    annotate!((-10, 8, text(string("p ", pv), font(8, "helvetica", :top, :left))))
    annotate!((-10, -10, text(string(round(c[2], digits=2), "x + ", round(c[1], digits=2)), font(8, "helvetica", :bottom, :left))))
    annotate!((10, -10, text(label, font(14, :bottom, :right))))
    p
end


function epi_cluster_proximity_enrichments(epiclusters, deresults)
    proxsum = @chain epiclusters.proxmcl begin
        @subset(:Direction .== "Up", occursin.("All", :ROI), .!occursin.(r"_(1|2|123|34)$", :ROI))
        groupby([:ROI, :Stage])
        combine(df -> df[argmin(df.pvalue), :])
        sort!([:Direction, order(:ROI, by = x -> ifelse(x == "MER11", "Z", x))])
    end
    codedict = Dict("12" => 1, "3" => 2, "4" => 3, "5" => 4)
    proxsum.Cluster = getindex.(Ref(codedict), last.(split.(proxsum.ROI, "_")))
    sort!(proxsum, :Cluster)


    fdrt = coalesce.(unstack(proxsum, [:Cluster, :Direction], :Stage, :pvalue), 1.0)
    ort  = coalesce.(unstack(proxsum, [:Cluster, :Direction], :Stage, :or), 1.0)

    F = min.(-log10.(Matrix(fdrt[!, deresults.stages])), 10)
    O = max.((Matrix(ort[!, deresults.stages])), 0)

    F = [F ones(size(F, 1))]
    O = [O range(40, 0, length=5)[1:end-1]] ## Add key

    bubbleplot(F, O, mso=0, msm=4, xticks=(1:size(F, 2), [deresults.stages ; "Odds Ratio\nKey"]), yticks=(1:size(F, 1), string.(fdrt.Cluster)), lab="", clims=(0, 10), fontfamily="helvetica", marker=(stroke(1)))
    annotate!([(size(F, 2), i, text(Int.(O[i, end]), font(6, :white, "helvetica"))) for i = 1:size(O, 1)])
    plot!(title="Maximal Proximity Enrichments\nMER11 clusters ↔ dysregulated genes", size=(600, 325), ylabel="Cluster")
end

firstmember(x) = x
firstmember(x::Tuple{T, V}) where {T, V} = first(x)

function enrichrplot(enfisher, tgs_up=[], tgs_dn=[], thr=10)
    
    upind = falses(size(enfisher, 1))
    dnind = falses(size(enfisher, 1))
    if !isempty(tgs_up)
        upind = enfisher.Direction .== "Up"
        upind .&= mapreduce(i -> termind(enfisher, i), +, tgs_up) .> 0
    end
    if !isempty(tgs_dn)
        dnind = enfisher.Direction .== "Down"
        dnind .&= mapreduce(i -> termind(enfisher, i), +, tgs_dn) .> 0
    end
    
    upterms = string.("Activated\n", firstmember.(tgs_up))
    dnterms = string.("Repressed\n", firstmember.(tgs_dn))
    
    ent_up = enfisher[upind, :]
    ent_dn = enfisher[dnind, :]
    
    # showwide(ent_up[!, [:GeneSet, :Term]] |> unique)
    # showwide(ent_dn[!, [:GeneSet, :Term]] |> unique)
    
    ent_up.Term = string.("Activated\n", ent_up.Term);
    ent_dn.Term = string.("Repressed\n", ent_dn.Term);
    ent = [ent_dn ; ent_up]
#     showwide(ent_up)
    ort  = unstack(ent, :Stage, :Term, :OR)
    fdrt = unstack(ent, :Stage, :Term, :FDR)
    dirvec = hcat(-ones(length(dnterms))', ones(length(upterms))')

    terms = [dnterms ; upterms]

    FDR = min.(-log10.(coalesce.(Matrix(fdrt[!, terms]), 1.0)), thr).*dirvec
    OR  = max.(coalesce.(Matrix(ort[!, terms]), 1.0), 1.0)
    
    # @show extrema(OR)
    
    FDR = [FDR fill(0, size(FDR, 1))]
    # OR = [OR [NaN ; NaN ; range(10, 0, length=size(OR, 1))[1:end-1] ;  NaN]]
    OR = [OR range(10, 2, length=size(OR, 1))]
    

    
    terms = replace.(terms, "PancreaticIslet" => "Pancreatic\nIslet")
    terms = replace.(terms, "regulation of transcription from RNA polymerase II promoter (GO:0006357)" => "Reg.\nTranscription")
    terms = replace.(terms, "TGF-beta regulation of extracellular matrix" => "Regulation\nTGF-beta")
    terms = replace.(terms, "negative regulation of cell differentiation (GO:0045596)" => "Neg. Reg.\nDifferentiation")
    terms = replace.(terms, r"Activated\n|Repressed\n" => "")
    bubbleplot(FDR, abs.(OR), mso=1, msm=7, c=cgrad([:steelblue, :white, :orange]), clims=(-thr, thr), yticks=(1:length(fdrt.Stage), getindex.(Ref(stagelabels),fdrt.Stage)), marker=stroke(1))
    # annotate!([(length(terms) + 1, i, text(OR[i, end], font(6, :black, "helvetica"))) for i = 3:(size(OR, 1)-1)])
    annotate!([(length(terms) + 1, i, text(Int(OR[i, end]), font(6, :black, "helvetica"))) for i = 1:size(OR, 1)])
    plot!(leg=false, xticks=(1:(length(terms) + 1), [terms ; "Odds Ratio\nKey"]), fontfamily="helvetica", ylabel="Stage", title="Enrichr Gene Set Enrichments", xrotation=30)
end

tt(x) = string(x)
tt(x::Float64) = string(round(x, digits=2))
function epi_bubble(epiclusters, deresults, clustercode="MI_CL6"; msm=4, npmax=15, omax=80, plottotals=false)
    proxsum = @chain epiclusters.proxmcl begin
        @subset(:Direction .== "Up", occursin.(clustercode, :ROI))
        groupby([:ROI, :Stage])
        combine(df -> df[argmin(df.pvalue), :])
        sort!([:Direction, order(:ROI, by = x -> ifelse(x == "MER11", "Z", x))])
    end
    proxsum.Cluster = last.(split.(proxsum.ROI, "_"))

    sort!(proxsum, :Cluster)

    fdrt = coalesce.(unstack(proxsum, [:Cluster, :Direction], :Stage, :pvalue), 1.0)
    ort  = coalesce.(unstack(proxsum, [:Cluster, :Direction], :Stage, :or), 1.0)
    ct   = coalesce.(unstack(proxsum, [:Cluster, :Direction], :Stage, :count), 0)
    
    C = Matrix(ct[!, deresults.stages])

    F = min.(-log10.(Matrix(fdrt[!, deresults.stages])), npmax)
    O = min.(max.((Matrix(ort[!, deresults.stages])), 0), omax)

    F = [F ones(size(F, 1))]
    @show extrema(O)
    O = [O range(omax, 0, length=size(O, 1))] ## Add key
    O[end, end] = max(O[end, end], 1)

    bubbleplot(F, O, mso=0, msm=msm, xticks=(1:size(F, 2), [deresults.stages ; "Odds Ratio\nKey"]), yticks=(1:size(F, 1), string.(fdrt.Cluster)), lab="", clims=(0, npmax), fontfamily="helvetica", marker=(stroke(1)))
    annotate!([(size(F, 2), i, text(tt(Int(O[i, end])), font(6, :white, "helvetica"))) for i = 1:size(O, 1)])
    plottotals && annotate!([(j, i .+ msm*sqrt(O[i, j])*.0125 .+ 0.02, text(tt(C[i, j]), font(7, ifelse(O[i, j] .> median(O), :red, :red), "helvetica"))) for i = 1:size(C, 1), j = 1:size(C, 2)][:])
    plot!(title="Maximal Proximity Enrichments\nMER11 clusters ↔ dysregulated genes", ylabel="Cluster")
end


function hryticks(p = plot!())
    ytm = yticks(p)[1]
    yticks!(p, (first(ytm), humanread.(first(ytm), string.(first(ytm)))))
end


function plotfreqkde!(bins, y; direction=1, lab="", kwargs...)
    
    
    fdf = length(y)*pdf(kde(y), bins)*(bins[2] - bins[1])
    
    ((direction == 1) || (direction == 3)) && plot!(bins, fdf,  lab=lab ;  kwargs...)
    (direction == 2) && plot!(bins, -fdf, lab=lab ;  kwargs...)
    (direction == 3) && plot!(bins, -fdf; lab="",  kwargs...)

end
RGB255(r, g, b) = RGB(r/255, g/255, b/255)
function plot_conservation_distributions(convdata)
    data = [@subset(convdata, :valid).maxdiff, @subset(convdata, :valid, :DDG).maxdiff, @subset(convdata, :valid, :OMIM).maxdiff]
    bins = range(0, 100, length=100)
    bins=0:.5:100
    p = plot()
    plotfreqkde!(bins, data[1], direction=1, fill=0, c=:black, alpha=0.1, line=stroke(0), lab="All protein coding genes")
    plotfreqkde!(bins, data[3], direction=1, fill=0, c=RGB255(228, 158, 115), alpha=0.65, line=stroke(0), lab="OMIM morbid genes")
    plotfreqkde!(bins, data[2], direction=1, fill=0, c=RGB255(86, 96, 233), line=stroke(0), lab="Developmental Disorder genes")
    plotfreqkde!(bins, @subset(convdata, :NoMam).maxdiff, direction=3, lab="Genes with no Ensemble non-primate ortholog", c=:red, fill=0, alpha=0.75)
    plot!(ylims=(0, 200), grid=false)
    znf808_md = @subset(convdata, :hgnc .== "ZNF808").maxdiff[1]
    vline!([znf808_md], lab="ZNF808", leg=:outertopright, c=:black, ls=:dash)
    plot!(xlims=(0, 100), bottom_margin=8mm, size=(1000, 190), xlabel="Sequence identity difference", ylabel="Gene Frequency", left_margin=5mm, fontfamily="helvetica")
end