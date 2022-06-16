# Primate-specific ZNF808 is essential for pancreatic development in humans
## https://doi.org/10.1101/2021.08.23.21262262
## Notebook to perform epigenetic and transcriptomic analysis to generate figure panels.
# This notebook clusters the H3K9me3 and HEK27ac epigenetic data over pancreatic differentiation time course and then analyses transcriptomic data collected at same stages in light of clusters.
# Note repository contains source data for all gene expression and quantifications of ChIP-seq at particular loci. Repository does not contain alignment files, but code is supplied to regenerate quantifications if suitable alignment files are present.

include(joinpath("..", "src", "project.jl"))
ENV["GKSwstype"] = "100" ### this needed for GR when running on server without a display terminal, comment out should you want GR to display figures rather than saving

## Differential Expression analysis to calculate dysregulation genes in ZNF808 KO with proximity and gene set enrichments
genetss, lrois, mer11, genemodels, laser_dorsal_hep, human_embryo_rnaseq, meta, counts, ncounts, enrichrgenes, peakmeta, FM, stagelabels = loadall();

sampleselection = @subset(meta, :Run .== "R4")
    localdict=Dict("D0" => ["D0", "DE"],
                    "DE" => ["D0", "DE", "S2"],
                    "S2" => ["DE", "S2", "S3"],
                    "S3" => ["S2", "S3", "S4"],
                    "S4" => ["S3", "S4", "S5"],
                    "S5" => ["S4", "S5", "S6"],
                    "S6" => ["S5", "S6", "S7"],
                    "S7" => ["S6", "S7"])
@time deresults =  enrichments(deseq_batch_local(sampleselection, counts, genetss, localdict), laser_dorsal_hep, enrichrgenes);
saveenrichr(deresults.enfisher, stagelabels);


### Clustering of Epigenetic data

# k-means clustering applied to H1 H3K9me3, H3K27ac WT and KO data.

epiclusters = clusterpeaks(deresults, mer11);

plot!(epiclusters.silhouetteplot)
savedisplay("epi_clustter_silhouette")
#Median Silhoette score suggests k = 5, offers a good compromise. At k = 5, clusters 1 and 2 reveal comparable H3K27ac behaviour, but differ by H3K9me3. Recode as 4 clusters with differing H3K27ac behaviour:

ph = heatmap(epiclusters.ZMS[epiclusters.mcl.Cluster5Order, :])
clusterlines!(epiclusters.kd["All_CL5"])
plot!(yflip=true)
@with @subset(peakmeta, :GenoType .!= "iPSC") plot!(xticks=(1:length(:Chip), string.(:Chip, "\n", :GenoType, "\n", :Stage)), size=(1000, 500), bottom_margin=10mm, title="Original clusters")

epiclusters.mcl[!, :Cluster5Map] = Int.(floor.(epiclusters.mcl.Cluster5Label[epiclusters.mcl.Cluster5Order]))
clstats = combine(groupby(epiclusters.mcl, :Cluster5Map, sort=true), nrow => :count)
pc = heatmap(epiclusters.mcl.Cluster5Map[:, :], yflip=true, colorbar=false, c=:blues, size=(100, 300), xticks=false, yticks=false, framestyle=:none, title="Recoded\nClusters")
hline!(cumsum(clstats.count), lab="", c=:white)
annotate!([(1, v, text(i, font("helvetica", :white))) for (i, v) in enumerate(cumsum(clstats.count) .- clstats.count/2)])

plot(ph, pc, layout=grid(1, 2, widths=[0.95, 0.05]), size=(900, 500), titlefont=font(10, "helvetica"), top_margin=5mm)

@with epiclusters.mcl annotatemerdist!(deresults.deg, Interval.(:chrom, :start, :stop), "dMER11");
@with @subset(epiclusters.mcl, :All_CL5 .== 5) annotatemerdist!(deresults.deg, Interval.(:chrom, :start, :stop), "dMER11_CL5");
deresults.deg[!, :pMER11] = coalesce.(deresults.deg.dMER11, 1e+8) .< 1e+6;

design_ips = [sampleselection ; @subset(meta, :Clone .== "Patient")];
design_ips_full = [design_ips ; @subset(meta, :Clone .== "H1", :Run .== "R3")]
design_ips.GenoTypeIPS = ifelse.(design_ips.Clone .== "Patient", "iPSC", design_ips.GenoType);
design_ips_full.GenoTypeIPS = ifelse.(design_ips_full.Clone .== "Patient", "iPSC", design_ips_full.GenoType);


# Figure 4 panels

fig4A = plottotalde(deresults.deg, stagelabels)


fig4B_inset = proxplot(@subset(epiclusters.proxmcl, :Stage .== "D0"), ["MER_ALL"], labs=["Activated" "Repressed"], dx=1)
plot!(size=(200, 150), leg=false, yticks=0:15:30, ylabel="", title="S0", xticks=([4, 6, 8], ["10kb", "1Mb", "10Mb"]))
savedisplay("fig4b_inset")

pp = proxplot(epiclusters.proxmcl, ["MER_ALL"], labs=["Activated" "Repressed"], dx=2)
plot!(framestyle=:default)
plot!(pp.subplots[1], yticks=0:15:30)
for p in pp.subplots[2:end]
    plot!(p, yaxis=false)
end
[plot!(p, title="", xlabel=stagelabels[s]) for (p, s) in zip(pp.subplots, deresults.stages)]
fig4B = plot!(pp, xticks=([4, 6, 8], ["", "1Mb", ""]), right_margin=-1mm, left_margin=-2mm)

fig4C = epiproxproption(deresults)
yticks!(0:.1:.4, string.(0:10:40, "%"))
plot(fig4A, fig4B, fig4C, layout=(1, 3), size=(1000, 300), fontfamily="helvetica", bottom_margin=5mm)

savedisplay("fig4_abc")


fig4D_left = plot_gene_snapshot(lrois["TMTC1"]..., @subset(peakmeta, :Stage .∈ Ref(["D0"])), FM, mer11, genemodels, label="epi_h1", merex=1000, dx=50, transcripts=["ENST00000539277.6_8"], plottype=:stack, yld=Dict("K9" => 2.9, "K27" => 1.95))
plot!(size=(600, 300), fontfamily="helvetica") 
savedisplay("fig4_d_left")

fig4D_right = plotgenesimple(r"\|TMTC1$", deresults, stagelabels, meta = @subset(design_ips, :Stage .∈ Ref(["D0", "DE", "S2", "S3"])), ncounts=ncounts)
plot!(size=(300, 287), fontfamily="helvetica", leg=:topright)
savedisplay("fig4_d_right")


fig4e = panc_hep_enrichment_plot(deresults, laser_dorsal_hep)

# Supplemental Figure Panels

sup_fig1b = plot_human_embryo_rnaseq(human_embryo_rnaseq, r"ZNF808|ZNF440|ZNF525|ZNF468|ZNF433|ZNF578");
plot!(fontfamily="helvetica", title="Expression of MER11-binding KRAB ZFPs in human embryo (RNA-seq)")

sup_fig7ab = geneheatmap_scatter(r"ZNF(808|525|440|468|433|578)$", deresults, meta = design_ips_full, ncounts=ncounts, size=(1000, 250), fontfamily="helvetica", titlefont=font(10, "helvetica"))


sup_fig7cd = geneheatmap_scatter(r"GATA(1|3|4|6)$|HNF4A|HNF4G|\|DUX$|\|NANOG$",  deresults, meta = design_ips_full, ncounts=ncounts, k=6, size=(1000, 250), fontfamily="helvetica", titlefont=font(10, "helvetica"))


supp_fig9 = ips_comparison(deresults, design_ips_full, ncounts)

supp_fig_10a = epi_cluster_proximity_enrichments(epiclusters, deresults)

supp_fig_10b = enrichrplot(deresults.enfisher, ["Fetalliver", "Liver", "TGF-beta regulation of extracellular matrix"],
    ["PancreaticIslet", "Pancreas", "regulation of transcription from RNA polymerase II promoter (GO:0006357)", "negative regulation of cell differentiation (GO:0045596)", "TGF-beta regulation of extracellular matrix"], 10)


pa = plotgenesimple(r"\|AFP$", deresults, stagelabels, meta = @subset(deresults.design, :Stage .∈ Ref(["D0", "DE", "S2", "S3", "S4", "S5", "S6", "S7"])), plotpair=true)#, ncounts=ncounts, trans = x -> x)
hryticks()
pd = plotgenesimple(r"\|PDX1$", deresults, stagelabels, meta = @subset(deresults.design, :Stage .∈ Ref(["D0", "DE", "S2", "S3", "S4", "S5", "S6", "S7"])), plotpair=true)#, ncounts=ncounts, trans = x -> x)
hryticks()
sup_fig_11b = plot(pa, pd, leg=:topleft, size=(700, 300))
