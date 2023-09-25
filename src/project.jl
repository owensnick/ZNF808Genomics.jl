## up ver

using RCall
using DataFrames, DataFramesMeta, CSV, Parsers, CodecZlib
using Statistics, StatsBase, GLM
import MultipleTesting
using ProgressMeter, Glob

using GenomicFeatures, ProximityEnrichment, GenomicIntersections, GenomeFragments
using Enrichr

using ClusterOrderTools, Clustering, MultivariateStats, Distances

using Plots, StatsPlots, Measures, Humanize, ImageFiltering, KernelDensity

include("deseq.jl")
include("enrichments.jl")
include("epigenetics.jl")
include("genemodels.jl")
include("pilesave.jl")
include("plots.jl")
include("rnaseq.jl")
include("gtex.jl")
include("conservation.jl")

theme(:wong2)
loaddeseq();


function showwl(lines=200, nc=10000)
    
    table -> begin
        l = ENV["LINES"]
        c = ENV["COLUMNS"]
        ENV["LINES"] = string(lines)
        ENV["COLUMNS"] = string(nc)
        display(table)
        ENV["LINES"] = l;
        ENV["COLUMNS"] = c;
        nothing;
    end
end


function showwide(table, nc=10000)
    c = ENV["COLUMNS"]
    ENV["COLUMNS"] = string(nc)
    display(table)
    ENV["COLUMNS"] = c;
    nothing;
end

function getprojectdir()
    d = pwd()
    if basename(d) == "notebooks"
        return dirname(d)
    else
        return d
    end
end



function loadall()
    

    genetss = loadtss();
       
    lrois = load_loci_of_interest();
    mer11 = loadmer11();
    genemodels = loadgenemodels()
    laser_dorsal_hep = loadlaserde();
    
    human_embryo_rnaseq = CSV.read(joinpath(getprojectdir(), "data", "RNAseq_multitissue.csv.gz"), DataFrame)
    
    meta, counts = loadrnaseq(joinpath(getprojectdir(), "data", "rnaseq_meta.tsv"), joinpath(getprojectdir(), "data", "rnaseq_counts.tsv.gz"));
    meta.ExpNum = replace.(meta.ExpNum, "#" => "E");
    setupdeseq(meta, counts);
    ncounts = normcounts();
    
    enrichrgenes =  loadenrichrgenesets();
    peakmeta, FM = loadepigenetics();
    

    stagelabels = Dict("D0" => "S0", "DE" => "S1", "S2" => "S2", "S3" => "S3", "S4" => "S4");
    gtex = load_gtex()

    convdata = load_conservation_data()

    genetss, lrois, mer11, genemodels, laser_dorsal_hep, human_embryo_rnaseq, meta, counts, ncounts, enrichrgenes, peakmeta, FM, stagelabels, gtex, convdata

end