

function loaddeseq()
    R"""
    suppressMessages(library(DESeq2))
    suppressMessages(library(data.table))
    """
end

function loaddesign(design)
    @rput design
    
    R"""
    design            <- as.data.frame(design)
    design[]          <- lapply(design, factor)
    row.names(design) <- design$Sample
    """

end

function loadcounts(counts)
    countdata     = Matrix{Int}(counts[:, 2:end]) 
    genes = counts.Gene
    samples = string.(names(counts)[2:end])
    @rput samples countdata genes
    R"""
        row.names(countdata) <- genes
        colnames(countdata)  <- samples
    """
end

function setupdeseq(design, counts, formula=R"~Group")
    loaddesign(design)
    loadcounts(counts)

    R"""
        dds <- DESeqDataSetFromMatrix(countData = countdata, colData = design, design = $formula);
    """
end

function normcounts()
    R"""
        dds <- estimateSizeFactors(dds)
        ncX <- counts(dds, normalized=TRUE)
        genes <- rownames(countdata)
        samples <- colnames(countdata)
    """;
    @rget ncX genes samples
    normcounts = [DataFrame(Gene=genes) DataFrame(ncX, samples)]
end
deseq(; obj=R"dds") = R"$obj <- DESeq($obj, quiet=TRUE)"

function local_sample_dict()
    localdict=Dict("D0" => ["D0", "DE"],
    "DE" => ["D0", "DE", "S2"],
    "S2" => ["DE", "S2", "S3"],
    "S3" => ["S2", "S3", "S4"],
    "S4" => ["S3", "S4"])

    localdict
end

function deseq_batch_local(meta, counts, genetss, sampledict=local_sample_dict(); fct=log2(1.25))
    design, filtcounts = filterdesign(meta, counts);
    stages = unique(design.Stage)

    setupdeseq(design, filtcounts, R"~ExpNum + Group")
    ncounts = normcounts()
    dfs = DataFrame[]
    println("Running DESeq sharing information between local stages, experiment batch effect...")
    @showprogress "DESeq by Stage: " for s in stages
        designstage, countstage = filterdesign(@subset(design, :Stage .∈ Ref(sampledict[s])), counts)
        setupdeseq(designstage, countstage,  R"~ExpNum + Group")
        # deseq()
        R"""
            dds <- DESeq(dds, quiet=TRUE)
        """
        
        degs = deseq_pairwisegroup_batch([s], independentFiltering=true, fct=fct);
        push!(dfs, degs)
    end
    deg = reduce(vcat, dfs)
    deg = annotatetss(deg, genetss)
    meancounts = calcmeancounts(design, ncounts);
    
    (design=design, filtcounts=filtcounts, stages=stages, ncounts=ncounts, meancounts=meancounts, deg=deg)
end



function collectresults(label="")
    R"""
    rdf <- as.data.frame(res)
    setDT(rdf, keep.rownames=TRUE)[]
    colnames(rdf)[1] <- "Gene"
    """
    
    @rget rdf
    
    if !isempty(label)
        rename!(rdf, ["Gene" ; string.(names(rdf)[2:end], "_", label)])
    end
    rdf
end


"""
    deseq_pairwisegroup_batch(meta)

Build pairwise comparisons data must of have loaded with ~ExpNum + Group design.
"""
function deseq_pairwisegroup_batch(stages; independentFiltering=false, padj_thresh=0.05, fct=log2(1.25), obj=R"dds")
    dfs = DataFrame[]
    @showprogress "Collecting DESeq2 results: " for stage in stages
        #reslabel = string("Group_WT_", stage, "_vs_KO_", stage)
        stlabel_WT = string("WT_", stage)
        stlabel_KO = string("KO_", stage)
        #res <- results($obj, name=$reslabel, independentFiltering=$independentFiltering)
        R"""
            res <- results($obj, contrast=c("Group", $stlabel_KO, $stlabel_WT), independentFiltering=$independentFiltering)
            disp <- dispersions($obj)
        """
        res = collectresults()
        #res.log2FoldChange .= .-res.log2FoldChange ## comparison is WT/KO, flip to KO/WT
        res[!, :Stage] .= stage
        @rget disp
        res[!, :Dispersion] = disp
        push!(dfs, res)
    end
    df = reduce(vcat, dfs)
    df[!, :sig]  = (coalesce.(df.padj, 1.0) .< padj_thresh) .& (abs.(df.log2FoldChange) .> fct)
    df[!, :dsig] = df.sig.*sign.(df.log2FoldChange)
    df[!, :GeneName] = last.(split.(df.Gene, "|"))
    df
end
    