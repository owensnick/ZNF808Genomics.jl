pile(chrom, loc, FM::FragMatrixSingle{T}, fl=120) where{T}  = 1e+6*fragregion(chrom, loc, "+", FM, inc_meta_frag!, fl, 4, one)/totalfrags(FM)
pile(chrom, loc, FM::FragMatrixPair{T}, minfrag=0, maxfrag=1000) where{T}  = 1e+6*fragregion(chrom, loc, FM, inc_meta_frag!, minfrag, maxfrag, 4, one)/totalfrags(FM)

function loadpiles(file, labels)
    df = CSV.read(file, DataFrame)

    # @assert names(df) == labels
    collect(eachcol(df[!, labels]))
end

function piles_saveload(chrom, loc, inds, labels, FM::Nothing, file; saveload=true)
    !isfile(file) && error("Error: precalcuated $file not found, load alignments to regenerate")
    if saveload == false
        @warn "Loading precalculated signal tracks set to false, but alignments not loaded, loading precalculated signals"
    end
    loadpiles(file, labels)
end
function piles_saveload(chrom, loc, inds, labels, FM, file; saveload=true)
    saveload && isfile(file) && return loadpiles(file, labels)
    #P = @showprogress map(ind -> mapreduce(i -> vec(pile(chrom, loc, FM[i])), +, ind), inds)
    P = @showprogress [pile(chrom, loc, FM[i]) for i in inds]

    if saveload
        mkpath(dirname(file))
        CSV.write(file, DataFrame(P, labels), delim='\t', compress=true)
    end
    P
end
function piles(chrom, loc, inds, grouplabel, labels, FM; saveload=true, projdir=getprojectdir())
    file = joinpath(projdir, "results", "piles", string("pile_", grouplabel, "_", chrom, "_", first(loc), "_", last(loc), ".tsv.gz"))

    piles_saveload(chrom, loc, inds, labels, FM, file, saveload=saveload)
end


piles(regionlabel, grouplabel, datagroup::DataFrame, regioncoords, FM ; saveload=false) = piles(regioncoords[regionlabel]..., datagroup.Index, grouplabel, datagroup.Label, FM, saveload=saveload)