### load

using Parsers

#### Functions for loading gene models
function parsecoords(c)
    f = split(c, r"[, ]", keepempty=false)
    cc = map(c -> Parsers.parse(Int, c), f)
    reshape(cc, 2, div(length(cc), 2))
end

function loadgenemodels(; file="gencode.v36lift37.annotation.sort.coords.tsv.gz", projdir=getprojectdir())
    filepath = joinpath(projdir, "data", file)

    genemodels = CSV.read(filepath, DataFrame, delim='\t')
    genemodels.coords = map(parsecoords, genemodels.coords)

    DataFrames.transform!(groupby(genemodels, :GeneName), [:start, :stop] => ((s, e) -> sortperm(e .- s, rev=true)) => :Rank)

    genemodels
end

function load_loci_of_interest()
    rois = Dict{String, Tuple{String, UnitRange{Int}}}()
    rois["TMTC1"] = ("chr12", 29915923:29975976)
    rois
end


#### Plot Gene models
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
function plotgenemodel!(trans, exoncoords;  y=1, off=1, h = 0.5, label="", fs=10, pos=:left, shiftstrand=false, c=:black, strandarrows=true, ns=15, strandarrowrate=0.1, kwargs...)
    
    g = @subset(exoncoords, :TranscriptID .== trans)
    p = plot!(size=(1200, 100) ; kwargs...)
    
    if shiftstrand    
        if g.strand[1] == "-"
            y = y - 1
            c = :red
        end
    end
    
    
    for (s, e) in eachcol(g.coords[1])
       plot!(rectangle(e - s + 1, h, s - off + 1, y), c=c, linecolor=c, lab="")
    end
    plot!([g.start[1], g.stop[1]] .- off .+ 1, y .+ h/2 .+ [0, 0], c=c, lab="")
    
    
    if ns == -1
        ns = Int(round(strandarrowrate*(g.stop[1] - g.start[1] + 1)/1000))
    end
    rp = range(g.start[1], g.stop[1], length=ns)[2:end-1]
    scatter!(rp, zeros(length(rp)) .+ y .+ h/2, marker=(:black, ifelse(g.strand[1] == "+", :rtriangle, :ltriangle), stroke(0)), lab="")
    
    if !isempty(label)
        if pos == :left
            annotate!([(g.start[1] - 2000, y + h/2, text(label, font(:right, fs, "helvetica", 6)))])
        elseif pos == :right
            annotate!([(g.stop[1] + 2000, y + h/2, text(label, font(:left, fs, "helvetica", 6)))])
        elseif pos == :bottom
            annotate!([((g.start[1] + g.stop[1])/2, y-0.5, text(label,font(:center, fs, "helvetica", 6)))])
        elseif pos == :top
            annotate!([((g.start[1] + g.stop[1])/2, y+0.5, text(label,font(:center, fs, "helvetica", 6)))])
        end
    end
    p
end