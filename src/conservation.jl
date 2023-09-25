

function load_conservation_data(; file="ZNF808_conservation_data_22Nov2022.csv.gz", projectdir=getprojectdir())
    convdata = CSV.read(joinpath(projectdir, "data", file), DataFrame);
    rename!(convdata, "(ape+prim)-npmam" => "maxdiff", "morbid-MW" => "morbidMW")
    convdata.valid = ((convdata.max_ape .!= 0) .| (convdata.max_primate .!= 0)) .& (convdata.max_npmam .!= 0) # .& (convdata.maxdiff .>= 0)
    convdata.NoMam = ((convdata.max_ape .!= 0) .| (convdata.max_primate .!= 0)) .& (convdata.max_npmam .== 0) 
    convdata.DDG = (convdata.ddg2p .!= "0") .& (convdata.hgnc .!= "ZNF808")
    convdata.OMIM = (convdata.morbidMW .!= "0") .& (convdata.hgnc .!= "ZNF808");

    convdata
end