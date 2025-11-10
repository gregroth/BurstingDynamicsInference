

include("perturbation.jl")
include("../data_analysis/survival_analysis.jl")
using .survival_analysis


function data_analysis(pmiss, minnbrna, config, parameterlist, pathtosave,pathtotrajectories)
    bstat = DataFrame(uid = Vector{Int64}(), parameters = Vector{Vector{Float64}}(), pmiss = Vector{Float64}(), config = Vector{NamedTuple{}}(),
    km_burst = Vector{KaplanMeier{Float64, Int64}}(), km_interburst = Vector{KaplanMeier{Float64, Int64}}(), 
            meanburstfreq = Vector{Float64}(),pburststar =  Vector{Float64}(), 
                    eventcorr_burst = Vector{Array{Int64}}(), eventcorr_interburst = Vector{Array{Int64}}()
                        )
    uid = 1
    trajfiles = glob("*.jld2",pathtotrajectories)
    for p in pmiss
        for i in eachindex(trajfiles)
            data = load(trajfiles[i])
            trajmodel_set = data["trajmodel"]
            trajnoise_set = remove_weakbursts(trajmodel_set,p, minnbrna)
            bstat_temp =survival_analysis.survival_wosinglet(trajnoise_set, 1);

            push!(bstat, (uid, parameterlist[i], p, config,
                                bstat_temp.km_burst[1], bstat_temp.km_interburst[1],  
                                    bstat_temp.meanburstfreq[1], bstat_temp.pburststar[1], 
                                        bstat_temp.eventcorr_burst[1], bstat_temp.eventcorr_interburst[1]));
            uid += 1
        end
        Arrow.write("$(pathtosave)$(simname).arrow", bstat)
        println("pmiss $p is done")
    end
end
