#--------------------------------------------
using Tables, FileIO, JLD2
using DataFrames
using Arrow
using Random
using Glob
using Survival

include("../StoThyLiveCell/StoThyLiveCell.jl") 
using .StoThyLiveCell

#set the random seed 
Random.seed!(1*11)

# model
modelname = :m_2s2r
model = replace(String(modelname),"m_" => "")

#set the minumum number of RNA detectable
minnbrna = 1


#simulation id
simname = "sim_$(model)_$(minnbrna)"
fitname = simname 

include("../model_fit_$(model)/modelBase.jl")


include("../model_simulation/m_$(model).jl")
@eval using .$modelname: simulate_bursttraj


parameters_data = load("./selected_parameters/parameterlist_filtered_$(model).jld2")
parameter_list = parameters_data["prameterlist"]
parameter_list = parameter_list[1:20]
#parameter range
#3s_1r
if model == "3s1r"
    SRange = [(0.00025,10.),(0.00025,10.),(0.00025,10.),(0.00025,10.),(0.01,10.0),(0.1,5.0),]
elseif model == "2s2r"
    SRange = [(0.00025,0.1),(0.00025,.1),(0.000025,.1),(0.01,10.),(0.01,10.),(0.01,10.0),(0.1,5.0),]
else
    println("SRange not specified")
end


#probabilties to miss a weak burst
pmiss = 0:.1:1

#total number of simulation
NbtotSim = string(length(parameter_list)*length(pmiss))

#upload the track length distribution from experimental data
tl = load("../../processed_data/tracklength_dist_5G7.jld2")
tracklengthlist = tl["tracklengths"]
NbTrack = length(tracklengthlist)


configvec = (model=model, paraRange=SRange, nbsim = NbtotSim, nbtrack=NbTrack)

println("starting simulations")
include("../model_simulation/simulation_generic.jl")

for i in eachindex(parameter_list)
    parameters = copy(parameter_list[i])
    trajmodel = simulation_track(tracklengthlist, parameters)
    save("./data_simulations/sim_$(model)/trajset_$i.jld2", "trajmodel", trajmodel)
end

println("simulations done")



println("adding perturbation, and preparing the data for fit")
include("data_analysis.jl")

pathtotrajectories = "./data_simulations/sim_$(model)/"
pathtosave = "./data_forfit/"
data_analysis(pmiss, minnbrna, configvec, parameter_list, pathtosave,pathtotrajectories)
println("done")



println("launching the fits")

run(`tmux new-window -t session_simulations -n script_run "bash /tungstenfs/scratch/ggiorget/Gregory/julia/BurstingDynamicsInference/src/modelassumption_robustness/run_fit_simbased_batch_generic.sh $simname $NbtotSim $model"`)



# Wait for the done file
donefile = "/tungstenfs/scratch/ggiorget/Gregory/julia/BurstingDynamicsInference/src/modelassumption_robustness/myscript_done.flag"
while !isfile(donefile)
    sleep(1)  # wait 1 second
end

# Clean up the flag (optional)
println("the done flag has been detected and removed")
rm(donefile)

println("the fit run files are removed")

for i=1:parse(Int, NbtotSim)
    rm("/tungstenfs/scratch/ggiorget/Gregory/julia/BurstingDynamicsInference/src/modelassumption_robustness/fitsim/runfit_$i.jl")
end

println("Script exp fit completed! Continuing...")

println("starting the fit summary")

include("summary_fits_generic.jl") 

