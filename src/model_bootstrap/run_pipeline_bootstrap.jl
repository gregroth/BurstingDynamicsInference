
runfitid = 1

cloneid = 1

simid = "runfit_$runfitid"
simname = "sim_2s2r_$(simid)"
fitname = simname



if cloneid == 1 
    include("config_5kb.jl")
else
    include("config_otherclones.jl")
end 


include("../data_analysis/survival_analysis.jl")
using .survival_analysis

include("../utiles/dataToFit.jl")



#upload the track length distribution from experimental data
tl = load("../../processed_data/tracklength_dist_5G7.jld2")
tracklengthlist = tl["tracklengths"]
NbTrack = length(tracklengthlist)


    #set up the data to fit
    include("./run_fit_generic/runfit.jl")

for i = 1:2
    println("starting the process for sim $i")

    println("starting simulations")
    include("../model_simulation/simulation_generic.jl")
    if !isdir("./data_simulations/sim_2s2r_$(simid)")
        mkdir("./data_simulations/sim_2s2r_$(simid)")
    end
    trajmodel = simulation_track(tracklengthlist, ref_parameters)
    #save("./data_simulations/sim_2s2r_$(simid)/trajset.jld2", "trajmodel", trajmodel)
    
    println("simulations done")



    println("preparing the data for fit")
    #data = load("./data_simulations/sim_2s2r_$(simid)/trajset.jld2")
    #trajmodel_set = data["trajmodel"]
    data_table= survival_analysis.survival_wosinglet(trajmodel, 1);
   
    data_all = dataToFit_sim(data_table,1)



    println("launching the fits")
    pathToLog = "./run_fit_log/"
    pathToBestfit = "./bestfit_list/sim_clone_$(cloneid)_runfit_$(runfitid)_$i.jld2"
  
    runfit(data_all, pathToLog, pathToBestfit)
end

