
"""
    dataToFit_exp(dataname, cloneid)

Prepare the experimental data for the fit
"""
function dataToFit_exp(dataname, cloneid)
    #load the data for the fit
    metadata = load(abspath(@__DIR__, "../../processed_data/$(dataname).jld2")) ;
    bstat = DataFrame(Arrow.Table(abspath(@__DIR__, "../../processed_data/$(dataname).arrow"))) ;
    leftlimit = metadata["leftlimit"]
    rightlimit = metadata["rightlimit"]



    #probability to observe a burst from LiveCell 
    data_pburst = mean(bstat.pburst[cloneid][leftlimit:rightlimit])


    #survival probabilities of the  BURST and INTER-BURST times, from LiveCell
    timevec_interburst = bstat.km_interburst[cloneid].events.time
    survivalinterburst_data = bstat.km_interburst[cloneid].survival
    timevec_burst = bstat.km_burst[cloneid].events.time
    survivalburst_data = bstat.km_burst[cloneid].survival

    
    #correlation between two consecutive inter-burst events, from LiveCell
    consecutive_evts = reshape(copy(bstat.eventcorr_interburst[cloneid]),Int(length(bstat.eventcorr_interburst[cloneid])/3),3)
    
    #mean of the first and second event times
    m = mean(consecutive_evts,dims=1)
    #sd of the first and second event timess
    sd = sqrt.(var(consecutive_evts,dims=1))
    global mutemp = 0
    for i in axes(consecutive_evts,1)
        global mutemp += (consecutive_evts[i,1] - m[1])*(consecutive_evts[i,2] - m[2])
    end
    data_corr_interburst =1/size(consecutive_evts,1) * mutemp/(sd[1]*sd[2])

    return ((timevec_burst,survivalburst_data,), (timevec_interburst,survivalinterburst_data,), data_pburst,  data_corr_interburst, )

end



"""
    dataToFit(pathtotracks)

Prepare the simulated data for the fit
"""
function dataToFit_sim(survivalanalysis_table,simidx)
    #probability to observe a burst from LiveCell 
    data_pburst = survivalanalysis_table.pburststar[simidx]

    #survival probabilities of the  BURST and INTER-BURST times, from LiveCell
    timevec_interburst = survivalanalysis_table.km_interburst[simidx].events.time
    survivalinterburst_data = survivalanalysis_table.km_interburst[simidx].survival
 
    timevec_burst = survivalanalysis_table.km_burst[simidx].events.time
    survivalburst_data = survivalanalysis_table.km_burst[simidx].survival
 
    #correlation between two consecutive inter-burst events, from LiveCell
    consecutive_evts = reshape(copy(survivalanalysis_table.eventcorr_interburst[simidx]),Int(length(survivalanalysis_table.eventcorr_interburst[simidx])/2),2)
    #mean of the first and second event times
    m = mean(consecutive_evts,dims=1)
    #sd of the first and second event timess
    sd = sqrt.(var(consecutive_evts,dims=1))
    global mutemp = 0
    for i in axes(consecutive_evts,1)
        global mutemp += (consecutive_evts[i,1] - m[1])*(consecutive_evts[i,2] - m[2])
    end
    data_corr_interburst =1/size(consecutive_evts,1) * mutemp/(sd[1]*sd[2])
    
    return ((timevec_burst,survivalburst_data,), (timevec_interburst,survivalinterburst_data,), data_pburst,  data_corr_interburst)
end


#= function dataToFit_simbased(simname, simidx)
    #load the data for the fit
    data_table = Arrow.Table("/tungstenfs/scratch/ggiorget/Gregory/julia/bursting_manuscript/src/model_simulations/data_forfit/$(simname).arrow") ;
    
  
    #probability to observe a burst from LiveCell 
    data_pburst = data_table.pburststar[simidx]

   
    #survival probabilities of the  BURST and INTER-BURST times, from LiveCell
    timevec_interburst = data_table.km_interburst[simidx].events.time
    survivalinterburst_data = data_table.km_interburst[simidx].survival
    survivalinterburst_data_std = data_table.km_interburst[simidx].stderr

    timevec_burst = data_table.km_burst[simidx].events.time
    survivalburst_data = data_table.km_burst[simidx].survival
    survivalburst_data_std = data_table.km_burst[simidx].stderr
  
    #correlation between two consecutive inter-burst events, from LiveCell
    consecutive_evts = reshape(copy(data_table.eventcorr_interburst[simidx]),Int(length(data_table.eventcorr_interburst[simidx])/2),2)
    #mean of the first and second event times
    m = mean(consecutive_evts,dims=1)
    #sd of the first and second event timess
    sd = sqrt.(var(consecutive_evts,dims=1))
    global mutemp = 0
    for i in axes(consecutive_evts,1)
        global mutemp += (consecutive_evts[i,1] - m[1])*(consecutive_evts[i,2] - m[2])
    end
    data_corr_interburst =1/size(consecutive_evts,1) * mutemp/(sd[1]*sd[2])

    return ((timevec_burst,survivalburst_data,survivalburst_data_std,), (timevec_interburst,survivalinterburst_data,survivalinterburst_data_std,), data_pburst,  data_corr_interburst)

end =#

