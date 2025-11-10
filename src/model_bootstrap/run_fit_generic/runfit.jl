function runfit(data_all, pathToLog, pathToBestfit)
    datatofit = data_all[dataForFit_idx]

    #set up the optimization based on the config file.
    dataForFit = StoThyLiveCell.DataFit{typeof(datatype),typeof(datatofit)}(datatype,datatofit,detectionLimitLC, detectionLimitNS, burstsinglet)
    optimStruct = StoThyLiveCell.OptimStruct{typeof(dataForFit), typeof(dist), typeof(model)}(dataForFit,dist,model, autodiff)

    #run the optimization
    if @isdefined method
        sol, bfparameters, minval, minidx, estimate_signal = StoThyLiveCell.optim_function(SRange, FRange, optimStruct; NbOptim=nboptims, maxtime=maxoptimtime, maxiters=maxoptimiters, fixedparameters=fixedparameters,  freeparametersidx=freeparametersidx, maxrnaLC=maxrnaLC, maxrnaFC=maxrnaFC, Method=method, pathToLog=pathToLog, initialparameters=initialparameters)
    else
        sol, bfparameters, minval, minidx, estimate_signal = StoThyLiveCell.optim_function(SRange, FRange, optimStruct; NbOptim=nboptims, maxtime=maxoptimtime, maxiters=maxoptimiters, fixedparameters=fixedparameters,  freeparametersidx=freeparametersidx, maxrnaLC=maxrnaLC, maxrnaFC=maxrnaFC, pathToLog=pathToLog,initialparameters=initialparameters)
    end


    objfunc = []
    nbfeval = []
    retcode = []
    timesec = []
    for i in eachindex(sol)
        push!(objfunc, sol[i].objective)
        push!(nbfeval, sol[i].stats.fevals)
        push!(retcode, sol[i].retcode)
        push!(timesec, sol[i].stats.time)
    end
    solstats = (objfct = objfunc, retcode=retcode, nbfeval=nbfeval, timesec=timesec,)
    println(bfparameters)
    println("fobj = $minval" )
    #save the optimization outputs
    save(pathToBestfit, "solstats", solstats, "bfparameters", bfparameters, "minobj", minval, "minidx", minidx, "estimate_signal", estimate_signal)
end
    