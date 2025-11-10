#Configuration file for the fit
#-------------------------------------
using DelimitedFiles

#MODEL configuration
model= modelInstance()

#DATA configuration
    #define the detection limit fot the Live Cell experiment and for the Fixed Cell experiment for the nascent RNA count
    detectionLimitLC = 1
    detectionLimitNS = 2
    
#FIT configuration
    #decide to take into account or not the burst singlets
    burstsinglet = :without
    #define the data type used in the fit
    datatype = (StoThyLiveCell.Survival_Burst(),StoThyLiveCell.Survival_InterBurst(),StoThyLiveCell.Survival_InterBurst(), StoThyLiveCell.Prob_Burst(), StoThyLiveCell.Correlation_InterBurst(),)
    #define the group of data type (i.e. :LiveCell if only live-cell data, :rna_dist if only FISH mRNA dist, :mixed if both)
    datagroup = StoThyLiveCell.LiveCellData()
    #define the distance function for each data type
    dist= (StoThyLiveCell.LsqSurvival(.01), StoThyLiveCell.LsqSurvivalLogLin(1.,.8),StoThyLiveCell.LsqSurvivalLogLin(10.,.8),  StoThyLiveCell.LsqProb(.1), StoThyLiveCell.LsqNumber(1.5))
    #define the data for each data type
    # 1= survival_burst, 2=survival_interburst, 3=survival_nextburst, 4=intensity_burst, 5=mean_nascent, 6=prob_burst, 7=correlation_interburst, 8=distribution_RNA)
    dataForFit_idx = [1,2,2,3,4]

 

    #define the maximum number of mRNA for the Live Cell model and the steady-state number of mRNA model (FC)
    maxrnaLC = 10
    maxrnaFC = 40 #ATTENTION this should be larger or equal to the FRange corresponding to distribution_RNA
    tmax_nextburst = 200
    tmax_intensity = 6

    #define range that you want to fit for each data type
    FRange = [(0,7),(0,290),(0,60),(0,0),(0,0)]
 

    #define the parameter range for the fit
    SRange = [(0.00025,0.01),(0.00025,.01),(0.000025,.1),(0.001,1.),(0.001,1.),(0.01,1.0),(0.1,1.0),]






    #indices of the free parameters
    freeparametersidx = [1,2,3,4,5,6,7]

    #define the value of the fixed parameters, -1 if there is no
    fixedparameters = [-1]
    #set the initial parameters
    initialparameters =[]
    function iniparam()
        pathTomodels = "/tungstenfs/scratch/ggiorget/Gregory/julia/FitStoThyLiveCell/src/models"
        pathToparameters = pathTomodels*"/m_2s2r/run_20250218_5/log_bestfits.txt"
        pathTofval = pathTomodels*"/m_2s2r/run_20250218_5/log_fval.txt"
        paratemp = readdlm(pathToparameters)
        fvaltemp = vec(readdlm(pathTofval))
        fvalsort = copy(fvaltemp)
        sort!(fvalsort)
        threshold = fvalsort[1]
        idxs = findall(fvaltemp .<= threshold)
        initialparameters = paratemp[idxs,:]
#=         for i =2:5
            pathToparameters = pathTomodels*"/m_2s2r/run_20250210_$i/log_bestfits.txt"
            pathTofval = pathTomodels*"/m_2s2r/run_20250210_$i/log_fval.txt"
            paratemp = readdlm(pathToparameters)
            fvaltemp = vec(readdlm(pathTofval))
            fvalsort = copy(fvaltemp)
            sort!(fvalsort)
            threshold = fvalsort[2]
            idxs = findall(fvaltemp .<= threshold)
            initialparameters = vcat(initialparameters,paratemp[idxs,:])
        end =#
        return initialparameters
    end
    initialparameters = iniparam()

    #optimization parameters
    maxoptimtime =100#2000
    maxoptimiters =100#2000
    nboptims = 1
    autodiff = true

    using Optimization
    using OptimizationOptimJL
    method = BFGS()




