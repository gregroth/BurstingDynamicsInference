#Configuration file for the fit
#-------------------------------------
using DelimitedFiles
using Tables, FileIO, JLD2
using DataFrames
using Arrow
using Random

using LinearAlgebra
using StatsBase

using Survival

include("../StoThyLiveCell/StoThyLiveCell.jl") 
using .StoThyLiveCell



include("../model_fit_2s2r/modelBase.jl")


include("../model_simulation/m_2s2r.jl")
using .m_2s2r


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
    datatype = (StoThyLiveCell.Survival_InterBurst(), StoThyLiveCell.Prob_Burst(),)
    #define the group of data type (i.e. :LiveCell if only live-cell data, :rna_dist if only FISH mRNA dist, :mixed if both)
    #define the distance function for each data type
    dist = (StoThyLiveCell.LsqSurvivalLogLin(1.,.8), StoThyLiveCell.LsqProb(.01),)
    #define the data for each data type
    # 1= survival_burst, 2=survival_interburst, 3=survival_nextburst, 4=intensity_burst, 5=mean_nascent, 6=prob_burst, 7=correlation_interburst, 8=distribution_RNA)
    dataForFit_idx = [2,3]
    #define the maximum number of mRNA for the Live Cell model and the steady-state number of mRNA model (FC)
    maxrnaLC = 10
    maxrnaFC = 55 #ATTENTION this should be larger or equal to the FRange corresponding to distribution_RNA
    tmax_nextburst = 200
    tmax_intensity = 6

    #define range that you want to fit for each data type
    FRange = [(0,290),(0,0),]
    #define the parameter range for the fit
    SRange = [(0.000025,0.1),(0.00025,.1),(0.000025,.01),(0.001,.1),(0.001,.5),(0.01,1.0),(0.1,2.0),]
    #indices of the free parameters
    freeparametersidx = [1]

    bestfitTrialID_ref = "fcfca7c6_p"*"_cl$cloneid"
    pathToBestfit_ref = "/tungstenfs/scratch/ggiorget/Gregory/julia/FitStoThyLiveCell/src/models/m_2s2r/run_allcl_20250218_11_p1/bestfits/$(bestfitTrialID_ref).jld2"
    bf= load(pathToBestfit_ref) ;
    ref_parameters = [bf["bfparameters"][1:7]];
    fixedparameters = bf["bfparameters"][2:7]

    initialparameters = []

    #optimization parameters
    maxoptimtime =5000
    maxoptimiters =5000
    nboptims = 2
    autodiff = false

    using Optimization
    using OptimizationOptimJL
    method = SAMIN()




