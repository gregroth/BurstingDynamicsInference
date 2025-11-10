#Configuration file for the fit
#-------------------------------------
using DelimitedFiles


#define the dataset that you want to use for the fit

simname = "SIMNAMEHERE"
fitname = "FITNAMEHERE"

include(joinpath(@__DIR__, "../../model_fit_MODELID/modelBase.jl"))
model_th = modelInstance()


#define the maximum number of mRNA for the Live Cell model and the steady-state number of mRNA model (FC)
maxrnaLC = 10
maxrnaFC = 40 #ATTENTION this should be larger or equal to the FRange corresponding to distribution_RNA
tmax_nextburst = 200
tmax_intensity = 6


# 1= survival_burst, 2=survival_interburst,  3=prob_burst, 4=correlation_interburst)
dataForFit_idx = [1,2,4]

#define range that you want to fit for each data type
FRange = [(0,50),(0,200), (0,0)]

#define the parameter range for the fit
SRange = [(0.00025,0.1),(0.00025,.1),(0.000025,.1),(0.01,10.),(0.01,10.),(0.01,10.0),(0.1,5.0),]



