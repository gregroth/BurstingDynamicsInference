# Fit of the Kaplan-Meier survival curve to the two-state modelexp
#
# This script fit a given KM-modelexp curve to the KM-data curve
#
# Gregory Roth
# First version:            25.04.2023
# Modified:                 25.04.2023
#
#-----------------------------------------------------------------------------------------------
using Tables, FileIO, JLD2
using StatsBase
using Survival
using DataFrames
using Random
using LsqFit
using Arrow
using LinearAlgebra


include(joinpath(@__DIR__, "../../StoThyLiveCell/StoThyLiveCell.jl")) 
using .StoThyLiveCell

#set the random seed 
Random.seed!(1*11)

include("config.jl")
#load the data for the fit
include(joinpath(@__DIR__, "../../utiles/dataToFit.jl"))

include(joinpath(@__DIR__, "../theoretical_timescales.jl"))


data_table = Arrow.Table(joinpath(@__DIR__,"../data_forfit/$(simname).arrow")) ;




simuid = SIMIDX
data_all = dataToFit_sim(data_table,simuid) 
datatofit = data_all[dataForFit_idx]

if minimum(data_all[1][2])<=0
    idx = findall(dataForFit_idx.==1)
    for i in idx
        FRange[i] = (FRange[i][1],minimum([FRange[i][2],datatofit[i][1][findfirst(data_all[1][2].<=0)-1]]))
    end
end

if minimum(data_all[2][2])<=0
    idx = findall(dataForFit_idx.==2)
    for i in idx
        FRange[i] =  (FRange[i][1],minimum([FRange[i][2],datatofit[i][1][findfirst(data_all[2][2].<=0)-1]]))
    end
end



#correlation interburst sim data
correlation_interburst_sim = datatofit[3]
#theoretical correlation
(survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst_th, intensity_burst) =    StoThyLiveCell.ModelOutput_wosinglet(model_th, copy(data_table.parameters[simuid]), 40, 1,1, 300,300,200,10) 




fit_df = DataFrame(uid = Vector{Int64}(), parameters = Vector{Vector{Float64}}(), pmiss = Vector{Float64}(), 
config = Vector{NamedTuple{}}(),
top_thy_ts_off = Vector{Vector{Tuple}}(),
top_thy_ts_on= Vector{Vector{Tuple}}(),
sim_ts_off = Vector{Vector{Tuple}}(),
sim_ts_on= Vector{Vector{Tuple}}(),
original_correlation = Vector{Float64}(),
infered_correlation = Vector{Float64}(),);

function findingbestfit(fit_df)
    (minbic, minidx) = findmin(fit_df.bic)
    deltabic = fit_df.bic .- minbic
    bbic = findfirst(deltabic .<0.1*abs(minbic))
    #bbic = findlast(deltabic .<=10)
    if isnothing(bbic)
        bfidx = 1
    else
        bfidx = bbic
    end
    return bfidx
end




function    fitallexp!(dataidx) #dataidx =1 if ON, =2 if OFF
    datatype = ["on", "off", "correlation"]
    fit_df = DataFrame(model= Vector{String}(), datafype = Vector{String}(), bestfit = Vector{Vector{Float64}}(), timescales = Vector{Vector{Tuple}}(), fittrace = Vector{Tuple}(), res = Vector{Float64}(), bic = Vector{Float64}())# DataFrame(clone=list_clone, meanTspot=Vector{Vector{Int}}, meanTdark=0)
    
    xdata=datatofit[dataidx][1][1:findlast(datatofit[dataidx][1] .<=FRange[dataidx][2])];
    ydata=datatofit[dataidx][2][1:findlast(datatofit[dataidx][1] .<=FRange[dataidx][2])];
    xdata = xdata[ydata.>1e-2]
    ydata = log.(ydata[ydata.>1e-2])

    # modelexp 1 exp
        ubu = [10.]
        lbu = [0.]
        p0 = [.05]
        modelexp1(t, p) = -p[1] .* t
    
        fit = LsqFit.curve_fit(modelexp1, xdata, ydata, p0; lower = lbu, upper = ubu, maxIter=10000, show_trace=false, x_tol=0, g_tol=0)
        res = sum(fit.resid.^2)#sum((modelexp1(xdata,fit.param) .- ydata).^2)
        fittrace = exp.(modelexp1(xdata,fit.param))
        #bic = length(p0)*log(length(xdata)) + 2*(res)
        bic = length(p0)*log(length(xdata)) + length(xdata)*log(res./length(xdata))
    
        push!(fit_df, ("1exp", datatype[dataidx],fit.param, [(1,fit.param[1])], (xdata,fittrace),res, bic ));
    
     # modelexp 2 exp
     
        ubu = [1.; 10.; 10.]
        lbu = [0.;0.;0.]
        p0 = [.5,.05,.05]
        modelexp2(t, p) = log.(p[1] .* exp.(-p[2] .* t) .+ (1-p[1]) .* exp.(-p[3] .* t));
        
        fit = LsqFit.curve_fit(modelexp2, xdata, ydata, p0; lower = lbu, upper = ubu, maxIter=10000, show_trace=false, x_tol=0, g_tol=0)
        res = sum(fit.resid.^2)#sum((modelexp2(xdata,fit.param) .- ydata).^2)
        fittrace = exp.(modelexp2(xdata,fit.param))
        #bic = length(p0)*log(length(xdata)) + 2*(res)
        bic = length(p0)*log(length(xdata)) + length(xdata)*log(res./length(xdata))
    
        push!(fit_df, ("2exp", datatype[dataidx],fit.param,[(fit.param[1], fit.param[2]), (1-fit.param[1], fit.param[3])], (xdata,fittrace),res, bic ));
    
    
    # modelexp 3 exp 
        ubu = [1.; 1.; 10.; 10.; 10.]
        lbu = [0.;0.;0.;0.;0.]
        p0 = [.05,.05,.05,.05,.05]
        modelexp3(t, p) = log.(p[1] .* exp.(-p[3] .* t) .+ (1 .-p[1]) .* p[2] .* exp.(-p[4] .* t)  .+ (1-p[1]).*(1-p[2]) .* exp.(-p[5] .* t));
    
        fit = LsqFit.curve_fit(modelexp3, xdata, ydata, p0; lower = lbu, upper = ubu, maxIter=10000, show_trace=false, x_tol=0, g_tol=0)
        res = sum(fit.resid.^2)#sum((modelexp3(xdata,fit.param) .- ydata).^2)
        fittrace = exp.(modelexp3(xdata,fit.param))
        #bic = length(p0)*log(length(xdata)) + 2*(res)
        bic = length(p0)*log(length(xdata)) + length(xdata)*log(res./length(xdata))
    
        push!(fit_df, ("3exp", datatype[dataidx],fit.param, [(fit.param[1], fit.param[3]), ((1 .-fit.param[1]) .* fit.param[2], fit.param[4]), ((1-fit.param[1]).*(1-fit.param[2]), fit.param[5])],  (xdata,fittrace),res, bic ));
    
    # modelexp 4 exp 
        ubu = [1.; 1.; 1.; 10.; 10. ;10.; 10.]
        lbu = [0.;0.;0.;0.;0.;0.;0.]
        p0 = [.5,.5,.5,.5,.05,.05,.05]
        modelexp4(t, p) = log.(p[1] .* exp.(-p[4] .* t) .+ (1 .-p[1]) .* p[2] .* exp.(-p[5] .* t)  .+ (1-p[1]).*(1-p[2]).*p[3] .* exp.(-p[6] .* t) .+ (1-p[1]).*(1-p[2]).*(1-p[3]).* exp.(-p[7] .* t));
    
        fit = LsqFit.curve_fit(modelexp4, xdata, ydata, p0; lower = lbu, upper = ubu, maxIter=10000, show_trace=false, x_tol=0, g_tol=0)
        res = sum(fit.resid.^2)#sum((modelexp4(xdata,fit.param) .- ydata).^2)
        fittrace = exp.(modelexp4(xdata,fit.param))
        #bic = length(p0)*log(length(xdata)) + 2*(res)
        bic = length(p0)*log(length(xdata)) + length(xdata)*log(res./length(xdata))
    
        push!(fit_df, ("4exp", datatype[dataidx],fit.param,[(fit.param[1], fit.param[4]), ((1 .-fit.param[1]) .* fit.param[2], fit.param[5]), ((1-fit.param[1]).*(1-fit.param[2])*fit.param[3], fit.param[6]), ((1-fit.param[1]).*(1-fit.param[2])*(1-fit.param[3]), fit.param[7])], (xdata,fittrace),res, bic ));

        bfidx = findingbestfit(fit_df)

        return (fit_df.timescales[bfidx])
end

#simulation based time scales
sim_ts_on = fitallexp!(1)
sim_ts_off = fitallexp!(2)
#theoretical time scales
top_thy_ts_on, top_thy_ts_off = top_timescales(model_th, copy(data_table.parameters[simuid])) 

push!(fit_df, (data_table.uid[simuid],  data_table.parameters[simuid], data_table.pmiss[simuid], data_table.config[simuid], top_thy_ts_off, top_thy_ts_on, sim_ts_off, sim_ts_on, correlation_interburst_th, correlation_interburst_sim ))

save(joinpath(@__DIR__,"./bestfits/$fitname/exp_$(fitname)_$(simuid).jld2"), "fit_df", fit_df)


