using Plots
using Tables, FileIO, JLD2
using StoThyLiveCell

include("config.jl")

include("../../../utility/dataToFit.jl")
include("../../../utility/dataToPrint.jl")

colorlist = ["#e74c3c", "#a569bd", "#5dade2", "#52be80", "#eb984e"]

cd("/tungstenfs/scratch/ggiorget/Gregory/julia/FitStoThyLiveCell/src/models/m_2s2r/run_allcl_20250218_11_p1/")

clidx = 2
bestfitTrialID = "fcfca7c6_p"*"_cl$(cloneid[clidx])"
pathToBestfit = pwd()*"/bestfits/$(bestfitTrialID).jld2"

#upload the bf model ouptuts 
bf = load(pathToBestfit) ;
bf_parameters = bf["bfparameters"];
bf_objfunc = bf["minobj"]
estimate_signal = bf["estimate_signal"]
solstats = bf["solstats"]

#upload the data
data_all = dataToFit(dataname, cloneid[clidx], tmax_nextburst, tmax_intensity)
(survival_interburst_print, survival_interburst_print_psd, survival_interburst_print_msd, survival_burst_print, survival_nextburst_print, intensityburst_print, mean_nascentrna_print, prob_burst_print, correlation_interburst_print, mrna_distribution_print) = dataToPrint(dataname, cloneid[clidx])
data_signal = (survival_burst = data_all[1], survival_interburst = data_all[2], survival_nextburst = data_all[3], prob_burst =data_all[6], mean_nascentrna = data_all[5], correlation_interburst = data_all[7], intensity_burst = data_all[4], mrna_distribution = data_all[8])
    

#print bestfits
println("bestfit parameters:")
println(bf_parameters)
println("Objective Function")
println(bf_objfunc)

#plots stats from the optimization process
p_optim_objfcts = histogram(solstats.objfct, color=:gray, legend=false)#xlims!(-5, 5)
xlabel!("Bestfit Objective Fct Eval")
ylabel!("Fraction of optim")
display(p_optim_objfcts)

p_optim_nbeval = histogram(solstats.nbfeval, color=:gray, legend=false)#xlims!(-5, 5)
xlabel!("Nb of Eval")
ylabel!("Fraction of optim")
display(p_optim_nbeval)

#calculate the model prediction

#plots
#= p_mrna_distribution = histogram(data, label="Experimental", bins=0:2:maximum(data_signal.mrna_distribution)+1, normalize=:pdf, color=:gray)
plot!(0:1:length(estimate_signal.mrna_distribution),estimate_signal.mrna_distribution, label="Fit", lw=3, color=:red)
#xlims!(-5, 5)
#ylims!(0, 0.4)
#title!("")
xlabel!("Number of mRNA")
ylabel!("Fraction of Cells")
display(p_mrna_distribution) =#

#surivival burst 
xdata = data_signal.survival_burst[1]
ydata = data_signal.survival_burst[2]
xmodel = 1:length(estimate_signal.survival_burst)
ymodel = estimate_signal.survival_burst
p_survival_burst = plot(xdata[ydata .> 0],ydata[ydata .>0], label="Data", lw=3, color=colorlist[cloneid[clidx]], yscale=:log)
plot!(xmodel[ymodel .> 0],ymodel[ymodel .>0], label="Fit", lw=3, color=:black)
xlabel!("Time (frame)")
ylabel!("Survival Prob. Burst Duration")

#surivival inter-burst 
xdata = data_signal.survival_interburst[1]
ydata = data_signal.survival_interburst[2]
xmodel = 1:length(estimate_signal.survival_interburst)
ymodel = estimate_signal.survival_interburst
p_survival_interburst = plot(xdata[ydata .> 0],ydata[ydata .>0], label="Data", lw=3, color=colorlist[cloneid[clidx]], yscale=:log)
plot!(xmodel[ymodel .> 0],ymodel[ymodel .>0], label="Fit", lw=3, color=:black)
plot!(xdata, survival_interburst_print_msd[:,2], fillrange = survival_interburst_print_psd[:,2], fillalpha = 0.35, c = colorlist[cloneid[clidx]], label = "Confidence band")
xlabel!("Time (frame)")
ylabel!("Survival Prob. Interburst Duration")

#surivival inter-burst ZOOM IN
xdata = data_signal.survival_interburst[1]
ydata = data_signal.survival_interburst[2]
xmodel = 1:length(estimate_signal.survival_interburst)
ymodel = estimate_signal.survival_interburst
p_survival_interburst_zoom = plot(xdata[ydata .> 0],ydata[ydata .>0], label="Data", lw=3, color=colorlist[cloneid[clidx]], yscale=:log, xlimits=(0,20),ylimits=(0.5,1))
plot!(xmodel[ymodel .> 0],ymodel[ymodel .>0], label="Fit", lw=3, color=:black)
plot!(xdata, survival_interburst_print_msd[:,2], fillrange = survival_interburst_print_psd[:,2], fillalpha = 0.35, c = colorlist[cloneid[clidx]], label = "Confidence band")
xlabel!("Time (frame)")
ylabel!("Survival Prob. Interburst Duration")

#surivival next_burst
xdata = data_signal.survival_nextburst[1]
ydata = data_signal.survival_nextburst[2]
xmodel = 1:length(estimate_signal.survival_nextburst)
ymodel = estimate_signal.survival_nextburst
p_survival_nextburst = plot(xdata[ydata .> 0],ydata[ydata .>0], label="Data", lw=3, color=colorlist[cloneid[clidx]], yscale=:log)
plot!(xmodel[ymodel .> 0],ymodel[ymodel .>0], label="Fit", lw=3, color=:black)
xlabel!("Time (frame)")
ylabel!("Survival Prob. Next Burst Duration")

#mean nascent rna
xdata = "data"
ydata = data_signal.mean_nascentrna
xmodel = "fit"
ymodel = estimate_signal.mean_nascentrna
p_mean_nascentrna = plot(bar([xdata,xmodel], [ydata,ymodel], legend =false, color=[colorlist[cloneid[clidx]],:grey]))
xlabel!("Time (frame)")
ylabel!("Mean # Nascent mRNA")

#proba burst
xdata = "data"
ydata = data_signal.prob_burst
xmodel = "fit"
ymodel = estimate_signal.prob_burst
p_prob_burst = plot(bar([xdata,xmodel], [ydata,ymodel], legend =false, color=[colorlist[cloneid[clidx]],:grey]))
xlabel!("Time (frame)")
ylabel!("Prob. Detect Burst")

#correlation inter burst
xdata = "data"
ydata = data_signal.correlation_interburst
xmodel = "fit"
ymodel = estimate_signal.correlation_interburst
p_correlation_interburst = plot(bar([xdata,xmodel], [ydata,ymodel], legend =false, color=[colorlist[cloneid[clidx]],:grey]))
xlabel!("Time (frame)")
ylabel!("Correlation Interburst Duration")


#intensity burst
xdata = data_signal.intensity_burst[1]
ydata = data_signal.intensity_burst[2]
xmodel = 1:length(estimate_signal.intensity_burst)
ymodel = estimate_signal.intensity_burst
p_intensity_burst = plot(xdata[ydata .> 0],ydata[ydata .>0], label="Data", lw=3, color=colorlist[cloneid[clidx]])
plot!(xmodel[ymodel .> 0],ymodel[ymodel .>0], label="Fit", lw=3, color=:black)
xlabel!("Time (frame)")
ylabel!("Normalized Intensity")


display(p_survival_burst)
display(p_survival_interburst)
display(p_survival_interburst_zoom)
display(p_survival_nextburst)
display(p_mean_nascentrna)
display(p_prob_burst)
display(p_correlation_interburst)
display(p_intensity_burst)