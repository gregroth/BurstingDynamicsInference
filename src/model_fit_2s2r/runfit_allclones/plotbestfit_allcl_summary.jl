using Plots
using StatsPlots
using Tables, FileIO, JLD2
using DataFrames
using StoThyLiveCell
using XLSX

include("config.jl")

include("../../../utility/dataToFit.jl")
include("../../../utility/dataToPrint.jl")

colorlist = ["#e74c3c", "#a569bd", "#5dade2", "#52be80", "#eb984e"]
colorlist_fit = ["#950606", "#9400D3", "#0000ff", "#05472A"]
colorbarplots =["#E5AA70" "#898989"]

pltsize = 160
labelfontsize = 6

cd("/tungstenfs/scratch/ggiorget/Gregory/julia/FitStoThyLiveCell/src/models/m_2s2r/run_allcl_20250218_11_p1/")

printid = 1
foldername = "2s2r_model_allclones_20250219_directStrategy_p1"


estimate_signal = []
bf_parameters = []


#upload summary of data
datasummary = load("/tungstenfs/scratch/ggiorget/Gregory/julia/jana_bursting/data/$(dataname).jld2") ;
bstat = datasummary["bstat"];

bs = load("/tungstenfs/scratch/ggiorget/Gregory/julia/jana_bursting/data/memory_bootstrap_summary.jld2") ;
memory_sy = bs["memory_sy"]

bs_pb = load("/tungstenfs/scratch/ggiorget/Gregory/julia/jana_bursting/data/pburst_bootstrap.jld2") ;
pburst_sy = bs_pb["pburst_sy"]

#best fit parameters for the reference clone
#pathToBestfitRef = "/tungstenfs/scratch/ggiorget/Gregory/julia/FitStoThyLiveCell/src/models/m_2s2r/run_allcl_20250218_11_p1/$(bestfitRefTrialID).jld2"
#bf_ref = load(pathToBestfitRef) ;
push!(bf_parameters, bf_ref["bfparameters"][1:7]);
push!(estimate_signal, bf_ref["estimate_signal"])

#best fit parameters for the other clones
for clidx = 2:4
    bestfitTrialID = "fcfca7c6_p"*"_cl$clidx"
    pathToBestfit = pwd()*"/bestfits/$(bestfitTrialID).jld2"
    bf= load(pathToBestfit) ;
    push!(bf_parameters, bf["bfparameters"][1:7]);
    push!(estimate_signal, bf["estimate_signal"])
end 

#upload the data
function uploaddata()
    survival_interburst_print_psd = []
    survival_interburst_print_msd = []
    data_signal = []
    for clidx = 1:4
        data_all = dataToFit(dataname, clidx, tmax_nextburst, tmax_intensity)
        (survival_interburst_print, survival_interburst_print_psd_temp, survival_interburst_print_msd_temp, survival_burst_print, survival_nextburst_print, intensityburst_print, mean_nascentrna_print, prob_burst_print, correlation_interburst_print, mrna_distribution_print) = dataToPrint(dataname, clidx)
        data_signal_temp = (survival_burst = data_all[1], survival_interburst = data_all[2], survival_nextburst = data_all[3], prob_burst =data_all[6], mean_nascentrna = data_all[5], correlation_interburst = data_all[7], intensity_burst = data_all[4], mrna_distribution = data_all[8])
        push!(survival_interburst_print_psd,survival_interburst_print_psd_temp)
        push!(survival_interburst_print_msd,survival_interburst_print_msd_temp)
        push!(data_signal,data_signal_temp)
    end
    return (data_signal,survival_interburst_print_psd,survival_interburst_print_msd)
end

if !@isdefined(data_signal)
    (data_signal,survival_interburst_print_psd,survival_interburst_print_msd) = uploaddata()
elseif !(typeof(data_signal) == Vector{Any})
    (data_signal,survival_interburst_print_psd,survival_interburst_print_msd) = uploaddata()
end


#surivival burst 
p_survival_burst = plot(yscale=:log,guidefontsize=labelfontsize, tickfontsize=labelfontsize, legend=false)
for clidx = 1:4
    xdata = data_signal[clidx].survival_burst[1]
    ydata = data_signal[clidx].survival_burst[2]
    xmodel = 1:length(estimate_signal[clidx].survival_burst)
    ymodel = estimate_signal[clidx].survival_burst
    plot!(xdata[ydata .> 0],ydata[ydata .>0], label="$(abs(round(Int,bstat.distance[clidx]/1000)))kb", lw=3, color=colorlist[clidx])
    plot!(xmodel[ymodel .> 0],ymodel[ymodel .>0], lw=2, label="Fit", color=colorlist_fit[clidx], linestyle=:dash)
end
xlabel!("Time (frame)")
ylabel!("Survival Probability")
#plot!(legend=:bottomleft)
plot!(size=(pltsize,pltsize))

#surivival inter-burst 
p_survival_interburst = plot(yscale=:log,guidefontsize=labelfontsize, tickfontsize=labelfontsize, legend=false)
for clidx = 1:4
    xdata = data_signal[clidx].survival_interburst[1]
    ydata = data_signal[clidx].survival_interburst[2]
    xmodel = 1:length(estimate_signal[clidx].survival_interburst)
    ymodel = estimate_signal[clidx].survival_interburst
    plot!(xdata[ydata .> 0],ydata[ydata .>0], label="$(abs(round(Int,bstat.distance[clidx]/1000)))kb", lw=3, color=colorlist[clidx])
    plot!(xmodel[ymodel .> 0],ymodel[ymodel .>0], lw=3, label="Fit",color=colorlist_fit[clidx], linestyle=:dash)
    plot!(xdata, survival_interburst_print_msd[clidx][:,2], fillrange = survival_interburst_print_psd[clidx][:,2], fillalpha = 0.35, c = colorlist[clidx], label = "Confidence band")
end
xlabel!("Time (frame)")
ylabel!("Survival Probability")
#plot!(legend=:bottomleft)
plot!(size=(pltsize,pltsize))

#surivival inter-burst ZOOM IN
clidx = 1
xdata = data_signal[clidx].survival_interburst[1]
ydata = data_signal[clidx].survival_interburst[2]
xmodel = 1:length(estimate_signal[clidx].survival_interburst)
ymodel = estimate_signal[clidx].survival_interburst
p_survival_interburst_zoom_cl1 = plot(xdata[ydata .> 0],ydata[ydata .>0], label="$(abs(round(Int,bstat.distance[clidx]/1000)))kb", lw=3, color=colorlist[clidx],yscale=:log, xlimits=(0,80),ylimits=(0.2,1),legend=false)
plot!(xmodel[ymodel .> 0],ymodel[ymodel .>0], label="Fit", lw=3, color=colorlist_fit[clidx], linestyle=:dash)
plot!(xdata, survival_interburst_print_msd[clidx][:,2], fillrange = survival_interburst_print_psd[clidx][:,2], fillalpha = 0.35, c = colorlist[clidx], label = "Confidence band")
xlabel!("Time (frame)")
ylabel!("Survival Probability")
#plot!(legend=:bottomleft)
plot!(size=(pltsize,pltsize))

clidx = 2
xdata = data_signal[clidx].survival_interburst[1]
ydata = data_signal[clidx].survival_interburst[2]
xmodel = 1:length(estimate_signal[clidx].survival_interburst)
ymodel = estimate_signal[clidx].survival_interburst
p_survival_interburst_zoom_cl2 = plot(xdata[ydata .> 0],ydata[ydata .>0], label="$(abs(round(Int,bstat.distance[clidx]/1000)))kb", lw=3, color=colorlist[clidx],yscale=:log, xlimits=(0,80),ylimits=(0.2,1),legend=false)
plot!(xmodel[ymodel .> 0],ymodel[ymodel .>0], label="Fit", lw=3, color=colorlist_fit[clidx], linestyle=:dash)
plot!(xdata, survival_interburst_print_msd[clidx][:,2], fillrange = survival_interburst_print_psd[clidx][:,2], fillalpha = 0.35, c = colorlist[clidx], label = "Confidence band")
xlabel!("Time (frame)")
ylabel!("Survival Probability")
#plot!(legend=:bottomleft)
plot!(size=(pltsize,pltsize))

clidx = 3
xdata = data_signal[clidx].survival_interburst[1]
ydata = data_signal[clidx].survival_interburst[2]
xmodel = 1:length(estimate_signal[clidx].survival_interburst)
ymodel = estimate_signal[clidx].survival_interburst
p_survival_interburst_zoom_cl3 = plot(xdata[ydata .> 0],ydata[ydata .>0], label="$(abs(round(Int,bstat.distance[clidx]/1000)))kb", lw=3, color=colorlist[clidx],yscale=:log, xlimits=(0,80),ylimits=(0.2,1),legend=false)
plot!(xmodel[ymodel .> 0],ymodel[ymodel .>0], label="Fit", lw=3, color=colorlist_fit[clidx], linestyle=:dash)
plot!(xdata, survival_interburst_print_msd[clidx][:,2], fillrange = survival_interburst_print_psd[clidx][:,2], fillalpha = 0.35, c = colorlist[clidx], label = "Confidence band")
xlabel!("Time (frame)")
ylabel!("Survival Probability")
#plot!(legend=:bottomleft)
plot!(size=(pltsize,pltsize))

clidx = 4
xdata = data_signal[clidx].survival_interburst[1]
ydata = data_signal[clidx].survival_interburst[2]
xmodel = 1:length(estimate_signal[clidx].survival_interburst)
ymodel = estimate_signal[clidx].survival_interburst
p_survival_interburst_zoom_cl4 = plot(xdata[ydata .> 0],ydata[ydata .>0], label="$(abs(round(Int,bstat.distance[clidx]/1000)))kb", lw=3, color=colorlist[clidx],yscale=:log, xlimits=(0,80),ylimits=(0.2,1),legend=false)
plot!(xmodel[ymodel .> 0],ymodel[ymodel .>0], label="Fit", lw=3, color=colorlist_fit[clidx], linestyle=:dash)
plot!(xdata, survival_interburst_print_msd[clidx][:,2], fillrange = survival_interburst_print_psd[clidx][:,2], fillalpha = 0.35, c = colorlist[clidx], label = "Confidence band")
xlabel!("Time (frame)")
ylabel!("Survival Probability")
#plot!(legend=:bottomleft)
plot!(size=(pltsize,pltsize))



#surivival next_burst
p_survival_nextburst = plot(yscale=:log,guidefontsize=labelfontsize, tickfontsize=labelfontsize, legend=false)
for clidx = 1:4
    xdata = data_signal[clidx].survival_nextburst[1]
    ydata = data_signal[clidx].survival_nextburst[2]
    xmodel = 1:length(estimate_signal[clidx].survival_nextburst)
    ymodel = estimate_signal[clidx].survival_nextburst
    plot!(xdata[ydata .> 0],ydata[ydata .>0], label="$(abs(round(Int,bstat.distance[clidx]/1000)))kb", lw=3, color=colorlist[clidx])
    plot!(xmodel[ymodel .> 0],ymodel[ymodel .>0], lw=3, label="Fit",color=colorlist[clidx], linestyle=:dash)
end
xlabel!("Time (frame)")
ylabel!("Survival Probability")
#plot!(legend=:bottomleft)
plot!(size=(pltsize,pltsize))


#mean nascent rna
xdata = "data"
ydata = [data_signal[1].mean_nascentrna,data_signal[2].mean_nascentrna,data_signal[3].mean_nascentrna,data_signal[4].mean_nascentrna]
xmodel = "fit"
ymodel = [estimate_signal[1].mean_nascentrna,estimate_signal[2].mean_nascentrna,estimate_signal[3].mean_nascentrna,estimate_signal[4].mean_nascentrna]
ctg = repeat(["Data", "Model"], inner = 4)
name = repeat(["CL1-5kb", "CL2-30kb", "CL3-149kb", "CL4-233kb"], outer = 2)
p_mean_nascentrna = groupedbar(name, [ydata ymodel], group = ctg, color=colorbarplots)
xlabel!("Clones")
ylabel!("Mean # Nascent mRNA")

#proba burst
xdata = "data"
ydata = [data_signal[1].prob_burst,data_signal[2].prob_burst,data_signal[3].prob_burst,data_signal[4].prob_burst]
xmodel = "fit"
ymodel = [estimate_signal[1].prob_burst,estimate_signal[2].prob_burst,estimate_signal[3].prob_burst,estimate_signal[4].prob_burst]
ctg = repeat(["Data", "Model"], inner = 4)
name = repeat(["CL1-5kb", "CL2-30kb", "CL3-149kb", "CL4-233kb"], outer = 2)
p_prob_burst = groupedbar(name, [ydata ymodel], group = ctg, color=colorbarplots,  yerror=[pburst_sy.sd_pb zeros(4,1)], guidefontsize=labelfontsize, tickfontsize=labelfontsize, legend=false)
xlabel!("Clones")
ylabel!("Probability")
plot!(size=(pltsize,pltsize))

#correlation inter burst
xdata = "data"
ydata = [data_signal[1].correlation_interburst,data_signal[2].correlation_interburst,data_signal[3].correlation_interburst,data_signal[4].correlation_interburst]
xmodel = "fit"
ymodel = [estimate_signal[1].correlation_interburst,estimate_signal[2].correlation_interburst,estimate_signal[3].correlation_interburst,estimate_signal[4].correlation_interburst]
ctg = repeat(["Data", "Model"], inner = 4)
name = repeat(["CL1-5kb", "CL2-30kb", "CL3-149kb", "CL4-233kb"], outer = 2)
p_correlation_interburst = groupedbar(name, [ydata ymodel], group = ctg, color=colorbarplots, yerror=[memory_sy.sd zeros(4,1)],guidefontsize=labelfontsize, tickfontsize=labelfontsize, legend=false)
xlabel!("Clones")
ylabel!("Correlation")
plot!(size=(pltsize,pltsize))


display(p_survival_burst)
display(p_survival_interburst)
display(p_survival_interburst_zoom_cl1)
display(p_survival_interburst_zoom_cl2)
display(p_survival_interburst_zoom_cl3)
display(p_survival_interburst_zoom_cl4)
display(p_survival_nextburst)
display(p_mean_nascentrna)
display(p_prob_burst)
display(p_correlation_interburst)
display(p_intensity_burst)



#summary of the bestfit parameter values

df_summary_frame = DataFrame(clone = Vector{String}(), kforward = Vector{Float64}(), kbackward = Vector{Float64}(), konLowR = Vector{Float64}(), konHighR = Vector{Float64}(), koff = Vector{Float64}(), kinitiation = Vector{Float64}(), delta = Vector{Float64}())
df_summary_minute = DataFrame(clone = Vector{String}(), kforward = Vector{Float64}(), kbackward = Vector{Float64}(), konLowR = Vector{Float64}(), konHighR = Vector{Float64}(), koff = Vector{Float64}(), kinitiation = Vector{Float64}(), delta = Vector{Float64}())
df_summary_waitingtime = DataFrame(clone = Vector{String}(), TlowRegime = Vector{Float64}(), ThighRegime= Vector{Float64}(), ToffLowR = Vector{Float64}(), ToffHighR = Vector{Float64}(), Ton = Vector{Float64}(), Tinitiation = Vector{Float64}(), Tsignal = Vector{Float64}())

bf_parameters =stack(bf_parameters)
bf_parameters_min = 2 .* bf_parameters
bf_waitingtime = 1 ./ bf_parameters_min
bf_waitingtime =round.(bf_waitingtime;digits=1)
bf_parameters_min =round.(bf_parameters_min;digits=5)

for i=1:4
    push!(df_summary_frame,(bstat.clone[i], bf_parameters[1,i], bf_parameters[2,i], bf_parameters[3,i], bf_parameters[4,i], bf_parameters[5,i], bf_parameters[6,i], bf_parameters[7,i]))
    push!(df_summary_minute,(bstat.clone[i], bf_parameters_min[1,i], bf_parameters_min[2,i], bf_parameters_min[3,i], bf_parameters_min[4,i], bf_parameters_min[5,i], bf_parameters_min[6,i], bf_parameters_min[7,i]))
    push!(df_summary_waitingtime,(bstat.clone[i], bf_waitingtime[1,i], bf_waitingtime[2,i], bf_waitingtime[3,i], bf_waitingtime[4,i], bf_waitingtime[5,i], bf_waitingtime[6,i], bf_waitingtime[7,i]))
end

wprint = 250
wprintzoom = 250
    # Print and save
    if printid == 1
        pathnewfolder = "/tungstenfs/scratch/ggiorget/Gregory/julia/bursting_manuscript/figures/$foldername"
        if !ispath(pathnewfolder)
        mkdir(pathnewfolder)
        end
        savefig(p_survival_burst, "$(pathnewfolder)/fig_survival_burst.pdf")
        savefig(p_survival_interburst, "$(pathnewfolder)/fig_survival_interburst.pdf")
        savefig(p_survival_interburst_zoom_cl1, "$(pathnewfolder)/fig_interburst_zoom_cl1.pdf")
        savefig(p_survival_interburst_zoom_cl2, "$(pathnewfolder)/fig_interburst_zoom_cl2.pdf")
        savefig(p_survival_interburst_zoom_cl3, "$(pathnewfolder)/fig_interburst_zoom_cl3.pdf")
        savefig(p_survival_interburst_zoom_cl4, "$(pathnewfolder)/fig_interburst_zoom_cl4.pdf")

        savefig(p_survival_nextburst, "$(pathnewfolder)/fig_survival_nextburst.pdf")
        savefig(p_mean_nascentrna, "$(pathnewfolder)/fig_mean_nascentrna.pdf")
        savefig(p_prob_burst, "$(pathnewfolder)/fig_prob_burst.pdf")
        savefig(p_correlation_interburst, "$(pathnewfolder)/fig_correlation_interburst.pdf")
        savefig(p_intensity_burst, "$(pathnewfolder)/fig_intensity_burst.pdf")
        
        if ispath("$(pathnewfolder)/bestfitparameters_frame.xlsx")
            rm("$(pathnewfolder)/bestfitparameters_frame.xlsx")
            rm("$(pathnewfolder)/bestfitparameters_minute.xlsx")
            rm("$(pathnewfolder)/bestfitwaitingtime_minute.xlsx")
        end
        XLSX.writetable("$(pathnewfolder)/bestfitparameters_frame.xlsx", df_summary_frame)
        XLSX.writetable("$(pathnewfolder)/bestfitparameters_minute.xlsx", df_summary_minute)
        XLSX.writetable("$(pathnewfolder)/bestfitwaitingtime_minute.xlsx", df_summary_waitingtime)
    end
