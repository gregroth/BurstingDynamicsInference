#Calculate the burst stats from a dataframe of tracks
#
#
# Gregory Roth
# First version:            25.04.2023
# Modified:                 25.04.2023
#
#-----------------------------------------------------------------------------------------------

module survival_analysis

using DataFrames
using Survival
using StatsBase

export survival

"""
    event_wcensored(track)
    call for events (sequence of ones)
TBW
"""
function    event_wcensored(track)
    #check if the track is constant
    eventtime = Vector{Int64}();
    eventcensor = Vector{Bool}();
    startevent = false;
    if minimum(track) == 1
        push!(eventtime, length(track)) 
        push!(eventcensor, false)
    else
        #check if the first event is left censored
        if track[1] == 1
            startevent = true
        end
        lastnum=0;
        currentnum=0
        runlength = 0
        starttemp = 0
        tracklengthtemp = length(track) 
        for i in axes(track,1)
            currentnum=track[i]
            #println("time is $i and current is $currentnum and last is $lastnum")
            if lastnum == 0 && currentnum == 1
                runlength = 1
                starttemp = i
            elseif lastnum == currentnum && currentnum == 1
                runlength = runlength + 1
            elseif lastnum == 1 && currentnum == 0
                push!(eventtime, runlength) 
                push!(eventcensor, true)
            end
            lastnum = currentnum
        end
        if track[end] == 1
            push!(eventtime,runlength)
            push!(eventcensor,false)
        end 
        if startevent #remove the first event if it is left censored
            deleteat!(eventtime,1)
            deleteat!(eventcensor,1)
        end
    end
    return eventtime, eventcensor
end

"""
    survival(df_track, list_clone, colname)
    calculate the survival curves for the burst and waiting times
TBW
"""
function survival(df_track, list_clone, lenggthfulltrack)
    # call events
    bstat = DataFrame(clone = Vector{String}(), burst_eventlist = Vector{Vector{Float64}}(),burst_censorlist = Vector{Vector{Int32}}(), burst_events = Vector{Vector{EventTime}}(), km_burst = Vector{KaplanMeier{Float64, Int64}}(), 
        interburst_eventlist = Vector{Vector{Float64}}(), interburst_censorlist = Vector{Vector{Int32}}(), interburst_events = Vector{Vector{EventTime}}(), km_interburst = Vector{KaplanMeier{Float64, Int64}}(),
                pfullinterburst = Vector{Float32}(), nbfulltrack = Vector{Float32}(), meanburstfreq = Vector{Float64}(), cvburstfreq = Vector{Float64}(),
                    meannbburst = Vector{Float64}(), cvnbburst = Vector{Float64}(), burst_tracklength = Vector{Vector{Float64}}(), interburst_tracklength = Vector{Vector{Float64}}(),nbfullon = Vector{Float32}(),
                       pburststar =  Vector{Float32}(), eventcorr_burst = Vector{Array{Int64}}(), eventcorr_interburst = Vector{Array{Int64}}(),nb_track = Vector{Float64}(), nb_interbursttrack = Vector{Float64}(),burstfreq = Vector{Vector{Float64}}()  );# DataFrame(clone=list_clone, meanTburst=Vector{Vector{Int}}, meanTinterburst=0)
                          

    for clone in list_clone
        eventtime_burst = Vector{Int64}[];
        eventcorr_burst = Array{Int64}(undef,1,3)
        tracklength_burst = Vector{Int64}[];
        eventcensor_burst = Vector{Bool}[];
    
        eventtime_interburst = Vector{Int64}[];
        tracklength_interburst = Vector{Int64}[];
        eventcensor_interburst = Vector{Bool}[];
        eventcorr_interburst = Array{Int64}(undef,1,3)

        track_selectedclone = df_track[df_track.clone .== clone, :];
        nb_fulldark = 0;
        nb_fulltrack = 0;
        nb_fullbursting = 0;


        burstfreq = zeros(0);
        nbburst = zeros(0);
        
        nbtimepts = 0
        nbtimeburst = 0
        
        nbtrack = 0
        nbinterbursttrack = 0
        for row in eachrow(track_selectedclone)
            nbtrack = nbtrack+1
            bursttemp = row.bursttrack;
            nbtimepts += length(bursttemp)
            nbtimeburst += sum((bursttemp.>0))
      
            eventtime_bursttemp, eventcensor_bursttemp = event_wcensored(bursttemp)
      
            if maximum(bursttemp) == 1 && minimum(bursttemp) == 0 # consider only track with at least one burst
                # BURST
                eventtime_burst = vcat(eventtime_burst, eventtime_bursttemp);
                eventcensor_burst = vcat(eventcensor_burst, eventcensor_bursttemp);
                # INTERBURST
                interbursttemp = .-bursttemp.+1;
                eventtime_interbursttemp, eventcensor_interbursttemp = event_wcensored(interbursttemp)
                eventtime_interburst = vcat(eventtime_interburst, eventtime_interbursttemp);
                eventcensor_interburst = vcat(eventcensor_interburst, eventcensor_interbursttemp);

                # Extract consecutive events
                if length(eventcensor_interbursttemp)>=2
                    if eventcensor_interbursttemp[end] == false
                        if length(eventcensor_interbursttemp)>=3
                            eventcorr_interburst_temp = [eventtime_interbursttemp[1:end-2] eventtime_interbursttemp[2:end-1] row.trackid .*ones(length(eventtime_interbursttemp)-2)]
                            eventcorr_interburst = vcat(eventcorr_interburst,eventcorr_interburst_temp)
                        end
                    else
                        eventcorr_interburst_temp = [eventtime_interbursttemp[1:end-1] eventtime_interbursttemp[2:end] row.trackid .*ones(length(eventtime_interbursttemp)-1)]
                        eventcorr_interburst = vcat(eventcorr_interburst,eventcorr_interburst_temp)
                    end
                    
                end
                if length(eventcensor_bursttemp)>=2
                    if eventcensor_bursttemp[end] == false
                        if length(eventcensor_bursttemp)>=3
                            eventcorr_burst_temp = [eventtime_bursttemp[1:end-2] eventtime_bursttemp[2:end-1] row.trackid .*ones(length(eventtime_bursttemp)-2)]
                            eventcorr_burst = vcat(eventcorr_burst,eventcorr_burst_temp)
                        end
                    else
                        eventcorr_burst_temp = [eventtime_bursttemp[1:end-1] eventtime_bursttemp[2:end] row.trackid .*ones(length(eventtime_bursttemp)-1)]
                        eventcorr_burst = vcat(eventcorr_burst,eventcorr_burst_temp)
                    end
                end
            else
                nbinterbursttrack = nbinterbursttrack +1
            end

            # counting the full-length track
            if length(bursttemp) >= lenggthfulltrack
                nb_fulltrack += 1;
                push!(burstfreq, length(eventtime_bursttemp)/length(bursttemp))   
                push!(nbburst, length(eventtime_bursttemp))   
                if maximum(bursttemp) == 0
                    nb_fulldark += 1;
                else
                    if minimum(bursttemp) == 1
                        nb_fullbursting += 1
                    end
                end
            end
        end
        if isempty(eventtime_burst)
            eventtime_burst = [0]
            eventcensor_burst = [0]
        end
        event_burst = EventTime.(eventtime_burst, eventcensor_burst);
        
        if isempty(eventtime_interburst)
            eventtime_interburst = [0]
            eventcensor_interburst = [0]
        end
        event_interburst = EventTime.(eventtime_interburst, eventcensor_interburst);
        

        prop_fullinterburst = nb_fulldark/nb_fulltrack;
        meanburstfreq = mean(burstfreq)
        cvburstfreq = std(burstfreq)/meanburstfreq
        meannbburst = mean(nbburst)
        cvnbburst = std(nbburst)/meannbburst
        # Kaplan Meier analysis
        km_burst = fit(KaplanMeier, event_burst);
        km_interburst = fit(KaplanMeier, event_interburst);
        push!(bstat, (clone, eventtime_burst, eventcensor_burst, event_burst, km_burst, 
        eventtime_interburst, eventcensor_interburst, event_interburst, km_interburst, 
                prop_fullinterburst, nb_fulltrack, meanburstfreq, cvburstfreq, 
                    meannbburst, cvnbburst, tracklength_burst, tracklength_interburst, nb_fullbursting, 
                        nbtimeburst/nbtimepts, eventcorr_burst[2:end,:], eventcorr_interburst[2:end,:], nbtrack, nbinterbursttrack, burstfreq));
    end
    return bstat
end



"""
    survival_wosinglet_fast(trajset, treshold)

    survival analysis + removing of the singlet, for set of trajectories that have the same length. It is mainly made for simulation based data.
    trajset = Vector of Vector = set of trajectories
    threshold = integer = rna number threshold above which a burst is called 
"""
function survival_wosinglet(trajset, treshold)
    # call events
    bstat = DataFrame(km_burst = Vector{KaplanMeier{Float64, Int64}}(), km_interburst = Vector{KaplanMeier{Float64, Int64}}(), 
                            meanburstfreq = Vector{Float64}(), pburststar =  Vector{Float64}(), 
                                eventcorr_burst = Vector{Array{Int64}}(), eventcorr_interburst = Vector{Array{Int64}}());# DataFrame(clone=list_clone, meanTspot=Vector{Vector{Int}}, meanTdark=0)

        eventtime_spot = Vector{Int64}[];
        eventcorr_burst = Array{Int64}(undef,1,2)
        eventcorr_interburst = Array{Int64}(undef,1,2)
        eventcensor_spot = Vector{Bool}[];

        eventtime_dark = Vector{Int64}[];
        eventcensor_dark = Vector{Bool}[];
        burstfreq = zeros(0);
  
        nbtimepts = 0
        nbtimespot = 0
        for k in eachindex(trajset)
            spottemp = (trajset[k].>=treshold);
            #REMOVING SINGLETS
            for i in eachindex(spottemp)
                if i==1
                    if (spottemp[i]==1) && (spottemp[i+1]==0)
                        spottemp[i] = 0
                    end
                elseif i==length(spottemp)
                    if (spottemp[i]==1) && (spottemp[i-1]==0)
                        spottemp[i] = 0
                    end
                else
                    if (spottemp[i]==1) && (spottemp[i-1]==0) && (spottemp[i+1]==0)
                        spottemp[i] = 0
                    end
                end
            end
            nbtimepts += length(spottemp)
            nbtimespot += sum((spottemp.>0))
  
            eventtime_spottemp, eventcensor_spottemp  = event_wcensored(spottemp)
            
            if maximum(spottemp) == 1 && minimum(spottemp) == 0 # consider only track with at least one burst
                # ON event
                eventtime_spot = vcat(eventtime_spot, eventtime_spottemp);
                eventcensor_spot = vcat(eventcensor_spot, eventcensor_spottemp);
                # OFF events
                darktemp = .-spottemp.+1;
                eventtime_darktemp, eventcensor_darktemp = event_wcensored(darktemp)
                if length(eventcensor_darktemp)>=2
                        eventcorr_interburst_temp = [eventtime_darktemp[1:end-2] eventtime_darktemp[2:end-1]]
                        eventcorr_interburst = vcat(eventcorr_interburst,eventcorr_interburst_temp) 
                end
                if length(eventcensor_spottemp)>=2
                        eventcorr_burst_temp = [eventtime_spottemp[1:end-2] eventtime_spottemp[2:end-1]]
                        eventcorr_burst = vcat(eventcorr_burst,eventcorr_burst_temp)
                end
                eventtime_dark = vcat(eventtime_dark, eventtime_darktemp);
                eventcensor_dark = vcat(eventcensor_dark, eventcensor_darktemp);
            end
                push!(burstfreq, length(eventtime_spottemp)/length(spottemp))     
        end
        if isempty(eventtime_spot)
            eventtime_spot = [0]
            eventcensor_spot = [0]
        end
        event_spot = EventTime.(eventtime_spot, eventcensor_spot);
        
        if isempty(eventtime_dark)
            eventtime_dark = [0]
            eventcensor_dark = [0]
        end
        event_dark = EventTime.(eventtime_dark, eventcensor_dark);
        
        meanburstfreq = mean(burstfreq)
        # Kaplan Meier analysis
        km_burst = fit(KaplanMeier, event_spot);
        km_interburst = fit(KaplanMeier, event_dark);
        push!(bstat, (km_burst, km_interburst,  
                        meanburstfreq, nbtimespot/nbtimepts, 
                            eventcorr_burst[2:end,:], eventcorr_interburst[2:end,:]));  
        return bstat
end

end

