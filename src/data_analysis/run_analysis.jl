using CSV
using Statistics
using Survival
using Tables, FileIO, JLD2, DataFrames
using Arrow

include("survival_analysis.jl")
using .survival_analysis

pathtodata = "/tungstenfs/scratch/ggiorget/Jana/transcrip_dynamic/E10_Mobi/live_imaging/data/306KI-ddCTCF-dpuro-MS2-HaloMCP-E10Mobi_JF549_30s_fiveclone_combined_threshold2800_curated.csv"
readout_colname = :spotdetected_filtered_curated


#load the data
df = CSV.read(pathtodata, DataFrame)
list_clone = unique(df.clone);
#remove missing 
deleteat!(df,ismissing.(df[:,readout_colname]))

#Information about the clones
df_clones = CSV.read("Positions_plate5.6.csv", DataFrame);
delete!(df_clones, ismissing.(df_clones.group_name_seq))


# build a dataframe of trakcs
track_gdf = groupby(df, :unique_id)
alltracks = DataFrame(clone = Vector{String}(), bursttrack = Vector{Vector{Bool}}(), frametrack = Vector{Vector{Int64}}(),intensitytrack = Vector{Vector{Float64}}(), trackid = Vector{Int64}())
for group in track_gdf
    push!(alltracks, (group.clone[1], group[:,readout_colname], group.frame,  group.corr_trace, group.unique_id[1]))
end

#select a frame window and cut the tracks accordingly
leftlimit = 120
rightlimit = 600
cuttracks = DataFrame(clone = Vector{String}(), bursttrack = Vector{Vector{Bool}}(), frametrack = Vector{Vector{Int64}}(),intensitytrack = Vector{Vector{Float64}}(), trackid = Vector{Int64}())
for row in eachrow(alltracks)
    bursttrack = row.bursttrack
    frametrack = row.frametrack
    intensitytrack =  row.intensitytrack
    if frametrack[end] > leftlimit && frametrack[1]< rightlimit
        if frametrack[1] < leftlimit
            deleteat!(bursttrack, (1:findfirst(frametrack .== leftlimit)-1))
            deleteat!(intensitytrack, (1:findfirst(frametrack .== leftlimit)-1))
            deleteat!(frametrack, (1:findfirst(frametrack .== leftlimit)-1))
        end
        if frametrack[end] > rightlimit
            deleteat!(bursttrack, (findfirst(frametrack .== rightlimit)+1:length(frametrack)))
            deleteat!(intensitytrack, (findfirst(frametrack .== rightlimit)+1:length(frametrack)))   
            deleteat!(frametrack, (findfirst(frametrack .== rightlimit)+1:length(frametrack)))
        end  
        push!(cuttracks, (row.clone, bursttrack, frametrack, intensitytrack, row.trackid))    
    end
end


bstat = survival_analysis.survival(cuttracks, list_clone, minimum([600,rightlimit-leftlimit]));


#probability to detect a spot
bstat.pburst = missings(Vector{Float64}, nrow(bstat))
track_gdfclone = groupby(df, :clone)
for (keyclone, subdfclone) in pairs(track_gdfclone)
    track_gdftemp = groupby(subdfclone, :frame)
    p = Vector{Float64}(undef,size(track_gdftemp,1))
    idx = 1
    for (key, subdf) in pairs(track_gdftemp)
        p[idx] = sum(subdf[:,readout_colname])/length(subdf[:,readout_colname])
        idx +=1
    end
    bstat[findfirst(item -> item == subdfclone.clone[1], bstat.clone),:pburst] = p
end


# add the distance and contact probability for each clone
bstat.distance = zeros(size(bstat,1));
bstat.absdistance = zeros(size(bstat,1));
bstat.cp = zeros(size(bstat,1));

for clone in list_clone
    cloneinfo = df_clones[occursin.(clone, df_clones.clone),:];
    if size(cloneinfo,1)>1
        println("two clones found")
    else
        bstat[bstat.clone .== clone, "distance"] = cloneinfo.distance;
        bstat[bstat.clone .== clone, "absdistance"] = abs.(cloneinfo.distance);
        bstat[bstat.clone .== clone, "cp"] = cloneinfo.cp;
    end
end
sort!(bstat,"absdistance");


Arrow.write("../../processed_data/data_$(leftlimit)_$(rightlimit).arrow", bstat)
save("../../processed_data/metadata_$(leftlimit)_$(rightlimit).jld2", "leftlimit", leftlimit, "rightlimit", rightlimit)

