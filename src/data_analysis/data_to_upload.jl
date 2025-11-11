using CSV
using Statistics
using Survival
using Tables, FileIO, JLD2, DataFrames
using Arrow

pathtodata = "/tungstenfs/scratch/ggiorget/Jana/transcrip_dynamic/E10_Mobi/live_imaging/data/306KI-ddCTCF-dpuro-MS2-HaloMCP-E10Mobi_JF549_30s_fiveclone_combined_threshold2800_curated.csv"


#load the data
df = CSV.read(pathtodata, DataFrame)
#remove missing 
deleteat!(df,ismissing.(df[:,readout_colname]))


df_toupload = select(df,:track_id, :unique_id, :clone,:frame, :x, :y, :mean_spot, :mean_spot_bleach, :corr_trace,:mean_localbackground, :mean_localbackground_bleach, :mean_completecell, :spotdetected_filtered_curated)

rename!(df_toupload, [:corr_trace => :corrected_intensity, :spotdetected_filtered_curated => :spotdetected])

CSV.write("/tungstenfs/scratch/ggiorget/Gregory/julia/BurstingDynamicsInference/processed_data/dataToupload.csv", df_toupload)