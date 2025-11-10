
fitstat = DataFrame(uid = Vector{Int64}(), parameters = Vector{Vector{Float64}}(), pmiss = Vector{Float64}(), 
config = Vector{NamedTuple{}}(),
top_thy_ts_off = Vector{Vector{Tuple}}(),
top_thy_ts_on= Vector{Vector{Tuple}}(),
sim_ts_off = Vector{Vector{Tuple}}(),
sim_ts_on= Vector{Vector{Tuple}}(),
original_correlation = Vector{Float64}(),
infered_correlation = Vector{Float64}(),);


bestfit_list = glob("*.jld2","./fitsim/bestfits/$fitname")
for i in eachindex(bestfit_list)
    bf= load(bestfit_list[i]);
    fit_df_temp = bf["fit_df"]
    push!(fitstat, fit_df_temp[1,:])
end

Arrow.write("./results/summary_$(fitname).arrow", fitstat)
