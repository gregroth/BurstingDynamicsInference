
lb = [SRange[i][1] for i in eachindex(SRange)]
ub = [SRange[i][2] for i in eachindex(SRange)]
db = ub - lb



function simulation_track(tracklength::Vector{Int64},parameters)
    trajmodel_set_temp = []
    for lg = tracklength
        traj_model = simulate_bursttraj(parameters, lg, 1)
        push!(trajmodel_set_temp, traj_model)
    end
    return trajmodel_set_temp
end




