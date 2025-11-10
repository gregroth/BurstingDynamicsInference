using Random
using Distributions


function intensity(c,M)
    return c*M
end

function noiseSDT(σ_b, β, I)
    return sqrt(σ_b^2 + β*I)
end


function runlengths(v::BitVector)
    # find transitions between 0 and 1
    diffv = diff([0; v; 0])  # pad with zeros
    starts = findall(==(1), diffv)       # indices where 0 → 1
    ends   = findall(==(-1), diffv)      # indices where 1 → 0
    return starts, ends .- starts
end

function runlengths(v::Vector{Int64})
    # find transitions between 0 and 1
    diffv = diff([0; v; 0])  # pad with zeros
    starts = findall(==(1), diffv)       # indices where 0 → 1
    ends   = findall(==(-1), diffv)      # indices where 1 → 0
    return starts, ends .- starts
end

#typical parameters
#   σ_b = 1.22
#   β = 0.16
#   c = 500
function noise(traj, σ_b, β, c)
    trajnoise = zeros(length(traj))
    for i in eachindex(traj)
        trajnoise[i] = rand(Normal(intensity(c,traj[i]), noiseSDT(σ_b, β, intensity(c,traj[i]))))
    end
    return trajnoise
end


function missedrna(traj, p, minnbrna)
    trajnoise = copy(traj)
    if p>0
        tidx =1 
        nextcanbe = false
        while(tidx<=length(traj))
            if traj[tidx] ∈ (1:minnbrna)
                k = 0
                if tidx < length(traj)
                    nextcanbe = true
                else
                    nextcanbe = false
                end
                while(nextcanbe)
                    if traj[tidx+k+1] ∈ (1:minnbrna)
                        k += 1
                        if tidx + k < length(traj)
                            nextcanbe = true
                        else
                            nextcanbe = false
                        end
                    else
                        nextcanbe = false
                    end
                end

                if rand() <= p
                    trajnoise[tidx:tidx+k] .= 0
                end
                tidx += k+1
            else
                tidx += 1
            end            
        end
    end
    return trajnoise
end

#= function missedrna_fast(traj, p, minnbrna)
    trajnoise = copy(traj)
    if p>0
        traj_lowburst = (traj .<=minnbrna)
        starts, durations = runlengths(traj_lowburst)
        for i in eachindex(starts)
            if rand().<= p
                trajnoise[starts[i]:starts[i]+durations[i]-1] .= 0
            end
        end
    end
    return trajnoise
end =#

function missedrnaperframe(traj, p)
    trajnoise = copy(traj)
    if p>0
        for i in eachindex(traj)
            if traj[i] == 1
                if rand() <= p
                    trajnoise[i] = 0
                end
            end           
        end
    end
    return trajnoise
end

function missedanyrnaperframe(traj, p)
    trajnoise = copy(traj)
    if p>0
        for i in eachindex(traj)
            trajnoise[i] = sum(rand(Int(traj[i])) .> p)
        end
    end
    return trajnoise
end