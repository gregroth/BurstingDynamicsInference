

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

function remove_weakbursts(trajmodel_set, pmiss, minnbrna)
    trajnoise_set = []
    for k in eachindex(trajmodel_set)
        traj_noise = missedrna(trajmodel_set[k], pmiss,minnbrna)
        push!(trajnoise_set, traj_noise)
    end
    return trajnoise_set
end
