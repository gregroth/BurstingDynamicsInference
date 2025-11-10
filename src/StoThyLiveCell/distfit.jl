#list of the distance functions used in the fit of the multi-state model
abstract type AbstractDistanceFitRNA end
abstract type AbstractDistanceFitBurst end
abstract type AbstractDistanceFitRNAandBurst end



struct LsqSurvival <: AbstractDistanceFitBurst 
    weight::Float32
end

"weighted average of the lsq on the linear data and the log data"
struct LsqSurvivalLogLin <: AbstractDistanceFitBurst 
    weight::Float32
    plog::Float32
end


struct LsqProb <: AbstractDistanceFitBurst 
    weight::Float32
end


struct LsqNumber <: AbstractDistanceFitBurst 
    weight::Float32
end




function (f::LsqSurvival)(estimate_signal::AbstractVector{T}, ref_signal::Tuple{Vector,Vector}; kwargs...) where T
    f.weight * sum((log.(estimate_signal .+ 1e-6) - log.(ref_signal[2].+ 1e-6)).^2)/length(estimate_signal)
end

function (f::LsqSurvivalLogLin)(estimate_signal::AbstractVector{T}, ref_signal::Tuple{Vector,Vector}; kwargs...) where T
    f.weight * (f.plog * sum((log.(estimate_signal .+ 1e-6) - log.(ref_signal[2].+ 1e-6)).^2) + (1-f.plog) * sum(((estimate_signal) - (ref_signal[2])).^2) )/length(estimate_signal)
end

function (f::LsqProb)(estimate_signal_tot::T, ref_signal::Float64; kwargs...)where T
    f.weight * (log.(estimate_signal_tot.+ 1e-6) - log.(ref_signal.+ 1e-6)).^2
end

function (f::LsqNumber)(estimate_signal_tot::T, ref_signal::Float64; kwargs...) where T
    f.weight * (estimate_signal_tot - ref_signal)^2
end




