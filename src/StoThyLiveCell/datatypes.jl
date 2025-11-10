#####################################
#
#Define the data types used to fit the multistate model
#
#####################################



abstract type AbstractData end
abstract type AbstractDataGroup end

abstract type AbstractDataOptim end

struct DataFit{DT, D} <: AbstractDataOptim
    datatypes::DT
    data::D
    detectionLimitLC::Int
    detectionLimitNS::Int
    burstsinglet::Symbol
end


struct Survival_Burst <: AbstractData
end

struct Survival_InterBurst <: AbstractData

end

struct Prob_Burst <: AbstractData

end

struct Correlation_InterBurst <: AbstractData

end




