#= using Tables, FileIO, JLD2
using Random
using LinearAlgebra
using StatsBase
 =#
#set the random seed and the maximum number of threads for the BLAS functions
# used in the paper Random.seed!(1*17)
Random.seed!(1*17)
BLAS.set_num_threads(10)
