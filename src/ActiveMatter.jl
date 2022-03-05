# module ActiveMatter

using LinearAlgebra
using Plots
using StaticArrays
using Parameters
const SV = SVector{2,Float64}

include("types.jl")
include("Particles.jl")
include("Interaction/ObstacleLattice.jl")
include("Interaction/Collision.jl")
include("Simulation.jl")
include("Parameter.jl")
include("Logger.jl")
include("Analysis.jl")
include("Plots.jl")

# end

