module ActiveMatter

using LinearAlgebra
using Plots
using StaticArrays
using Parameters
using LoopVectorization
const SV = SVector{2,Float64}

include("types.jl")
include("Parameter.jl")
include("Particles.jl")
include("Interaction/ObstacleLattice.jl")
include("Interaction/Collision.jl")
include("Interaction/Chemotaxis.jl")
include("Simulation.jl")
include("Logger.jl")
include("Analysis.jl")
include("Plots.jl")

end

