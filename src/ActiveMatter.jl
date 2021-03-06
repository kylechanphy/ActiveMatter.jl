module ActiveMatter

using LinearAlgebra
using Plots
using StaticArrays
using Parameters
using LoopVectorization

include("types.jl")
include("Parameter.jl")
include("Particles.jl")
include("Interaction/ObstacleLattice.jl")
include("Interaction/Collision.jl")
include("Interaction/Chemotaxis.jl")
include("Interaction/Chemotaxis3D.jl")
include("Simulation.jl")
include("Logger.jl")
include("Analysis.jl")
include("Plots.jl")

end

