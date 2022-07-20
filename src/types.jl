"""
Define type 
"""

export
    AbsAbstractParicles,
    AbstractParticle,
    AbstractLogger,
    ObstacleLattice,
    Interaction,
    Parameter,
    System,
    AbstractSystem,
    SV,
    SV3,
    CustomLogger
###
abstract type AbstractParicles end
abstract type ObstacleLattice end
abstract type Interaction end
abstract type Parameter end
abstract type AbstractSystem end
abstract type AbstractLogger end
abstract type CustomLogger <: AbstractLogger end

const SV = SVector{2,Float64}
const SV3 = SVector{3,Float64}

"""
System interface
"""
mutable struct System{T,T2,T3,T4,T5} <: AbstractSystem
    particles::T
    interactions::T2
    parameter::T3
    loggers::T4
    saving::T5
end


# @with_kw mutable struct System <: AbstractSystem
#     particles
#     interactions
#     parameter
#     saving

# end

function System(;
    particles,
    interactions,
    parameter,
    loggers,
    saveing)

    p = particles
    in = interactions
    para = parameter
    log = loggers
    save = saveing
    return System(p, in, para, log, save)
end

# function setSystem(;
#                     particles::AbstractParicles,
#                     interactions,
#                     parameter,
#                     logger)
#     p = particles
#     in = interactions
#     para = parameter
#     log = logger

#     return System(p,in,para,log )
# end
