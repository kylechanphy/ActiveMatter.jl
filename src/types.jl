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
    CustomLogger
###
abstract type AbstractParicles end
abstract type ObstacleLattice end
abstract type Interaction end
abstract type Parameter end
abstract type AbstractSystem end
abstract type AbstractLogger end
abstract type CustomLogger <: AbstractLogger end

"""
System interface
"""
mutable struct System{T,T2,T3,T4} <: AbstractSystem
    particles::T
    interactions::T2
    parameter::T3
    loggers::T4
end

function System(;
    particles,
    interactions,
    parameter,
    loggers)

    p = particles
    in = interactions
    para = parameter
    log = loggers
    return System(p, in, para, log)
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
