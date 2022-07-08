export 
    logging!,
    TimeLogger,
    TrajLogger,
    PosLogger,
    CollideLogger,
    runLogger


"""
pre-allocate the size of logger
"""
function setLogger!(loggers::Dict, para::Parameter; every::Bool = false)

    for (key, logger) in loggers
        loggers[key] = preAllocate(logger, para.n_step)
    end
end
function setLogger!(loggers::CustomLogger, para::Parameter; every::Bool = false)
    nothing
end

function setLogger!(loggers::AbstractArray, para::Parameter; every::Bool = false)
    for i in eachindex(loggers)
    loggers[i] = preAllocate(loggers[i], para.n_step)
    end
end
function setLogger!(loggers::NTuple, para::Parameter; every::Bool=false)
    nothing
end


"""
run logger
"""
function runLogger!(loggers::Dict, p, step, para::Parameter; every::Bool = false)
    for logger in values(loggers)
        logging!(logger, p::AbstractParicles, para, step)
    end
end

function runLogger!(loggers::CustomLogger, p, step, para::Parameter, inter::Interaction; every::Bool = false)
    loggers.logfunc(loggers, p::AbstractParicles, para, inter, step)
    # logging!(loggers, p::AbstractParicles, para, step)
end


function runLogger!(loggers::AbstractArray, p, step, para::Parameter; every::Bool = false)
    for i in eachindex(loggers)
        logging!(loggers[i], p::AbstractParicles, para, step)
    end
end

"""
record time step
"""
mutable struct TimeLogger <: AbstractLogger
    t::Vector{Float64}
end
TimeLogger(; every = false) = TimeLogger([])
TimeLogger(n_take::Int; every = false) = TimeLogger(Vector{Float64}(undef, n_take))

function logging!(logger::TimeLogger, p::AbstractParicles, para, step)
    logger.t[step] = para.dt * step

end
function preAllocate(logger::TimeLogger, n_take; every = false)
    logger.t = Vector{Float64}(undef, n_take)
    return logger
end


"""
record position
"""
mutable struct PosLogger <: AbstractLogger
    coord::Vector{SV}
    # t::Vector{Float64}
end
PosLogger(n_take::Int; every = false) = PosLogger(Vector{SV}(undef, n_take))
PosLogger(; every = false) = PosLogger([])

function preAllocate(logger::PosLogger, n_take; every = false)
    logger.coord = Vector{SV}(undef, n_take)

    return logger
end

function logging!(logger::PosLogger, p::AbstractParicles, para, step)
    logger.coord[step] = p.pos
    # log.t[step] = step * dt
end


"""
record traj
"""
mutable struct TrajLogger <: AbstractLogger
    coord::Vector{SV}
    vel::Vector{SV}
    t::Vector{Float64}
end
TrajLogger(n_take::Int; every = false) = TrajLogger(Vector{SV}(undef, n_take),
    Vector{SV}(undef, n_take),
    Vector{Float64}(undef, n_take))
TrajLogger(; every = false) = TrajLogger([], [], [])

function preAllocate(logger::TrajLogger, n_take; every = false)
    logger.coord = Vector{SV}(undef, n_take)
    logger.vel = Vector{SV}(undef, n_take)
    logger.t = Vector{Float64}(undef, n_take)

    return logger
end

function logging!(logger::TrajLogger, p::AbstractParicles, para, step)
    logger.coord[step] = p.pos
    logger.vel[step] = p.vel
    logger.t[step] = step * para.dt
end


"""
collision
"""

mutable struct CollideLogger <: AbstractLogger
    collide::Vector{Int}
end
CollideLogger() = CollideLogger([])

function preAllocate(logger::CollideLogger, n_take; every=false)
    logger.collide = Vector{Int}(undef, n_take)
    return logger
end

function logging!(logger::CollideLogger, p::AbstractParicles, para, step)
    logger.collide[step] = p.collide

end

"""
Save data to local
"""

function outputdata()

end