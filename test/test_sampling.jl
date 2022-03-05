using BenchmarkTools
using Pkg
Pkg.activate(".")
include("../src/ActiveMatter.jl")

# using BenchmarkTools
import Base.Threads.@threads, Base.Threads.@sync, Base.Threads.@spawn
tid = Threads.threadid
nth = Threads.nthreads()

"""
sys = System(
    particles = p,
    interactions = inters,
    parameter = para,
    loggers = loggers
)
"""


d = 1
r = 0.45

pos0 = [0 * d, 0 * d]
v0 = 1
Dr = 1
ω0 = 50
p = Particle(pos0, v0)

dt = 0.001
n_step = 100_000

n_take = n_step
# loggers = Dict("pos" => PosLogger(),
# "t" => TimeLogger())

# loggers = Dict("traj" => TrajLogger())

para = ParaCAP(Dr = Dr, v0 = v0, ω0 = ω0, n_step = n_step, dt = dt)
# ob = SquareLattice(d, r)
ob = FreeSpace()
inters = ObstacleCollision(ob)

@with_kw mutable struct Logger <: CustomLogger
    coord::Vector{SV} = Vector{SV}(undef, n_take)
    vel::Vector{SV} = Vector{SV}(undef, n_take)
    t::Vector{Float64} = Vector{Int}(undef, n_take)
    coliide::Vector{Int} = Vector{Int}(undef, n_take)
end

function logging!(logger::CustomLogger, p::Particle, para::Parameter, step)
    logger.coord[step] = p.pos
    logger.vel[step] = p.vel
    logger.t[step] = step * para.dt
    logger.coliide[step] = p.collide
end

function sampling(N_sample)
    syss = []
    for i in 1:nth
        push!(syss, System(p, inters, para, Logger()))
    end

    ### sampleing   
    Deff = zeros(nth)
    msd = [zeros(n_take) for _ in 1:nth]

    @sync for i in 1:N_sample
        @spawn begin
            sys = syss[tid()]
            ## reset System
            ϕ = rand() * 2π
            sys.particles = Particle(pos0, v0, ϕ)

            runSim(sys)
            # result[tid()] += DiffusConst(sys.loggers["traj"], sys.loggers["t"])
            Deff[tid()] += DiffusConst(sys.loggers.coord, sys.loggers.t)
            msd[tid()] .+= MSD(sys.loggers.coord)
            # @show result[tid()],tid()
        end
    end
    result = [Deff, msd]
    return syss[1].loggers.t, sum.(result) ./ N_sample
end

# @time t, result = sampling(10)
plot(t, result[2], xlabel="t", ylabel="MSD",label = "")
