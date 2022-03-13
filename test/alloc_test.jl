using Pkg
Pkg.activate("./test")

using BenchmarkTools
using ActiveMatter


d = 1
R = 0.45

pos0 = [0.5, 0.5] .+ randn(2) * d / 2
v0 = 1
ϕ0 = 2π * rand()

Dr = 1
ω0 = 50
p = Particle(pos0, v0, ϕ0)

dt = 0.001
n_step = 10_000

n_take = n_step


@with_kw mutable struct Logger <: CustomLogger
    coord::Vector{SV} = Vector{SV}(undef, n_take)
    vel::Vector{SV} = Vector{SV}(undef, n_take)
    t::Vector{Float64} = Vector{Int}(undef, n_take)
    coliide::Vector{Int} = Vector{Int}(undef, n_take)
end


loggers = Logger()
function logging!(logger::CustomLogger, p::Particle, para::Parameter, step)
    logger.coord[step] = p.pos
    logger.vel[step] = p.vel
    logger.t[step] = step * para.dt
    logger.coliide[step] = p.collide
end

# loggers = [CollideLogger(), TrajLogger()]
# loggers = Dict("collide" => CollideLogger(), "trak" => TrajLogger())
# loggers = Dict("traj" => TrajLogger())

para = ParaCAP(Dr = Dr, v0 = v0, ω0 = ω0, n_step = n_step, dt = dt)
ob = SquareLattice(d, R)
# ob = FreeSpace()
inters = ObstacleCollision(ob)

sys = System(p, inters, para, loggers)

@show ob

@btime runSim($sys)


