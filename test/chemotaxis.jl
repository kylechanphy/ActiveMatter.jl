using Pkg
Pkg.activate(".")
# include("../src/ActiveMatter.jl")
using ActiveMatter

Pkg.activate("./test")
using Revise
using Parameters
using Plots

v0 = 1
lp = 1
r = 1

srctype = "const_src"
# srctype = "free"

dt = 0.05*lp/v0
n_step = 500
# dt = 0.01
# nx = (2v0*dt*n_step)/ (v0*dt)
ny = nx =  n_step

ω0 = v0 / r;
Dr = v0 / lp;
para = ParaChemoDroplet(dt = dt, nx = nx, ny = ny, ω0 = ω0, Dr = Dr, n_step = n_step)
# para.α = 0

pos = SV(nx / 2, ny / 2) .* (para.dx, para.dy)
p = ChemoDroplet(pos = pos, srctype = srctype, src = 1)
inter = Chemotaxis(nx, ny)

function logging!(logger::CustomLogger, p::ChemoDroplet, para::Parameter, step)
    logger.coord[step] = p.pos
    logger.F[step] = p.force
    logger.t[step] = step * para.dt
end
n_take = para.n_step
@with_kw mutable struct Logger <: CustomLogger
    coord::Vector{SV} = Vector{SV}(undef, n_take)
    F::Vector{SV} = Vector{SV}(undef, n_take)
    t::Vector{Float64} = Vector{Int}(undef, n_take)

    log::Function = logging!
end
logger = Logger()


loggers = Dict("traj" => TrajLogger())

sys = sys = System(p, inter, para, logger)

println("--- simulation started ---")
@show sys.parameter
@time runSim(sys)
println("--- all done ---")