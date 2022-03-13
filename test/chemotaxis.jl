using Pkg
Pkg.activate("./test")

using ActiveMatter
using Parameters
using Plots

v0 = 1
lp = 1
r = 1

srctype = "const_src"
# srctype = "free"

n_step = 1000
dt = 0.01
# nx = (2v0*dt*n_step)/ (v0*dt)
ny = nx = 2*n_step

ω0 = v0 / r;
Dr = v0 / lp;
para = ParaChemoDroplet(nx = nx, ny = ny, ω0 = ω0, Dr = Dr)

pos = SV(nx / 2, ny / 2) .* (para.dx, para.dy)
p = ChemoDroplet(pos = pos, srctype = srctype, src = 1)
inter = Chemotaxis(nx, ny)

n_take = para.n_step
@with_kw mutable struct Logger <: CustomLogger
    coord::Vector{SV} = Vector{SV}(undef, n_take)
    F::Vector{SV} = Vector{SV}(undef, n_take)
    t::Vector{Float64} = Vector{Int}(undef, n_take)
end
logger = Logger()
function logging!(logger::CustomLogger, p::ChemoDroplet, para::Parameter, step)
    logger.coord[step] = p.pos
    logger.F[step] = p.force
    logger.t[step] = step*para.dt
end

# loggers = Dict("traj" => TrajLogger())

sys = sys = System(p, inter, para, logger)

println("--- simulation started ---")
@time runSim(sys)
println("--- all done ---")