using Pkg
Pkg.activate(".")
# include("../src/ActiveMatter.jl")
using ActiveMatter

Pkg.activate("./test")
using Revise
using Parameters
using Plots
using BenchmarkTools 

v0 = 1
lp = 5
r = 1

srctype = "const_src"
# srctype = "free"

dt = 0.05*lp/v0
n_step = 120
# dt = 0.01
# nx = (2v0*dt*n_step)/ (v0*dt)
ny = nx = nz =  Int(round(2*n_step))

ω0 = v0 / r;
Dr = v0 / lp;
para = ParaChemoDroplet(dt = dt, nx = nx, ny = ny, nz = nz, ω0 = ω0, Dr = Dr, n_step = n_step)
# para.α = 0

pos = SV3(nx / 2, ny / 2, nz /2) .* (para.dx, para.dy, para.dz)
p = ChemoDroplet3D(pos = pos, srctype = srctype, src = 1)
inter = Chemotaxis3D(nx,ny,nz)

function logging!(logger::CustomLogger, p::ChemoDroplet3D, para::Parameter, step)
    logger.coord[step] = p.pos
    logger.F[step] = p.force
    logger.t[step] = step * para.dt
end
n_take = para.n_step
@with_kw mutable struct Logger3D <: CustomLogger
    coord::Vector{SV3} = Vector{SV3}(undef, n_take)
    F::Vector{SV3} = Vector{SV3}(undef, n_take)
    t::Vector{Float64} = Vector{Int}(undef, n_take)

    logfunc::Function = logging!
end
logger = Logger3D()


# logger = Dict("traj" => TrajLogger())

sys = sys = System(p, inter, para, logger)

println("--- simulation started ---")
@show sys.parameter
@time runSim(sys)
println("--- all done ---")