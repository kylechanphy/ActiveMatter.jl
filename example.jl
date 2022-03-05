using Pkg
Pkg.activate(".")
include("src/ActiveMatter.jl")
include("setPara.jl")

using Plots;
gr(show = true);



### input parameter

lp = 100
r = 1
pack = 0.7
ob = SquareLattice()

p, inters, para = setPara(pack, lp, r, ob)


### custom logger -- higer performence than Dict
@with_kw mutable struct Logger <: CustomLogger
    coord::Vector{SV} = Vector{SV}(undef, para.n_step)
    vel::Vector{SV} = Vector{SV}(undef, para.n_step)
    t::Vector{Float64} = Vector{Int}(undef, para.n_step)
    coliide::Vector{Int} = Vector{Int}(undef, para.n_step)
end
loggers = Logger()

### define what data want to log 
function logging!(logger::CustomLogger, p::Particle, para::Parameter, step)
    logger.coord[step] = p.pos
    logger.vel[step] = p.vel
    logger.t[step] = step * para.dt
    logger.coliide[step] = p.collide
end

### crete sys object 
sys = System(p, inters, para, loggers)

### run simulation 
runSim(sys)

### plot trajactory 
pic = plot(ob, 4, loggers.coord)
display(pic)