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

p, inters, para, ob = setPara(pack, lp, r, ob)

loggers = Dict("traj" => TrajLogger(), "coll"=>CollideLogger())
### crete sys object 
sys = System(p, inters, para, loggers)

### run simulation 
runSim(sys)

### plot trajactory 
pic = plot(ob, 3, loggers["traj"].coord)
display(pic)