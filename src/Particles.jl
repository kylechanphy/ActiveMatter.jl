"""
Define the type of particles 
"""

export
    SV,
    AbstractParicles,
    Particle


###




"""
Define a noramal particle
"""
mutable struct Particle <: AbstractParicles
    pos::SV
    vel::SV
    cell_pos::SV
    collide::Int
end
Particle(x0 = 0.0, y0 = 0.0, v0 = 0, ϕ0 = 0.0) = Particle(SV(x0, y0), v0 * SV(cos(ϕ0), sin(ϕ0)), SV(x0, y0), 0)
Particle(pos::AbstractArray, v0 = 0, ϕ0 = 0.0) = Particle(SV(pos[1], pos[2]), v0 * SV(cos(ϕ0), sin(ϕ0)), SV(pos[1], pos[2]), 0)
