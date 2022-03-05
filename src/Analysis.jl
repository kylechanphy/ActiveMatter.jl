"""
"""
function MSD(coord::AbstractArray)
    coord .-= (coord[1],)
    return norm.(coord) .^ 2
end

function slope(data, time::AbstractArray)
    return time \ data
end

"""
Diffusion constant
"""

# function DiffusConst(PosLogger, TimeLogger, ;dims=2)
#     msd = MSD(PosLogger)
#     return slope(msd, TimeLogger) ./ (2*dims)
# end

# function DiffusConst(TrajLogger; dims=2)
#     msd = MSD(TrajLogger)
#     t = TrajLogger.t
#     slope = t \ msd
#     return slope ./ 2*(dims)
# end

function DiffusConst(coord, time)
    dims = length(coord[1])
    msd = MSD(coord)
    return slope(msd, time) ./ (2 * dims)
end
