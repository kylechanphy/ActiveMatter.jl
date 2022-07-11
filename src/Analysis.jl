export
    MSD,
    slope,
    DiffusConst,
    FFT

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


function FFT(signal, para::ParaChemoDroplet)
    N = length(signal)
    dt = para.dt
    t = dt: dt : dt*N
    freq, F = FFT(signal, t)

    return freq, F
end

function FFT(signal, t)
    dt = t[2] - t[1]
    N = length(t)
    freq = rfftfreq(N, 1/dt)
    F = abs.(rfft(signal))
    
    return freq, F
end