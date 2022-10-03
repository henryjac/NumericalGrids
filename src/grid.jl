import Plots

function grid(gridfile, boundaryfile)
    p = Plots.plot()

    open("bin/data/"*boundaryfile*".bin") do io
        # Read first byte (int) that is how many points we have on each boundary
        precision = read(io, UInt8) 

        seek(io, 4) # Read from byte4, size of int
        input = Array{Float64}(undef, 2, (precision[1]+1)*4)
        read!(io, input)

        x = input[1,:]
        y = input[2,:]

        Plots.plot!(p,x,y,color="blue",linewidth=4,label=false)
    end

    binfile = "bin/data/"*gridfile*".bin"
    open(binfile) do io
        # # Read first byte (int) that is how many points we have on each boundary
        n = read(io, UInt8) 
        seek(io, 4) # Start from next byte for size of m
        m = read(io, UInt8)
        seek(io, 8) # pos is set to 0, so seek forward two UInts

        n = convert(Int64, n[1]); m = convert(Int64, m[1])
        input = Array{Float64}(undef, 2, n*m)
        read!(io, input)

        x = input[1,:]
        y = input[2,:]

        Plots.plot!(p,x,y,color="red",markersize=5,label=false,seriestype=:scatter)
    end
    imgfile = "img/"*gridfile*".svg"
    Plots.savefig(imgfile)
end

function save_grids()
    grid("grid","boundary")
    grid("grid_stretched","boundary")
end
