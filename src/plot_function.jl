import Plots
import LinearAlgebra

function func_plot(gridfile, valfile) 
    Z = 0; 
    X = 0; Y = 0;
    open("bin/data/"*valfile*".bin") do io
        n = read(io, UInt8)
        seek(io, 4)
        m = read(io, UInt8)
        seek(io, 8)

        n = convert(Int64, n[1]); m = convert(Int64, m[1])
        Z = Array{Float64}(undef, m,n)
        read!(io, Z)
    end
    n = 0; m = 0;
    open("bin/data/"*gridfile*".bin") do io
        n = read(io, UInt8)
        seek(io, 4)
        m = read(io, UInt8)
        seek(io, 8)

        n = convert(Int64, n[1]); m = convert(Int64, m[1])
        XY = Array{Float64}(undef, 2, n*m)
        read!(io, XY)
        X = reshape(XY[1,:],m,n)
        Y = reshape(XY[2,:],m,n)
    end
    # 3D style
    p = Plots.plot();
    Plots.plot!(X',Y',Z',color="blue",label=false, camera=(10,30)) 
    # Plot in both directions for grid: uncomment next line
    # Plots.plot!(X,Y,Z,color="blue",label=false, camera=(10,30))
    p
end

function save_plots()
    bins = ["function","pdx","pdy","laplace","analyticpdx","analyticpdy","analyticlaplace"]
    extras = ["","_stretch"]
    for extra in extras
        for bin in bins
            bin *= extra
            p = func_plot("grid"*extra, bin)
            Plots.savefig(p, "img/"*bin*".svg")
        end
    end
end
