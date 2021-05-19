#!/usr/local/bin/julia
# last tested Julia version: 1.6.1
# Fith-Hugh Nagumo model on a 2D lattice
# FvW 03/2018

using NPZ
using PyCall
using PyPlot
using Statistics
using VideoIO
@pyimport matplotlib.animation as anim

function fhn2d(N, T, t0, dt, s, D, a, b, c, I0, stim, blocks)
    c1 = 1/c
    # initialize Fitz-Hugh Nagumo system
    v = zeros(Float64,N,N)
    w = zeros(Float64,N,N)
    dv = zeros(Float64,N,N)
    dw = zeros(Float64,N,N)
    s_sqrt_dt = s*sqrt(dt)
    X = zeros(Float64,T,N,N)
    # stimulation protocol
    I = zeros(Float64,t0+T,N,N)
    for st in stim
        t_on, t_off = st[1]
        x0, x1 = st[2]
        y0, y1 = st[3]
        I[t0+t_on:t0+t_off, x0:x1, y0:y1] .= I0
    end
    # iterate
    for t in range(2, stop=t0+T, step=1)
        (t%100 == 0) && print("    t = ", t, "\r")
        # FHN equations
        dv = c1*(v - 1/3*v.*v.*v - w + I[t,:,:]) + D*L(v)
        dw = c*(v - a*w .+ b)
        # Ito stochastic integration
        v += (dv*dt + s_sqrt_dt*randn(N,N))
        w += (dw*dt)
        # dead block(s):
        for bl in blocks
            #x0, x1 = bl[1]
            #y0, y1 = bl[2]
            v[bl[1][1]:bl[1][2], bl[2][1]:bl[2][2]] .= 0.0
            w[bl[1][1]:bl[1][2], bl[2][1]:bl[2][2]] .= 0.0
        end
        (t > t0) && (X[t-t0,:,:] = v)
    end
    println("\n")
    return X
end

function animate_pyplot(fname, data)
    """
    Animate 3D array as .mp4 using PyPlot, save as `fname`
    array dimensions:
        1: time
        2, 3: space
    """
    nt, nx, ny = size(data)
    vmin = minimum(data)
    vmax = maximum(data)
    # setup animation image
    fig = figure(figsize=(6,6))
    axis("off")
    t = imshow(data[1,:,:], origin="lower", cmap=ColorMap("gray"),
			   vmin=vmin, vmax=vmax)
    tight_layout()
    # frame generator
    println("[+] animate")
    function animate(i)
        (i%100 == 0) && print("    t = ", i, "/", nt, "\r")
        t.set_data(data[i+1,:,:])
    end
    # create animation
    ani = anim.FuncAnimation(fig, animate, frames=nt, interval=10)
    println("\n")
    # save animation
    ani[:save](fname, bitrate=-1,
               extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
    show()
end

function animate_video(fname, data)
    """
    Animate 3D array as .mp4 using VideoIO, save as `fname`
    array dimensions:
        1: time
        2, 3: space
    """
    # BW
    y = UInt8.(round.(255*(data .- minimum(data))/(maximum(data)-minimum(data))))
    # BW inverted
    #y = UInt8.(round.(255 .- 255*(data .- minimum(data))/(maximum(data)-minimum(data))))
    encoder_options = (color_range=2, crf=0, preset="medium")
    framerate=30
    T = size(data,1)
    open_video_out(fname, y[1,end:-1:1,:], framerate=framerate,
                   encoder_options=encoder_options) do writer
        for i in range(2,stop=T,step=1)
            write(writer, y[i,end:-1:1,:])
        end
    end
end

function L(x)
    # Laplace operator
    # periodic boundary conditions
    xU = circshift(x, [-1 0])
    xD = circshift(x, [1 0])
    xL = circshift(x, [0 -1])
    xR = circshift(x, [0 1])
    Lx = xU + xD + xL + xR - 4x
    # non-periodic boundary conditions
    Lx[1,:] .= 0.0
    Lx[end,:] .= 0.0
    Lx[:,1] .= 0.0
    Lx[:,end] .= 0.0
    return Lx
end

function main()
    println("FitzHugh-Nagumo (FHN) lattice model\n")
    N = 128
    T = 1000
    t0 = 0
    dt = 0.1
    s = 0.02 # 0.02 # 0.10
    D = 1.0
    # FitzHugh-Nagumo parameters
    a = 0.5
    b = 0.7
    c = 0.3
    I = 1.0 # 1.0, 0.5
    println("[+] Lattice size N: ", N)
    println("[+] Time steps T: ", T)
    println("[+] Warm-up steps t0: ", t0)
    println("[+] Integration time step dt: ", dt)
    println("[+] Noise std. dev.: ", s)
    println("[+] Diffusion coefficient D: ", D)
    println("[+] FHN parameter a: ", a)
    println("[+] FHN parameter b: ", b)
    println("[+] FHN parameter c: ", c)
    println("[+] Stimulation current I: ", I)

    # stim protocol, array of elements [[t0,t1], [x0,x1], [y0,y1]]
    n_2 = Int(round(N/2)) # 1/2 lattice size
    n_4 = Int(round(N/4)) # 1/4 lattice size
    n_5 = Int(round(N/5)) # 1/5 lattice size
    stim = [ [[25,50], [1,N], [3,8]],
             [[130,150], [n_2-2,n_2+2], [10,25]] ]
    #stim = []

    # dead blocks, array of elementy [[x0,x1], [y0,y1]]
    blocks = [ [[2*n_4,3*n_4], [15,20]], [[2*n_4+10,3*n_4+10], [40,45]] ]
    #blocks = []

    # run simulation
    data = fhn2d(N, T, t0, dt, s, D, a, b, c, I, stim, blocks)
    println("[+] Data dimensions: ", size(data))

    # plot mean voltage
    m = mean(reshape(data, (T,N*N)), dims=2)
    plot(m, "-sk"); show()

    # save data
    I_str = rpad(I, 4, '0') # stim. current amplitude as 4-char string
    s_str = rpad(s, 4, '0') # noise as 4-char string
    D_str = rpad(D, 4, '0') # diffusion coefficient as 4-char string
    fname1 = string("fhn2d_I_", I_str, "_s_", s_str, "_D_", D_str, ".npy")
    #npzwrite(data_filename, data)
    #println("[+] Data saved as: ", data_filename)

    # video
    fname2 = string("fhn2d_I_", I_str, "_s_", s_str, "_D_", D_str, ".mp4")
    #animate_pyplot(fname2, data) # slow
    animate_video(fname2, data) # fast
    println("[+] Data saved as: ", fname2)
end

main()
