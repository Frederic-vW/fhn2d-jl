#!/usr/local/bin/julia
# last tested Julia version: 1.6.0
# Fith-Hugh Nagumo model on a 2D lattice
# FvW 03/2018

using NPZ
using PyPlot
using Statistics
using VideoIO

function array_to_mp4(filename, x)
	# BW
	y = UInt8.(round.(255*(x .- minimum(x))/(maximum(x)-minimum(x))))
	# BW inverted
	#y = UInt8.(round.(255 .- 255*(x .- minimum(x))/(maximum(x)-minimum(x))))
	encoder_options = (color_range=2, crf=0, preset="medium")
	framerate=30
	T = size(x,1)
	open_video_out(filename, y[1,end:-1:1,:], framerate=framerate,
				   encoder_options=encoder_options) do writer
		for i in range(2,stop=T,step=1)
			z = y[i,end:-1:1,:]
			write(writer, z)
		end
	end
end

function fhn2d(N, T, t0, a, b, c, I0, sd, D, dt, stim, blocks)
    c1 = 1/c
	# initialize Fitz-Hugh Nagumo system
	v, w = zeros(N,N), zeros(N,N)
	dv, dw = zeros(N,N), zeros(N,N)
	sqrt_dt = sqrt(dt)
	X = zeros(T,N,N)
	X[1,:,:] = v
	#offset = 0 # Int(round(1*nt))
	# stimulation protocol
	I = zeros(t0+T,N,N)
	for s in stim
		t_on, t_off = s[1]
		x0, x1 = s[2]
		y0, y1 = s[3]
		I[t0+t_on:t0+t_off, x0:x1, y0:y1] .= I0
	end
	# iterate
	for t in range(2, stop=t0+T, step=1)  # run...
	    (t%100 == 0) && print("    t = ", t, "\r")
		# FHN equations
		dv = D*L(v) + c1*(v - 1/3*v.*v.*v - w + I[t,:,:])
		dw = c*(v - a*w .+ b)
		# Ito stochastic integration
        v += (dv*dt + sd*sqrt_dt*randn(N,N))
		w += (dw*dt)
		# dead block(s):
		for bl in blocks
			#x0, x1 = bl[1]
			#y0, y1 = bl[2]
			v[bl[1][1]:bl[1][2], bl[2][1]:bl[2][2]] .= 0.0
			w[bl[1][1]:bl[1][2], bl[2][1]:bl[2][2]] .= 0.0
		end
		(t > t0) && (X[t-t0,:,:] = v)
        #X[t,:,:] = v
	end
	println("\n")
	return X
end

function L(x)
	# Laplace operator
    xU = circshift(x, [-1 0])
    xD = circshift(x, [1 0])
    xL = circshift(x, [0 -1])
    xR = circshift(x, [0 1])
    Lx = xU + xD + xL + xR - 4x
	nx = size(x)[1]
	ny = size(x)[2]
    Lx[1,:] .= 0.0
	Lx[nx,:] .= 0.0
	Lx[:,1] .= 0.0
	Lx[:,ny] .= 0.0
    return Lx
end

function main()
	#run(`clear`)
	println("FitzHugh-Nagumo (FHN) lattice model\n")
    N = 128
	T = 1000
	t0 = 0
    a = 0.5
    b = 0.7
    c = 0.3
    I = 1.0 # 0.5 # 1.0
    sd = 0.02
	D = 1.0
    dt = 0.1
	println("[+] FitzHugh-Nagumo 2D lattice:")
    println("[+] Lattice size N: ", N)
	println("[+] Time steps T: ", T)
	println("[+] warm-up steps t0: ", t0)
    println("[+] FHN parameter a: ", a)
    println("[+] FHN parameter b: ", b)
    println("[+] FHN parameter c: ", c)
    println("[+] Stimulation current I: ", I)
	println("[+] Diffusion coefficient D: ", D)
	println("[+] Noise std. dev. sd: ", sd)
    println("[+] Integration time step dt: ", dt)
	# auxiliary variables
	n_2 = Int(round(N/2)) # 1/2 lattice size
	n_4 = Int(round(N/4)) # 1/4 lattice size
	n_5 = Int(round(N/5)) # 1/5 lattice size
	# stim protocol, array of elements [[t0,t1], [x0,x1], [y0,y1]]
	stim = [ [[25,50], [1,N], [3,8]], [[130,150], [n_2-2,n_2+2], [10,25]] ]
	#stim = []
	# dead blocks, array of elementy [[x0,x1], [y0,y1]]
	blocks = [ [[2*n_4,3*n_4], [15,20]], [[2*n_4+10,3*n_4+10], [40,45]] ]
	# run simulation
	data = fhn2d(N, T, t0, a, b, c, I, sd, D, dt, stim, blocks)
    println("[+] data dimensions: ", size(data))
	# plot mean voltage
    m = mean(reshape(data, (T,N*N)), dims=2)
    plot(m, "-sk"); show()
	# save data
    s = rpad(I,4,'0') # stim. current amplitude as 4-char string
    data_filename = string("FitzHughNagumo2D_I_", s, ".npy")
    #npzwrite(data_filename, data)
	#println("[+] Data saved as: ", data_filename)
    mov_filename = string("FitzHughNagumo2D_I_", s, ".mp4")
	array_to_mp4(mov_filename, data[1:end,:,:])
	println("[+] Data saved as: ", mov_filename)
end


main()
