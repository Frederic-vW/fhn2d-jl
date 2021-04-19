# FitzHugh-Nagumo lattice model in Julia

This repo contains a simple implementation of the FitzHugh-Nagumo model of cellular excitability on a 2D lattice.

## FitzHugh-Nagumo model

$$ 
\frac{dv}{dt} = \frac{1}{c} \left( v - \frac{1}{3}v^3 + w + I_t \right) + D \nabla v \\
\frac{dw}{dt} = c \left( v - a w + b \right) \\
$$

The command is `fhn2d(N, T, t0, a, b, c, I, sd, D, dt, stim, blocks)`:
Examples 1, 2 use
- `stim = [ [[25,50], [1,N], [3,8]], [[130,150], [n_2-2,n_2+2], [10,25]] ]`
- `blocks = [ [[2*n_4,3*n_4], [15,20]], [[2*n_4+10,3*n_4+10], [40,45]] ]`

### Example-1
Parameters:  
`N = 128, T = 1000, t0 = 0, a = 0.5, b = 0.7, c = 0.3, I = 0.5, sd = 0.02, D = 1.0, dt = 0.1`

<p align="center">
<video src="videos/fhn2d_I_0.50_sd_0.02_D_1.00.webm" width="256" height="256" controls preload></video>
</p>

### Example-2
Parameters:  
`N = 128, T = 1000, t0 = 0, a = 0.5, b = 0.7, c = 0.3, I = 1.0, sd = 0.02, D = 1.0, dt = 0.1`

<p align="center">
<video src="videos/fhn2d_I_1.00_sd_0.02_D_1.00.webm" width="256" height="256" controls preload></video>
</p>


### Example-3
Large background noise, no stimulation current `stim = []`.
Parameters on the left:  
`N = 128, T = 1000, t0 = 0, a = 0.5, b = 0.7, c = 0.3, I = 0.5, sd = 0.10, D = 1.0, dt = 0.1`  
On the right, with reduced diffusion constant `D=0.25`, other parameters identical.

<p align="center">
<video src="videos/fhn2d_I_0.50_sd_0.10_D_1.00.webm" width="256" height="256" controls preload></video>
<span>     </span>
<video src="videos/fhn2d_I_0.50_sd_0.10_D_0.25.webm" width="256" height="256" controls preload></video>
</p>

**Conclusions:**
The FHN lattice can produce:
- stimulation-induced travelling waves
- stimulation-induced spiral waves
- noise-induced spiral waves
