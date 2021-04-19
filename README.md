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

**Example-1**
`N = 128, T = 1000, t0 = 0, a = 0.5, b = 0.7, c = 0.3, I = 0.5, sd = 0.02, D = 1.0, dt = 0.1`
<video src=videos/FitzHughNagumo2D_I_0.50.mp4>

**Example-2**
`N = 128, T = 1000, t0 = 0, a = 0.5, b = 0.7, c = 0.3, I = 1.0, sd = 0.02, D = 1.0, dt = 0.1`
<video src=videos/FitzHughNagumo2D_I_1.00.mp4>
