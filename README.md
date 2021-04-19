# FitzHugh-Nagumo lattice model in Julia

This repo contains a simple implementation of the FitzHugh-Nagumo model of cellular excitability on a 2D lattice.

## FitzHugh-Nagumo model

$$ 
\frac{dv}{dt} = \frac{1}{c} \left( v - \frac{1}{3}v^3 + w + I_t \right) + D \nabla v \\
\frac{dw}{dt} = c \left( v - a w + b \right) \\
$$
