### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ e3dbdee8-7f62-4301-b848-6393954121a8
begin
	import Pkg
	Pkg.activate(".")
end

# ╔═╡ 6075d462-d1ff-4b30-9de8-2f7a1eb9ecdf
md"""#### Task 1:
1.1. From small amplitude vibrations assumption: $sin(θ) ≈ θ$. \
1.2. Then from the given two balance of linear momentum equations, one can formulate the $\boldsymbol{M}$ & $\boldsymbol{K}$ matrices, as follows:

$\begin{align}
\boldsymbol{M} = 
\begin{bmatrix}
\frac{ρAL^3}{3} & 0 \\
0 & \frac{ρAL}{3}
\end{bmatrix}   \space , \space 
\boldsymbol{K} = 
\begin{bmatrix}
\frac{ρAgL^2}{2} & 0 \\
0 & \frac{EA}{L}
\end{bmatrix}
\end{align}$

with $\boldsymbol{u}$ & $\ddot{\boldsymbol{u}}$ are vectors, and there components are as follows:

$\begin{align}
\boldsymbol{u} = 
\begin{bmatrix}
θ \\ 
u
\end{bmatrix} \space , \space

\ddot{\boldsymbol{u}} = 
\begin{bmatrix}
\ddot{θ} \\
\ddot{u}
\end{bmatrix}
\end{align}$


1.3. For given $α_1$ & $α_2$, one can calculate the Rayleigh damping matrix $\boldsymbol{M}$ from the following equation:
$\begin{align}
\boldsymbol{C} = α_1 \boldsymbol{M} + α_2 \boldsymbol{K}
\end{align}$

1.4. Finally, the semi-discrete equation can be written as follows:


$\begin{align}
\boldsymbol{M} ⋅ \ddot{\boldsymbol{u}} + \boldsymbol{C} ⋅ \dot{\boldsymbol{u}} + \boldsymbol{K} ⋅ \boldsymbol{u} = \boldsymbol{R}
\end{align}$

with, 

$\begin{align}
\boldsymbol{R} = 
\begin{bmatrix} 
0 \\
0
\end{bmatrix} \space , \space 
\dot{\boldsymbol{u}} = 
\begin{bmatrix}
\dot{θ} \\
\dot{u}
\end{bmatrix}
\end{align}$

!!! note
	Semi-discrete means that the equation is discrete in space only, but 		continous in time.
"""

# ╔═╡ 774a153f-df2e-4fd0-a8db-ba370f8057c1
begin
	# Data inputs
	Δt = 0.01 
	L = 1.0 
	ρA = 1.0
	EA = 500 
	g = 9.8
end

# ╔═╡ 0b86fb31-1e75-4c0d-8704-a95ec5318218
begin
	# Global variables
	M = [(ρA*L^3)/3 0; 0  (ρA*L)/3]
	K = [(ρA*g*L^2)/2 0; 0 EA/L]
	R = [0 , 0]
	(M = M, K=K)
end

# ╔═╡ 380f1415-27af-4b21-a93e-9fbabdd27a30
begin
	# Calculate Rayleigh damping
	α_1 = 1.0
	α_2 = 0.0
	C = α_1*M + α_2 * K
end

# ╔═╡ 67cf6c46-7e1e-4f20-bba7-d6c2b9d53cf0
md"""#### Task 2:
2.1. Calculate the initial conditions $\boldsymbol{u}_o, \dot{\boldsymbol{u}}_o, \ddot{\boldsymbol{u}}_o$ as follows:

$\begin{align}
 \dot{\boldsymbol{u}}_o = 
\begin{bmatrix}
0 \\
-\frac{L}{5}
\end{bmatrix} \space , \space 
\dot{\boldsymbol{u}}_o = 
\begin{bmatrix}
\sqrt{\frac{g}{6L}} \\
0
\end{bmatrix} \rightarrow \text{\small{both are given in the problem prompt}}
\end{align}$

whereas, $\ddot{\boldsymbol{u}}_o$ can be calculated from the semi-discrete equation of motion as follows:

$\begin{gather}
\boldsymbol{M} ⋅ \ddot{\boldsymbol{u}}_o + \boldsymbol{C} ⋅ \dot{\boldsymbol{u}}_o + \boldsymbol{K} ⋅ \boldsymbol{u}_o = \boldsymbol{R} \\
\Rightarrow 
\ddot{\boldsymbol{u}}_o = \boldsymbol{M}^{-1} ⋅ (\boldsymbol{R} - \boldsymbol{C} ⋅ \dot{\boldsymbol{u}}_o + \boldsymbol{K} ⋅ \boldsymbol{u}_o)
\end{gather}$

2.2. For given $ρ_∞$, calculate the Newmark parameters according to *Chung and Hulbert (1993)* method:

$\begin{gather}
α_m = \frac{2ρ_∞ - 1}{ρ_∞ + 1}
 \space , \space 
α_f = \frac{ρ_∞}{ρ_∞ + 1} 
\space , \space 
\beta = \frac{1}{4} [1 - α_m + α_f]^2 
\space , \space 
γ = \frac{1}{2} - α_m + α_f
\end{gather}$

2.2. Calculate $\boldsymbol{K}_{eff}$ from the following equation:

$\begin{gather}
\boldsymbol{K}_{eff} = \boldsymbol{M} \frac{1 - α_m}{βΔt^2} + \boldsymbol{C} \frac{γ(1 - α_f)}{βΔt} + \boldsymbol{K} (1 - α_f) 
→ \text{\small{(slides pg. 63)}}
\end{gather}$

2.3. Loop over time from $t_o = 0$ upto $t_{final}$ with time step = $Δt$ and for each step $n$ do: \
2.3.1 calculate the effective right hand side from the following equation:

$\begin{align}

& \boldsymbol{r}_{eff} = -  \boldsymbol{K} ⋅ α_f \boldsymbol{u}_n \\ 
& + \boldsymbol{C} ⋅ \left[
\frac{γ(1-α_f)}{βΔt} \boldsymbol{u}_n + 
\frac{γ - γα_f - β}{β} \dot{\boldsymbol{u}}_n +
\frac{(γ - 2β)(1-α_f)}{2β} Δt \ddot{\boldsymbol{u}}_n
\right] \\ 
&+ \boldsymbol{M} ⋅ \left[ 
\frac{1-α_m}{βΔ^2} \boldsymbol{u}_n
\right] → \text{\small{(slides pg. 63)}}
\end{align}$

2.3.2. solve for $\boldsymbol{u}_{n+1}$ as follows:

$\begin{align}
\boldsymbol{u}_{n+1} = \boldsymbol{K}^{-1}_{eff} ⋅ \boldsymbol{r}_{eff}
\end{align}$

2.3.3 update velocity and acceleration (i.e. $\dot{\boldsymbol{u}}_{n+1}, \ddot{\boldsymbol{u}}_{n+1}$) from the following equations: 

$\begin{align}
&\dot{\boldsymbol{u}}_{n+1} = \frac{γ}{βΔt} (\boldsymbol{u}_{n+1} - \boldsymbol{u}_n) - \frac{γ- β}{β} \dot{\boldsymbol{u}}_n - \frac{γ - 2β}{2β} Δt \ddot{\boldsymbol{u}}_n
 → \text{\small{(slides pg. 60)}} \\

&\ddot{\boldsymbol{u}}_{n+1} = \frac{1}{βΔt^2} (\boldsymbol{u}_{n+1} - \boldsymbol{u}_n) - \frac{1}{βΔt} \dot{\boldsymbol{u}}_n - 
\frac{1-2β}{2β} \ddot{\boldsymbol{u}}_n → \text{\small{(slides pg. 60)}}
\end{align}$

"""

# ╔═╡ 3cea61f0-bc04-4092-8e3b-7803932865d5
begin
	# Initial conditions
	u_o = [0 , - L /5 ]
	ud_o = [sqrt(g/(6*L)),0]
end

# ╔═╡ 9d7e83e5-4d75-4d96-913f-a656d9866735
md"""
!!! note
	Test.
"""

# ╔═╡ 20193589-4189-4347-ac5d-0e4370152ce6


# ╔═╡ 4ac90089-eee4-4ce3-a2e8-0152877f4414
md"""

$\begin{align}
-\bigl(e^{-t} (\cos t + \sin t)\bigr)^2 \\
+ \bigl(e^{-t} (\cos t + \sin t)\bigr)^2 \\
+ \bigl(-e^{-t}\bigr)^2

\end{align}$
"""

# ╔═╡ Cell order:
# ╠═e3dbdee8-7f62-4301-b848-6393954121a8
# ╟─6075d462-d1ff-4b30-9de8-2f7a1eb9ecdf
# ╠═774a153f-df2e-4fd0-a8db-ba370f8057c1
# ╠═0b86fb31-1e75-4c0d-8704-a95ec5318218
# ╠═380f1415-27af-4b21-a93e-9fbabdd27a30
# ╠═67cf6c46-7e1e-4f20-bba7-d6c2b9d53cf0
# ╠═3cea61f0-bc04-4092-8e3b-7803932865d5
# ╠═9d7e83e5-4d75-4d96-913f-a656d9866735
# ╠═20193589-4189-4347-ac5d-0e4370152ce6
# ╠═4ac90089-eee4-4ce3-a2e8-0152877f4414
