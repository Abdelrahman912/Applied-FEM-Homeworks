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

1.4. Finally, the semi-discrete equation of motion will be as follows
"""

# ╔═╡ 9d7e83e5-4d75-4d96-913f-a656d9866735
md"""
!!! note
	Test.
"""

# ╔═╡ 20193589-4189-4347-ac5d-0e4370152ce6


# ╔═╡ Cell order:
# ╠═e3dbdee8-7f62-4301-b848-6393954121a8
# ╠═6075d462-d1ff-4b30-9de8-2f7a1eb9ecdf
# ╠═9d7e83e5-4d75-4d96-913f-a656d9866735
# ╠═20193589-4189-4347-ac5d-0e4370152ce6
