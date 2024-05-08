### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ e3dbdee8-7f62-4301-b848-6393954121a8
begin
	import Pkg
	Pkg.activate(".")
	Pkg.add("Plots")
end

# ╔═╡ c5bf3248-2ff1-4257-ba34-8b322b667abb
begin
	using LinearAlgebra
	using Plots
	using LaTeXStrings
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

# ╔═╡ 45b7da6f-7b81-452c-8274-fa797698d35f
struct DynamicResponse
	u::Matrix{Float64} # displacement
	ud::Matrix{Float64} # velocity
	udd::Matrix{Float64} # acceleration
end

# ╔═╡ b1057602-76fe-469c-b382-05c65fb14e60
struct Error
	e_abs::Vector{Float64} # absolute error (Zienkiewicz and Xie)
	η::Vector{Float64} # relative error
	e_cum::Vector{Float64} # cummulative error
end

# ╔═╡ 7b31aad9-d4f3-40d7-a7c7-18d25cb19139
struct TimeData
	times::Vector{Float64}
	steps::Vector{Float64}
end

# ╔═╡ 3cea61f0-bc04-4092-8e3b-7803932865d5
function generalized_alpha(tf,Δt, α₁, α₂, ρ_∞)
	# calculate damping
	C = α₁*M + α₂ * K
	
	# Initial conditions 
	u_o = [0 , - L /5 ] 
	ud_o = [sqrt(g/(6*L)),0]
	udd_o = inv(M) * (R - C * ud_o - K * u_o)

	# Chung and Hulbert (1993)
	α_m = (2*ρ_∞ - 1)/(ρ_∞ + 1)
	α_f = (ρ_∞)/(ρ_∞+1)
	β = 0.25 * (1 - α_m + α_f)^2
	γ = 0.5 - α_m + α_f

	K_eff = M * ((1-α_m)/(β*Δt^2)) + C * (γ * (1 - α_f))/(β*Δt) + K * (1 - α_f)

	# create time interval range to loop over later
	t = 0:Δt:tf
	n = length(t)
	steps = fill(Δt,n)

	# define our respone containers
	n = length(t)
	u = zeros(2,n)
	ud = zeros(2,n)
	udd = zeros(2,n)

	# set initial conditions in our response containers
	u[:,1] = u_o
	ud[:,1] = ud_o
	udd[:,1] = udd_o

	# define the containers for errors
	e_abs = zeros(n) # absolute error
	η = zeros(n) # relative error
	e_cum = zeros(n) # cummulative error
	for i = 1:n-1

		r_eff = -K * α_f * u[:,i] + 
				C * (((γ * (1-α_f))/(β*Δt)) * u[:,i] + ((γ - γ*α_f - β)/(β)) * ud[:,i] + (((γ-2*β)*(1-α_f))/(2*β)) * Δt*udd[:,i]) +
				M * (((1-α_m)/(β*Δt^2))*u[:,i] + ((1-α_m)/(β*Δt))* ud[:,i] + ((1-α_m - 2 *β)/(2*β)) * udd[:,i])

		# solve for u.
		u[:,i+1] = K_eff \ r_eff

		# update velocity and acceleration
		ud[:,i+1] = ((γ)/(β*Δt)) * (u[:,i+1] - u[:,i]) - ((γ - β)/(β))*ud[:,i] - ((γ - 2*β)/(2*β)) * Δt * udd[:,i]

		udd[:,i+1] = (1/(β*Δt^2))*(u[:,i+1] - u[:,i]) - (1/(β*Δt)) * ud[:,i] - ((1-2*β)/(2*β)) * udd[:,i]

		# calculate errors 
			e_abs[i+1] = norm(((6*β - 1)/(6))*(udd[:,i+1] - udd[:,i]) * Δt^2)
		η[i+1] = e_abs[i+1] / norm(u[i+1] - u[i])
		e_cum[i+1] = sum(e_abs)
	end

	(response = DynamicResponse(u,ud,udd), errors = Error(e_abs,η,e_cum),time = TimeData(t,steps))
end

# ╔═╡ 9d7e83e5-4d75-4d96-913f-a656d9866735
generalized_alpha(5.0,0.01, 1.0, 0.0, 1.0)

# ╔═╡ e0ced7bc-baae-4c15-8b6e-9e50be04ced9
md"""#### Task 3a: 
Given: $t_{final} = 5 \space s, α_1 = 1, α_2 = 0, ρ_∞ = 1.0$ 
"""

# ╔═╡ f27a727f-f80c-4dd2-bdef-63c7f3ca587f
resp_a , err_a ,time_a= generalized_alpha(5.0,0.01, 1.0, 0.0, 1.0)

# ╔═╡ 2a801b3c-49a3-4aa7-be75-452711c4d4a0
begin
	plot(time_a.times,resp_a.u[1,:],title="Displacement response with time (Task 3a)",xlabel=L"time",ylabel=L"displacement",label=L"$\theta$")
	plot!(time_a.times,resp_a.u[2,:],label=L"u")

end

# ╔═╡ 20193589-4189-4347-ac5d-0e4370152ce6
md"""#### Task 3b: 
Given: $t_{final} = 5 \space s, α_1 = 0, α_2 = 0, ρ_∞ = 0.1$ 
"""

# ╔═╡ 7cb8658f-9a03-4807-8e18-3775b3f675aa
resp_b , err_b, time_b= generalized_alpha(5.0,0.01, 0.0, 0.0, 0.1)

# ╔═╡ 4ac90089-eee4-4ce3-a2e8-0152877f4414
begin
	plot(time_b.times,resp_b.u[1,:],title="Displacement response with time (Task 3b)",xlabel=L"time",ylabel=L"displacement",label=L"$\theta$")
	plot!(time_b.times,resp_b.u[2,:],label=L"u")
end

# ╔═╡ 5309e3a3-ab35-4cc4-88cc-119813ab62a5
md"""#### Task 3 (Explanation):
##### Observtion:
In *Task 3a* **both** degree of freedoms (i.e. $u$ & $θ$) are being damped with time, whereas in *Task 3b* **only $u$** is being damped with time and the amplitude of $\theta$ is constant through time.

##### Explanation:
The reason behind such observation is that $α_1$ & $α_2$ are called **physical damping parameters** because they contribute to the formulation of the *Rayleigh damping* (i.e. $\boldsymbol{C} = α_1 \boldsymbol{M} + α_2 \boldsymbol{K}$) and as given in *Task 3a* $α_1 \neq 0 → C \neq 0$, consequently, all degree of freedoms are being damped by $\boldsymbol{C}$ regardless their frequencies. \

However, $ρ_∞$ is called **numerical damping parameter** due to the numerical error arised from this parameter that leads to damping. Additionally, it only damps degree of freedoms that has higher frequencies and this obvious in *task 3b* which has **no physical** damping but has **numerical** damping. \
Finaly, in *task 3a* there is no numerical damping because, $γ = \frac{1}{2} → ρ_∞ = 1$ (slides pg. 59) and numerical damping only occurs if $γ > \frac{1}{2} → ρ_∞ < 1$
"""

# ╔═╡ e3ad5cc3-4843-4d80-a0ba-99278ca335be
md"""#### Task 4 (Error Calculation):

!!! note
	Error calculation is already implemented in `generalized_alpha`. Accordingly, in this section, I will only show the equations that were used, some notes about regarding the implementation and the error graphs as well.

##### 4.1. Zienkiewicz and Xie (i.e. absolute error):
$\begin{align}
\boldsymbol{e}^{ZX} = \frac{6β-1}{6} (\ddot{\boldsymbol{u}}_{n+1} - \ddot{\boldsymbol{u}}_n)Δt^2 → \text{\small{(slides pg. 90)}}
\end{align}$

##### 4.2. Relative error:

$\begin{align}
η = \frac{||\boldsymbol{e}||}{||\boldsymbol{u}_{n+1} - \boldsymbol{u}_n||} → \text{\small{(slides pg. 91)}}
\end{align}$

##### 4.3. Cummulative error:

$\begin{align}

e_{cm} = \sum_{i=1}^{n}{||\boldsymbol{e}||}

\end{align}$

where $n$ is the index of the current time step.


"""

# ╔═╡ 590a6c46-c9ee-4976-afdb-e873bd0af809
begin
	plot(time_a.times,err_a.e_abs,title="Absolute Error (Task 3a)",xlabel=L"t",ylabel=L"||e||",label=L"$||e||$")
end

# ╔═╡ f81d3f7a-3543-4cda-baad-198410bb9aef
begin
	plot(time_a.times,err_a.η,title="Relative Error (Task 3a)",xlabel=L"t",ylabel=L"η",label=L"$\eta$")
end

# ╔═╡ 31592de4-c0f4-40b8-a940-06a1c83053c2
begin
	plot(time_a.times,err_a.e_cum,title="Cumulative Error (Task 3a)",xlabel=L"t",ylabel=L"e_{cm}",label=L"$e_{cm}$")
end

# ╔═╡ 8b058921-b41b-470b-b45a-35af6149428f
begin
	plot(time_b.times,err_b.e_abs,title="Absolute Error (Task 3b)",xlabel=L"t",ylabel=L"||e||",label=L"$||e||$")
end

# ╔═╡ 85eda22b-ac0f-49ae-aaf2-2aca6a6e7070
begin
	plot(time_b.times,err_b.η,title="Relative Error (Task 3b)",xlabel=L"t",ylabel=L"η",label=L"$\eta$")
end

# ╔═╡ eb4a70b8-e07f-4311-8d5c-b1bb879ba726
begin
	plot(time_b.times,err_b.e_cum,title="Cumulative Error (Task 3b)",xlabel=L"t",ylabel=L"e_{cm}",label=L"$e_{cm}$")
end

# ╔═╡ 1286dcc1-d323-4a5f-8dbc-f874c2ad5377
md"""#### Task 5 (Adaptive time stepping algorithim):

5.1. We start by setting the current time to the initial time $t ← t_o$, and the current time step with the initil time step $Δt ← Δt_o$.
5.2. Loop as long as current time is less or equal final time $(t ≤ t_f)$:
5.2.1. calculate the same objects as in the `generalized_alpha`.

!!! note
	Unlike `generalized_alpha`, $K_{eff}$ now has to be calculated inside the loop, because $Δt$ is now changing.

5.2.2. At the end of the loop we check the relative error $(η)$, whether it is inside the boundary $[ν_1η_e,ν_2η_e]$ or outside, and if it is outside, update $Δt$:

$\begin{gather}
Δt ←  
\begin{cases}
Δt \sqrt{\frac{η_e}{η}}, & \text{for } η ≤ ν_1η_e \text{ or } η ≥ ν_2η_e \\
Δt, & \text{otherwise}
\end{cases}
\end{gather}$

5.2.3. Update current time $(t)$:

$\begin{align}
t ← t + Δt
\end{align}$

"""

# ╔═╡ 205bd7ff-ab41-4ce2-bec8-6bfa411f21d0
struct AdaptiveTimeBoundary
	ν₁::Float64
	ν₂::Float64
	ηₑ::Float64
end

# ╔═╡ e070700c-5bb9-4c62-9fe8-ac7b4eed7284
function generalized_alpha_adaptive(tf,Δtₒ, α₁, α₂, ρ_∞,boundary::AdaptiveTimeBoundary)
	# calculate damping
	C = α₁*M + α₂ * K
	
	# Initial conditions 
	u_o = [0 , - L /5 ] 
	ud_o = [sqrt(g/(6*L)),0]
	udd_o = inv(M) * (R - C * ud_o - K * u_o)

	# Chung and Hulbert (1993)
	α_m = (2*ρ_∞ - 1)/(ρ_∞ + 1)
	α_f = (ρ_∞)/(ρ_∞+1)
	β = 0.25 * (1 - α_m + α_f)^2
	γ = 0.5 - α_m + α_f

	

	# create time interval range to loop over later
	ν₁ = boundary.ν₁
	ν₂ = boundary.ν₂
	ηₑ = boundary.ηₑ
	lb = ν₁ * ηₑ 
	ub = ν₂ * ηₑ
	tₒ = 0.0 # initial time
	t_current = tₒ # current time

	# define our respone containers
	u = zeros(2,0)
	ud = zeros(2,0)
	udd = zeros(2,0)
	t = [tₒ]
	tstep = [Δtₒ]
	Δt = Δtₒ # current time step 

	# set initial conditions in our response containers
	u = hcat(u,u_o)
	ud =hcat( ud,ud_o)
	udd = hcat(udd,udd_o)

	# define the containers for errors
	e_abs = [0.0] # absolute error
	η = [0.0] # relative error
	e_cum = [0.0] # cummulative error
	while t_current <= tf
		
		K_eff = M * ((1-α_m)/(β*Δt^2)) + C * (γ * (1 - α_f))/(β*Δt) + K * (1 - α_f)

		uₙ = u[:,end] # current displacement
		udₙ = ud[:,end] # current velocity
		uddₙ = udd[:,end] # current acceleration
		
		r_eff = -K * α_f * uₙ + 
				C * (((γ * (1-α_f))/(β*Δt)) * uₙ + ((γ - γ*α_f - β)/(β)) * udₙ + (((γ-2*β)*(1-α_f))/(2*β)) * Δt*uddₙ) +
				M * (((1-α_m)/(β*Δt^2))*uₙ + ((1-α_m)/(β*Δt))* udₙ + ((1-α_m - 2 *β)/(2*β)) * uddₙ)

		# solve for u.
		uₙ₁ = K_eff \ r_eff # u_{n+1} next displcement
		u = hcat(u, uₙ₁) 

		# update velocity and acceleration
		udₙ₁ = ((γ)/(β*Δt)) * (uₙ₁ - uₙ) - ((γ - β)/(β))*udₙ - ((γ - 2*β)/(2*β)) * Δt * uddₙ
		ud = hcat(ud,udₙ₁) 

		uddₙ₁ = (1/(β*Δt^2))*(uₙ₁ - uₙ) - (1/(β*Δt)) * udₙ - ((1-2*β)/(2*β)) * uddₙ
		udd = hcat(udd,uddₙ₁)

		# calculate errors 
			e_absₙ₁ = norm(((6*β - 1)/(6))*(uddₙ₁ - uddₙ) * Δt^2)
			push!(e_abs,e_absₙ₁)
		ηₙ₁ = e_absₙ₁ / norm(uₙ₁ - uₙ)
		push!(η,ηₙ₁)
		push!(e_cum, sum(e_abs))

		# adapting the time step
		if ηₙ₁  > ub || ηₙ₁ < lb 
			# here means Δt is either to big () 
			Δt = Δt * sqrt(ηₑ/ηₙ₁)
		end
		t_current += Δt
		push!(t,t_current)
		push!(tstep,Δt)
		
	end

	(response =  DynamicResponse(u,ud,udd), error =  Error(e_abs,η,e_cum),time = TimeData(t,tstep))
end

# ╔═╡ 0fdf106a-c40d-420d-a848-bf727fb49794
md"""#### Task 6:
Given $α_1 = 1$, $α_2 = 0$, and $ρ_∞ = 1$, It's required to:
* Solve the system for $5.0$ s.
* Plot the dynamic response, error, and the evolution of time step.
* Compare the results with those from task 3a.
* Cumulative error for both cases after $5.0$ s and the number of time steps.
"""

# ╔═╡ 6c54aec1-0c49-4f3f-979c-7fbf964ee22f
time_bound = AdaptiveTimeBoundary(1.0,10.0,1e-3)

# ╔═╡ a8a6e77e-1d13-4214-8927-0c6f1d1af222
ad_resp , ad_err,ad_time = generalized_alpha_adaptive(5.0,0.01, 1.0, 0.0, 1.0,time_bound)

# ╔═╡ a563b1e3-3013-4059-a33c-90c97f5b07fc
begin
	plot(ad_time.times,ad_resp.u[1,:],title="Displacement Response (Adaptive)",xlabel=L"time",ylabel=L"displacement",label=L"$\theta$")
	plot!(ad_time.times,ad_resp.u[2,:],label=L"u")

end

# ╔═╡ 34a2be2a-1256-471b-a6e3-7234702e796d
begin
	plot(ad_time.times,ad_err.e_abs,title="Absolute Error (Task 6)",xlabel=L"t",ylabel=L"||e||",label=L"$||e||$")
end

# ╔═╡ 4b3d6d84-07dc-4ccc-b26e-322946749c74
begin
	plot(ad_time.times,ad_err.η,title="Relative Error (Adaptive)",xlabel=L"t",ylabel=L"η",label=L"$\eta$")
end

# ╔═╡ 5050aadf-e39a-42aa-8640-6c0c2be826b3
begin
	plot(ad_time.times,ad_err.e_cum,title="Cumulative Error (Adaptive)",xlabel=L"t",ylabel=L"e_{cm}",label=L"$e_{cm}$")
end

# ╔═╡ 1cdbbff5-6d98-4b64-9a16-a4828daf05e8
begin
	plot(ad_time.times,ad_time.steps,title="Time step evolution (adaptive)",xlabel=L"time",ylabel=L"Δt",label=L"$Δt_{adaptive}$")
	plot!(time_a.times,time_a.steps,label=L"Δt_{3a}")
end

# ╔═╡ 23aa93b5-949e-4418-bfe0-b9033ff91627
begin
	plot(ad_time.times,ad_err.e_cum,title="Cumulative Error (Adaptive & Task 3a)",xlabel=L"t",ylabel=L"e_{cm}",label=L"$e_{cm, adaptive}$")
	plot!(time_a.times,err_a.e_cum,label=L"$e_{cm, 3a}$")
end

# ╔═╡ df7d4a60-cf6a-4a43-b2ec-dc795f2fdce0
md""" ##### Cumulative error after $5$ s:

"""

# ╔═╡ 3a8dabac-cfdd-49f1-a929-ad6899800bb0
(adapt = ad_err.e_cum[end],t3a =err_a.e_cum[end] )

# ╔═╡ 39d90f92-47dd-4046-8f99-1148d11381c1
md"""
!!! note
	From the previous result, we can observe that the cumulative error in the adaptive method is less than the fixed one and this was anticipated.
"""

# ╔═╡ 7a2e7714-ec5f-4ed5-8fd9-8061cce04678
md""" ##### Number of steps:

"""

# ╔═╡ e7c943c0-52eb-46a2-83ba-b2719458ab53
(adapt = length(ad_time.steps),t3a = length(time_a.steps) )

# ╔═╡ Cell order:
# ╠═e3dbdee8-7f62-4301-b848-6393954121a8
# ╠═c5bf3248-2ff1-4257-ba34-8b322b667abb
# ╟─6075d462-d1ff-4b30-9de8-2f7a1eb9ecdf
# ╠═774a153f-df2e-4fd0-a8db-ba370f8057c1
# ╠═0b86fb31-1e75-4c0d-8704-a95ec5318218
# ╠═380f1415-27af-4b21-a93e-9fbabdd27a30
# ╟─67cf6c46-7e1e-4f20-bba7-d6c2b9d53cf0
# ╠═45b7da6f-7b81-452c-8274-fa797698d35f
# ╠═b1057602-76fe-469c-b382-05c65fb14e60
# ╠═7b31aad9-d4f3-40d7-a7c7-18d25cb19139
# ╠═3cea61f0-bc04-4092-8e3b-7803932865d5
# ╠═9d7e83e5-4d75-4d96-913f-a656d9866735
# ╟─e0ced7bc-baae-4c15-8b6e-9e50be04ced9
# ╠═f27a727f-f80c-4dd2-bdef-63c7f3ca587f
# ╠═2a801b3c-49a3-4aa7-be75-452711c4d4a0
# ╟─20193589-4189-4347-ac5d-0e4370152ce6
# ╠═7cb8658f-9a03-4807-8e18-3775b3f675aa
# ╟─4ac90089-eee4-4ce3-a2e8-0152877f4414
# ╟─5309e3a3-ab35-4cc4-88cc-119813ab62a5
# ╟─e3ad5cc3-4843-4d80-a0ba-99278ca335be
# ╠═590a6c46-c9ee-4976-afdb-e873bd0af809
# ╠═f81d3f7a-3543-4cda-baad-198410bb9aef
# ╠═31592de4-c0f4-40b8-a940-06a1c83053c2
# ╠═8b058921-b41b-470b-b45a-35af6149428f
# ╠═85eda22b-ac0f-49ae-aaf2-2aca6a6e7070
# ╠═eb4a70b8-e07f-4311-8d5c-b1bb879ba726
# ╟─1286dcc1-d323-4a5f-8dbc-f874c2ad5377
# ╠═205bd7ff-ab41-4ce2-bec8-6bfa411f21d0
# ╠═e070700c-5bb9-4c62-9fe8-ac7b4eed7284
# ╠═0fdf106a-c40d-420d-a848-bf727fb49794
# ╠═6c54aec1-0c49-4f3f-979c-7fbf964ee22f
# ╠═a8a6e77e-1d13-4214-8927-0c6f1d1af222
# ╠═a563b1e3-3013-4059-a33c-90c97f5b07fc
# ╠═34a2be2a-1256-471b-a6e3-7234702e796d
# ╠═4b3d6d84-07dc-4ccc-b26e-322946749c74
# ╠═5050aadf-e39a-42aa-8640-6c0c2be826b3
# ╠═1cdbbff5-6d98-4b64-9a16-a4828daf05e8
# ╠═23aa93b5-949e-4418-bfe0-b9033ff91627
# ╟─df7d4a60-cf6a-4a43-b2ec-dc795f2fdce0
# ╠═3a8dabac-cfdd-49f1-a929-ad6899800bb0
# ╟─39d90f92-47dd-4046-8f99-1148d11381c1
# ╟─7a2e7714-ec5f-4ed5-8fd9-8061cce04678
# ╠═e7c943c0-52eb-46a2-83ba-b2719458ab53
