module AE546HW2

using Plots
using LaTeXStrings
using Match
using Printf
using LinearAlgebra
#=using ModelingToolkit=#
#=using DifferentialEquations=#

greet() = print("Hello World!")

@enum GuidanceLaw Persuit FixedLead ConstantBearing ProportionalNavigation

struct GuidanceProblem
    # Missile velocity vector
    Velocity_missile::AbstractVector
    # Target velocity vector
    Velocity_target::AbstractVector
    R₀
    β₀
    t₀
    δt
end

function GuidanceProblem(v_M::Real, v_T::Real, θ₀, θ_T₀, R₀, β₀, t₀, δt)
    Velocity_missile = v_M * [cos(θ₀), sin(θ₀), 0]
    Velocity_target = v_T * [cos(θ_T₀), sin(θ_T₀), 0]

    GuidanceProblem(Velocity_missile, Velocity_target, R₀, β₀, t₀, δt)
end

function R_diff1(Velocity_Target, Velocity_Missile, β, θ_Missile, θ_Target)
    Velocity_Target * cos(β - θ_Target) - Velocity_Missile * cos(β - θ_Missile)
end
function β_diff1(Velocity_Target, Velocity_Missile, R, β, θ_Missile, θ_Target)
    -(Velocity_Target * sin(β - θ_Target) - Velocity_Missile * sin(β - θ_Missile)) / R
end

function p1b()
    Velocity_missile = 300
    Velocity_target = 200
    R₀  = 6000
    β₀  = deg2rad(30)
    θ   = deg2rad(20)
    θ_T = 0
    δt  = 0.02
    T   = 20
    N   = round(Int, T / δt)

    # Allocate memory. Do not assign value (increases speed, but cannot assume that each value is a 0)
    Rlog = Matrix{Float64}(undef, N, 1)
    βlog = Matrix{Float64}(undef, N, 1)
    R_dotlog = Matrix{Float64}(undef, N, 1)
    β_dotlog = Matrix{Float64}(undef, N, 1)

    R = R₀
    β = β₀
    R_dot = 0.0
    β_dot = 0.0

    # * KRIS!!!
    # * we have to kill the fucking titan
    # Shoot missile           ♡            Do nothing

    for k in range(1, N)
        R_dot = Velocity_target * cos(β - θ_T) - Velocity_missile * cos(β - θ)
        β_dot = -(Velocity_target * sin(β - θ_T) - Velocity_missile * sin(β - θ)) / R

        R = R + R_dot * δt
        β = β + β_dot * δt

        Rlog[k] = R
        βlog[k] = β
        R_dotlog[k] = R_dot
        β_dotlog[k] = β_dot

    end

    @info "Initial values (meters, degrees)" R=R₀ β=rad2deg(β₀) R_dot=R_dotlog[1] β_dot=rad2deg(β_dotlog[1])
    @info "Final values (meters, degrees)" R β=rad2deg(β) R_dot=R_dotlog[end] β_dot=rad2deg(β_dotlog[end])

    t = range(0, N - 1) * δt

    plot(t, [Rlog,rad2deg.(βlog)];
         layout = (2,1),
         xlabel=["" L"$t$ [s]"],
         ylabel=[L"$R$ [m]" L"$\beta$ [deg]"],
         xlims=(0, T),
         ylims=[(3500,6000) (15,30)],
         #=legend=false,=#
         gridstyle=:dash,
         label=[L"R(t)" L"\beta(t)"],
        )
end

round_up_nearest(x, step) = ceil(x / step) * step

function guidancesim(guidance_law::AE546HW2.GuidanceLaw;
        Velocity_missile = 300,
        Velocity_target = 200,
        R₀ = 6000,
        β₀ = deg2rad(30),
        θ  = deg2rad(20),
        θ_T= 0,
        δt = 0.02,
        maxiter::Int = 100000,
        tol = 0.01)

    @debug "Integration parameters" maxiter tol

    if !(guidance_law in (Persuit, FixedLead, ConstantBearing))
        @error "Implimentation for this guidance law does not exist! Implemented laws are AE546HW2.Persuit, AE546HW2.FixedLead, AE546HW2.ConstantBearing." guidance_law
        return 1
    end

    guidancelawtext = @match guidance_law begin
        $Persuit         => "guidance-persuit"
        $FixedLead       => "guidance-fixedlead"
        $ConstantBearing => "guidance-constbearing"
    end

    #=T   = 20=#
    #=N   = round(Int, T / δt)=#

    # Allocate memory. Do not assign value (increases speed, but cannot assume that each value is a 0)
    Rlog = Vector{Float64}(undef, 0)
    βlog = Vector{Float64}(undef, 0)
    R_dotlog = Vector{Float64}(undef, 0)
    β_dotlog = Vector{Float64}(undef, 0)

    # Save the x and y of the missile and target
    missile_path = Vector{Tuple{Float64,Float64}}(undef, 0)
    target_path = Vector{Tuple{Float64,Float64}}(undef, 0)

    # Thing I shall be recording this time around
    θ_dotlog = Vector{Float64}(undef, 0)

    R = R₀
    β = β₀
    R_dot = 0.0
    β_dot = 0.0
    θ_dot = 0.0
    missile_position = [0.0, 0.0]
    target_position =  missile_position .+ [R₀*cos(β₀), R₀*sin(β₀)]
    @debug "target position initial" target_position

    # * KRIS!!!
    # * we have to kill the fucking titan
    # Shoot missile           ♡            Do nothing

    iters = 0

    # For FixedLead
    Θ₀ = asin((Velocity_target / Velocity_missile) * sin(β₀))
    if FixedLead == guidance_law
        @debug "Fixed lead guidance law selected." Θ₀
    end

    while true
        iters += 1
        # Guidance law:
        θ = @match guidance_law begin
            $Persuit   => β
            $FixedLead => β - Θ₀
            $ConstantBearing => β - asin((Velocity_target / Velocity_missile) * sin(β))
        end

        R_dot = Velocity_target * cos(β - θ_T) - Velocity_missile * cos(β - θ)
        β_dot = -(Velocity_target * sin(β - θ_T) - Velocity_missile * sin(β - θ)) / R

        # Guidance law, θ_dot computed afterwards because R_dotβ_dot are here
        θ_dot = @match guidance_law begin
            $Persuit || $FixedLead => β_dot
            $ConstantBearing      => β_dot - (Velocity_target * cos(β) * β_dot)/(Velocity_missile * sqrt(1 - (Velocity_target^2 * sin(β)^2)/(Velocity_missile^2)))
        end

        R += R_dot * δt
        β += β_dot * δt

        missile_position += (Velocity_missile * δt) .* [cos(θ), sin(θ)]
        target_position  += (Velocity_target * δt) .* [cos(θ_T), sin(θ_T)]

        push!(Rlog, R)
        push!(βlog, β)
        push!(R_dotlog, R_dot)
        push!(β_dotlog, β_dot)

        push!(missile_path, Tuple{Float64,Float64}(missile_position))
        push!(target_path,  Tuple{Float64,Float64}(target_position))

        # Guidance law also goes and says:
        push!(θ_dotlog, θ_dot)

        if R <= tol
            @debug "R value less than tolerance" iters
            break
        end

        if iters >= maxiter
            @warn "Max iterations reached." iters
            break
        end
    end

    @debug "Missile path" missile_path
    @debug "Target path" target_path

    @debug "Theta and Thetadot" θ θ_dotlog


    t = (iters - 1) * δt
    trange = range(0, iters - 1) * δt

    @info "Initial values (meters, degrees)" R=R₀ β=rad2deg(β₀) R_dot=R_dotlog[1] β_dot=rad2deg(β_dotlog[1]) t=0

    @info "Final values (meters, degrees)" R β=rad2deg(β) R_dot=R_dotlog[end] β_dot=rad2deg(β_dotlog[end]) t

    #=plot(t, [Rlog,rad2deg.(βlog)];=#
    #=     layout = (2,1),=#
    #=     xlabel=["" L"$t$ [s]"],=#
    #=     ylabel=[L"$R$ [m]" L"$\beta$ [deg]"],=#
    #=     xlims=(0, round_up_nearest(iters * δt, 10)),=#
    #=     ylims=[(0,R₀) (0,round_up_nearest(rad2deg(β₀), 5))],=#
    #=     #=legend=false,=#=#
    #=     gridstyle=:dash,=#
    #=     label=[L"R(t)" L"\beta(t)"],=#
    #=    )=#
    pxy = plot(missile_path,
        label = "Missile Path",
        xlabel= L"$x$ [m]",
        ylabel= L"$y$ [m]",
        #=aspect_ratio=:equal,=#
        gridstyle=:dash,
        xlims=(0,round_up_nearest(missile_position[1], 1000)),
        ylims=(0,round_up_nearest(missile_position[2], 1000)))
    plot!(target_path,
         label = "Target Path")
    #=plot(θ_dotlog)=#
    pθ_dot = plot(trange, θ_dotlog;
                xlims  = (0,round_up_nearest(t, 10)),
                label  = L"\dot{\theta}(t)",
                xlabel = L"$t$ [s]",
                gridstyle=:dash,
                ylabel = L"$\theta$ [rad]")

    savefig(pxy,    @sprintf "plot_path_%s_.pdf" guidancelawtext)
    savefig(pθ_dot, @sprintf "plot_thetadot_%s_.pdf" guidancelawtext)

    return 0
end

struct PNGuidanceIC2D
    position_missile::Tuple{Real,Real}
    position_target::Tuple{Real,Real}
    velocity_missile::Real
    velocity_target::Real
    θ_missile::Real
    θ_target::Real
    λ::Real
    T::Real
    hit_tolerance::Real
    target_acceleration::Real
end



function PNtraj(ic::PNGuidanceIC2D; maxiter::UInt=10000, δt = 0.1)
    pos_missile = collect(ic.position_missile)
    pos_target  = collect(ic.position_target)

    trajectory_missile = [ic.position_missile]
    trajectory_target  = [ic.position_target]

    β = 0.0
    β_diff1 = 0.0
    θ_missile = ic.θ_missile
    θ_missile_diff1 = 0.0

    β_lastiter = 0.0
    θ_missile_lastiter = 0.0

    θ_missile_command = 0.0
    θ_missile_diff1_command = 0.0

    velocity_target = ic.velocity_target

    iterations::UInt = 0
    while true
        iterations += 1

        β = atan(pos_target[2] - pos_missile[2],  pos_target[1] - pos_missile[1])
        β_diff1 = (β - β_lastiter) / δt

        θ_missile_diff1_command = ic.λ * β_diff1

        if ic.T != 0
            θ_missile_command += θ_missile_diff1_command * δt
            θ_missile += ((θ_missile_command - θ_missile)/ic.T) * δt
        else
            θ_missile += θ_missile_diff1_command * δt
        end

        pos_missile += ic.velocity_missile .* [cos(θ_missile), sin(θ_missile)] .* δt
        pos_target  += velocity_target  .* [cos(ic.θ_target), sin(ic.θ_target)] .* δt

        velocity_target += ic.target_acceleration * δt

        push!(trajectory_missile, Tuple(pos_missile))
        push!(trajectory_target,  Tuple(pos_target))

        β_lastiter = β

        if norm(pos_target - pos_missile) < ic.hit_tolerance
            @debug "Breaking PN guidance loop: target reached" missile_target_distance=norm(pos_target - pos_missile) ic.hit_tolerance pos_target pos_missile
            break
        end

        if iterations >= maxiter
            @error "Maximum iteratins reached." maxiter
            break
        end
    end

    time_final = iterations * δt

    trajectory_missile, trajectory_target, time_final
end


#=function missdist(T, N; maxiter::UInt=10000, δt = 0.1, n_T = 0.0, θ₀ = deg2rad(3))=#
#=    x = [0, 0]=#
#=    x_diff = x=#
#==#
#=    θ = θ₀=#
#==#
#=    λVm = N * V_C=#
#==#
#=    y_t = velocit=#
#==#
#=    #=pos_target = collect(pos_target_ic)=#=#
#==#
#=    # Assumptions: V_C, V_M do not change.=#
#==#
#=    # this value shall be used as N in the function.=#
#=    # thought process: we are given an initial value for θ₀. Thus, we can compute the value of the terms that comprise N.=#
#==#
#==#
#=    iterations:UInt = 0=#
#=    while true=#
#=        iterations += 1=#
#==#
#=        pos_target  += velocity_target  .* [cos(ic.θ_target), sin(ic.θ_target)] .* δt=#
#==#
#=        x_diff = [ x[2],=#
#=                  -x[2]/T - (λVm / T) * atan(x[1] - y_t=#
#=end=#


function computeyandacfunctions(N;
    Velocity_closing = 120,
    Velocity_missile = 300,
    θ₀ = deg2rad(5),
    R₀ = 6000,
    n_T = 0,
    samples::Int = 3000)

    @debug "values" N n_T

    # Compute the values of λ for this problem
    λ = @. N * Velocity_closing / Velocity_missile

    # Final time: eq 3.26
    time_final = R₀/Velocity_closing

    # y(t), eq 3.41
    y(t) = @match (n_T, N) begin
        # n_T == 0  AND  N != 1
        (nt, Nset) where (nt == 0 && Nset != 1)  =>  begin
            @. (Velocity_missile * θ₀)/(N - 1) * (1 - 1/(time_final^(N-1)) * (time_final - t)^(N-1))
        end

        # n_T != 0  AND  N is 1
        (nt, Nset) where (nt != 0 && Nset == 1) => begin
            @. -n_T * (time_final - t) * (t - time_final*log(time_final - t) + time_final * log(time_final))
        end

        (nt, Nset) where (nt != 0 && Nset == 2) => begin
            @. -n_T * ((1/time_final) * (time_final - t)^2 - (time_final - t))
        end

        # n_T != 0  AND  N is neither 1 nor 2
        (nt, Nset) where (nt != 0 && !(Nset in (1,2))) => begin
            @. n_T * (time_final - t)^N * ( (time_final - t)^(2 - N)/(2 - N) - time_final * (time_final - t)^(1-N)/(1-N) + (time_final^(2-N))/((1-N) * (2-N)))
        end
    end

    a_c(t) = @match (n_T, N) begin
        (nt, Nset) where nt == 0 => @. - (Velocity_missile * θ₀)/(time_final^(N-1)) * (time_final - t)^(N-2)
        (nt, Nset) where (nt != 0  && Nset == 1)  => @. Inf
        # NOTE: the "+ (t * 0)" at the end makes julia output a series instead of a single value
        (nt, Nset) where (nt != 0  && Nset == 2)  => @. -(1 + 2/time_final) * n_T + (t * 0)
        (nt, Nset) where (nt != 0  && !(Nset in (1,2)))  => @. (n_T * N) / (2 - N) * (1 - (1 - t/time_final)^(N-2))
    end

    y, a_c, time_final
end

function problem3at_ntsetting(N_values, n_T::Float64)
    funcs = computeyandacfunctions.(N_values; n_T=n_T)

    yPlot   = plot(; xlabel=L"$t$ [s]", ylabel=L"$y(t)$ [m]", gridstyle=:dash)
    a_cPlot = plot(; xlabel=L"$t$ [s]", ylabel=L"$a_c(t)$ [m/s]", gridstyle=:dash)

    dashpatterncyle = [:solid, :dashdotdot, :dash, :dashdot, :dot]

    for (idx, N) in enumerate(N_values)
        tbounds = (0, funcs[idx][3])
        @debug "The N value being passed to the plotting functions is:" N
        problem3plot_y!(yPlot, funcs[idx][1], tbounds, 500; label=L"N = %$(N)", ls=dashpatterncyle[((idx-1) % length(dashpatterncyle)) + 1])
        problem3plot_a_c!(a_cPlot, funcs[idx][2], tbounds, 500; label=L"N = %$(N)", ls=dashpatterncyle[((idx-1) % length(dashpatterncyle)) + 1])
    end

    savefig(yPlot, @sprintf "plot_p3-y_nT%g.pdf" n_T)
    savefig(a_cPlot, @sprintf "plot_p3-a_c_nT%g.pdf" n_T)
end



function problem3plot_y!(plotref::Plots.Plot, y::Function, time_bounds::Tuple{Real,Real}, samples::Int; kwargs...)
    trange = LinRange(time_bounds[1], time_bounds[2], samples)
    @debug "Adding y(t) to plot" plotref y time_bounds trange

    plot!(plotref,
          trange, y(trange);
          kwargs...)
          #=legend = false,=#
          #=xlabel=L"$t$ [s]",=#
          #=ylabel=L"$y$ [m]",=#
          #=gridstyle=:dash,=#
          #=xlims=(0,round_up_nearest(time_final,10)+1),=#
          #=ylims=(0,round_up_nearest(y(time_final), 10)))=#
end

function problem3plot_a_c!(plotref::Plots.Plot, a_c::Function, time_bounds::Tuple{Real,Real}, samples::Int; kwargs...)
    trange = LinRange(time_bounds[1], time_bounds[2], samples)
    @debug "Adding a_c(t) to plot" plotref a_c time_bounds trange

    plot!(plotref,
          trange, a_c(trange);
          kwargs...)
          #=legend = false,=#
          #=xlabel=L"$t$ [s]",=#
          #=ylabel=L"$y$ [m]",=#
          #=gridstyle=:dash,=#
          #=xlims=(0,round_up_nearest(time_final,10)+1),=#
          #=ylims=(0,round_up_nearest(y(time_final), 10)))=#
end

struct PNGuidanceIC3D
    VM
    VT
    N
    δt
    t_max
    tolR
    rm
    rt
    ψT
    γT
end

function PN3Dsim(ic::PNGuidanceIC3D, γ_d=false)
    missile_position = collect(ic.rm)
    target_position  = collect(ic.rt)

    vhatT = [cos(ic.γT) * sin(ic.ψT),
             cos(ic.γT) * cos(ic.ψT),
             sin(ic.γT)]

    velocity_target = ic.VT * vhatT

    R₀ = target_position - missile_position
    dist_xy₀ = hypot(R₀[1], R₀[2])
    ψM₀   = atan(R₀[1], R₀[2])
    γM₀   = atan(R₀[3], dist_xy₀)

    ψ     = ψM₀
    γ     = γM₀
    # number of datapoints:
    K = Int(ceil(ic.t_max/ic.δt)) + 1
    @debug "What is the value of K?" K ic.t_max ic.δt
    # data accumulators:
    RM = Vector{Tuple{Float64,Float64,Float64}}(undef, K)
    RT = Vector{Tuple{Float64,Float64,Float64}}(undef, K)

    Ψ = Vector{Float64}(undef, K)
    Γ = Vector{Float64}(undef, K)

    Ψ_L = Vector{Float64}(undef, K)
    Γ_L = Vector{Float64}(undef, K)

    tlog = Vector{Float64}(undef, K)

    azimuthElevation(v) = (atan(v[1], v[2]),
                           atan(v[3], hypot(v[1], v[2])))

    ψL_lastiter, γL_lastiter = azimuthElevation(R₀)

    iters = 0;

    for t in 0:ic.δt:ic.t_max
        iters += 1
        Rvec = target_position - missile_position
        R    = norm(Rvec)

        if R < ic.tolR
            @debug "R is lesser than tolR. breaking loop." R ic.tolR
            break
        end

        ψL, γL = azimuthElevation(Rvec)

        dψL = angle(exp(1im*ψL) / exp(1im*ψL_lastiter)) / ic.δt
        dγL = (γL - γL_lastiter) / ic.δt

        # PN commands on missile yaw/elevation
        dψ_cmd = ic.N * dψL
        dγ_cmd = @match γ_d begin
            (γd) where γd == false  =>  ic.N * dγL
            Real => ic.N * dγL + (γ_d - γL)
        end


        # update missile pointing
        ψ += dψ_cmd * ic.δt
        γ += dγ_cmd * ic.δt

        # missile and target velocity vectors
        vhatM = [cos(γ) * sin(ψ),
                 cos(γ) * cos(ψ),
                 sin(γ)]
        vm = ic.VM * vhatM

        # integrate positions
        missile_position += vm * ic.δt
        target_position  += velocity_target * ic.δt

        # Log
        RM[iters] = Tuple(missile_position)
        RT[iters] = Tuple(target_position)
        Ψ[iters] = ψ
        Γ[iters] = γ
        Ψ_L[iters] = ψL
        Γ_L[iters] = γL
        tlog[iters] = t

        # prepare next tick
        ψL_lastiter = ψL
        γL_lastiter = γL
    end

    send_values = range(1,iters-1)

    (iters = iters,
     RM = RM[send_values],
     RT = RT[send_values],
     Ψ = Ψ[send_values],
     Γ = Γ[send_values],
     Ψ_L = Ψ_L[send_values],
     Γ_L = Γ_L[send_values],
     tlog = tlog[send_values])
end

end # module AE546HW2
