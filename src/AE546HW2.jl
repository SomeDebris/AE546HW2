module AE546HW2

using Plots, ModelingToolkit, DifferentialEquations, LaTeXStrings

greet() = print("Hello World!")

struct GuidanceProblem
    # Missile velocity vector
    V_M::AbstractVector
    # Target velocity vector
    V_T::AbstractVector
    R₀
    β₀
    t₀
    δt
end

function GuidanceProblem(v_M::Real, v_T::Real, θ₀, θ_T₀, R₀, β₀, t₀, δt)
    V_M = v_M * [cos(θ₀), sin(θ₀), 0]
    V_T = v_T * [cos(θ_T₀), sin(θ_T₀), 0]

    GuidanceProblem(V_M, V_T, R₀, β₀, t₀, δt)
end

function R_diff1(Velocity_Target, Velocity_Missile, β, θ_Missile, θ_Target)
    Velocity_Target * cos(β - θ_Target) - Velocity_Missile * cos(β - θ_Missile)
end
function β_diff1(Velocity_Target, Velocity_Missile, R, β, θ_Missile, θ_Target)
    -(Velocity_Target * sin(β - θ_Target) - Velocity_Missile * sin(β - θ_Missile)) / R
end

function p1b()
    V_M = 300
    V_T = 200
    R0  = 6000
    β0  = deg2rad(30)
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

    R = R0
    β = β0
    R_dot = 0.0
    β_dot = 0.0

    # * KRIS!!!
    # * we have to kill the fucking titan
    # Shoot missile           ♡            Do nothing

    for k in range(1, N)
        R_dot = V_T * cos(β - θ_T) - V_M * cos(β - θ)
        β_dot = -(V_T * sin(β - θ_T) - V_M * sin(β - θ)) / R

        R = R + R_dot * δt
        β = β + β_dot * δt

        Rlog[k] = R
        βlog[k] = β
        R_dotlog[k] = R_dot
        β_dotlog[k] = β_dot

    end

    @info "Initial values (meters, degrees)" R=R0 β=rad2deg(β0) R_dot=R_dotlog[1] β_dot=rad2deg(β_dotlog[1])
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

function p2a()
    V_M = 300
    V_T = 200
    R0  = 6000
    β0  = deg2rad(30)
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

    # Thing I shall be recording this time around
    θ_dotlog = Matrix{Float64}(undef, N, 1)

    R = R0
    β = β0
    R_dot = 0.0
    β_dot = 0.0

    # * KRIS!!!
    # * we have to kill the fucking titan
    # Shoot missile           ♡            Do nothing

    for k in range(1, N)
        # Guidance law:
        θ = β

        R_dot = V_T * cos(β - θ_T) - V_M * cos(β - θ)
        β_dot = -(V_T * sin(β - θ_T) - V_M * sin(β - θ)) / R

        R = R + R_dot * δt
        β = β + β_dot * δt

        Rlog[k] = R
        βlog[k] = β
        R_dotlog[k] = R_dot
        β_dotlog[k] = β_dot

        # Guidance law also goes and says:
        θ_dotlog[k] = β_dot
    end

    @info "Initial values (meters, degrees)" R=R0 β=rad2deg(β0) R_dot=R_dotlog[1] β_dot=rad2deg(β_dotlog[1])
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

end # module AE546HW2
