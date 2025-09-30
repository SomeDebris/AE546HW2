module AE546HW2

using Plots
using LaTeXStrings
using Match

greet() = print("Hello World!")

@enum GuidanceLaw Persuit FixedLead ConstantBearing

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
        R_dot = V_T * cos(β - θ_T) - V_M * cos(β - θ)
        β_dot = -(V_T * sin(β - θ_T) - V_M * sin(β - θ)) / R

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

function p2a(guidance_law::GuidanceLaw;
        V_M = 300,
        V_T = 200,
        R₀ = 6000,
        β₀ = deg2rad(30),
        θ  = deg2rad(20),
        θ_T= 0,
        δt = 0.02,
        maxiter::Int = 100000,
        tol = 0.01)

    @debug "Integration parameters" maxiter tol

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

    # For GuidanceLaw.FixedLead
    Θ₀ = asin((V_T / V_M) * sin(β₀))
    if GuidanceLaw.FixedLead == guidance_law
        @debug "Fixed lead guidance law selected." Θ₀
    end

    while true
        iters += 1
        # Guidance law:
        @match guidance_law begin
            GuidanceLaw.Persuit   => begin
                θ = β
                θ_dot = β_dot
            end
            GuidanceLaw.FixedLead => begin
                θ = β - Θ₀
                θ_dot = β_dot
            end
            GuidanceLaw.ConstantBearing => begin
                θ = β - asin((V_T / V_M) * sin(β₀))
                θ_dot = β_dot - (V_T * cos(β) * β_dot)/(V_M * sqrt(1 - (V_T^2 * sin(β)^2)/(V_M^2)))
            end
        end

        R_dot = V_T * cos(β - θ_T) - V_M * cos(β - θ)
        β_dot = -(V_T * sin(β - θ_T) - V_M * sin(β - θ)) / R

        R += R_dot * δt
        β += β_dot * δt

        missile_position += (V_M * δt) .* [cos(θ), sin(θ)]
        target_position  += (V_T * δt) .* [cos(θ_T), sin(θ_T)]

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


    t = (iters - 1) * δt

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
    plot(missile_path,
        label = "Missile Path",
        xlabel= L"$x$ [m]",
        ylabel= L"$y$ [m]",
        #=aspect_ratio=:equal,=#
        gridstyle=:dash,
        xlims=(0,round_up_nearest(missile_position[1], 1000)),
        ylims=(0,round_up_nearest(missile_position[2], 1000)))
    plot!(target_path,
         label = "Target Path")
end



end # module AE546HW2
