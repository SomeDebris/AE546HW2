module AE546HW2

using Plots, ModelingToolkit, DifferentialEquations

greet() = print("Hello World!")

function p1b()
    V_M = 300
    V_T = 200
    R   = 6000
    β   = deg2rad(30)
    θ   = deg2rad(20)
    θ_T = 0
    δt  = 0.02
    T   = 20
    N   = round(Int, T / δt)

    # Allocate memory. Do not assign value (increases speed, but cannot assume that each value is a 0)
    Rlog = Matrix{Float64}(undef, 1, N)
    βlog = Matrix{Float64}(undef, 1, N)

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
    end

    t = range(0, N - 1) * δt
end

end # module AE546HW2
