mass(U, q, w, v, x, t, μ, η) = [1.0, 1.0, 1.0]
source(U, q, w, v, x, t, μ, η) = [0.0, 0.0, 0.0]
ubou(U, q, w, v, x, t, μ, η, û, n, τ) = [0.0, 0.0, 0.0]

function flux(U, q, w, v, x, t, μ, η)
    h, hu, hv = U

    h⁻¹ = 1/h

    uv = hu * h⁻¹ # u
    vv = hv * h⁻¹ # v

    g = μ[1] # gravitational acceleration

    p = g * h^2 / 2 # g * h^2 / 2

    f = [hu, hu * uv + p, hv * uv, hv, hu * vv, hv * vv + p]

    f = reshape(f, (3, 2))

    return f
end

function fbou(U, q, w, v, X, t, μ, η, Û, n, τ)
    f = flux(U, q, w, v, X, t, μ, η)
    f̂ = f[:, 1] * n[1] + f[:, 2] * n[2] + τ[1] * (U - Û)
    return f̂
end

function initu(X, μ, η)
    ϵ = 0.1 # Perturbation magnitude
    l = 0.5 # Gaussian width
    k = 0.5 # Sinusoidal wavenumber
    
    x = X[1]
    y = X[2]
    
    # The Bickley jet
    U = 1 / cosh(y)^2

    # Slightly off-center vortical perturbations
    ψ′ = exp(-(y + l/10)^2 / 2l^2) * cos(k * x) * cos(k * y)

    # Vortical velocity fields (ũ, ṽ) = (-∂_y, +∂_x) ψ̃
    u′ = + ψ′ * (k * tan(k * y) + y / l^2)
    v′ = - ψ′ * k * tan(k * x)

    Uᵢ₁ = 1.0        # h
    Uᵢ₂ = U + ϵ * u′ # h * u
    Uᵢ₃ = ϵ * v′     # h * v
   
    return [Uᵢ₁, Uᵢ₂, Uᵢ₃]
end
