function sor_self(A, b, ω, tol = 1e-6, maxiter = 10000)
    n = size(A,2)
    assert(size(A)==(n,n))

    ϕ = b  # initial guess
    iter = 0

    while (norm(b - A*ϕ, 2) > tol) && (iter <= maxiter)

        if iter == maxiter
            println("Maximum number of iterations reached: $(iter)")
            return ϕ
        end

        for i = 1:n

            σ = 0
            for j = 1:n
                if i != j
                   σ += A[i, j]*ϕ[j]
                end
            end

            ϕ[i] .= (1-ω)*ϕ[i] + ((ω/A[i,i]) * (b[i] - σ))
            #ϕ[i] .+= ω*((b[i] - σ)/A[i,i] - ϕ[i])   # saves one multiplication

        end

        iter += 1
    end

    #println("Converged in $(iter) iterations")
    return ϕ
end

# params
n = 5
ω = 1.25
tol = 1e-6
println("Running for n = $n and ω = $ω")

A = rand(Float64, n, n)
b = rand(Float64, n)
x = A \ b
ϕ = sor_self(A, b, ω)

println("x ≈ ϕ: $(x ≈ ϕ)")
