function sor(A, b, ω, tol = 1e-10, maxiter = 100000)
    n = size(A,2)
    assert(size(A)==(n,n))

    ϕ = zeros(n)  # initial guess
    iter = 0

    while (norm(b - A*ϕ, 2) > tol) && (iter <= maxiter)

        if iter == maxiter
            println("Maximum number of iterations reached: $(iter)")
            return ϕ
        end

        for i = 1:n

            σ = 0
            for j = 1:n
                if j != i
                   σ += A[i, j]*ϕ[j]
                end
            end

            # ALTERNATIVE IMPLEMENTATION
            # ϕ[i] .= (1-ω)*ϕ[i] + ((ω/A[i,i]) * (b[i] - σ))

            # saves one multiplication
            ϕ[i] .+= ω*((b[i] - σ)/A[i,i] - ϕ[i])

        end

        iter += 1
    end

    return ϕ, iter
end
