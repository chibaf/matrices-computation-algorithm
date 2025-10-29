using LinearAlgebra
function MCA_sqrtm_newton(A; p=2,X=copy(A),iter=20,tol=1e-12)
    Anorm = norm(A)
    for k = 1:iter
        E = (X - X^(1-p) * A) / p
        X = X - E

        if norm(E) / Anorm < tol  # 収束判定
            break
        end
    end

    return X
end