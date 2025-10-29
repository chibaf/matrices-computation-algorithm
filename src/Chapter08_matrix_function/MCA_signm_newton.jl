using LinearAlgebra
function MCA_signm_newton(A; X=copy(A),iter=20,tol=1e-14)
    Anorm = norm(A)
    for k = 1:iter
        E = (X - inv(X)) / 2
        X = X - E

        if norm(E) / Anorm < tol  # 収束判定
            break
        end
    end

    return X
end