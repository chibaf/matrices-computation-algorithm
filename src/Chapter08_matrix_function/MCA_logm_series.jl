using LinearAlgebra
function MCA_logm_series(A; r=20,m=4)
    S = copy(A)
    for i = 1:m    # A^1/s (s=2^m) の計算
        S = MCA_sqrtm_newton(S)
    end
    X = Ak = As = I(size(A,1)) - S
    for k = 2:r
        Ak = Ak * As
        X = X + Ak / k
    end
    X = -2^m * X

    return X
end