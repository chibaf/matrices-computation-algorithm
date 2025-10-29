using LinearAlgebra
function MCA_expm_series(A; r=20,m=4)
    s = 2^m
    ck = 1/s
    Ak = copy(A) 
    X = I(size(A,1)) + ck * Ak
    for k = 2:r
        ck = ck / (k*s)
        Ak = Ak * A
        X = X + ck * Ak
    end
    for k = 1:m         # X^s (s=2^m) の計算
        X = X * X
    end

    return X
end