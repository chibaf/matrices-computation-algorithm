using LinearAlgebra
function MCA_trigm_series(A; r=20,m=4)
    s = 2^m
    ck = im/s
    Ak = copy(A) 
    X = I(size(A,1)) + ck * Ak
    for k = 2:r
        ck = im * ck / (k*s)
        Ak = Ak * A
        X = X + ck * Ak
    end
    for k = 1:m         # X^s (s=2^m) の計算
        X = X * X
    end

    return (C = real(X), S = imag(X))
end