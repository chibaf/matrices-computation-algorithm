function MCA_expmv_cauchy(A,b; N=128,gamma=0,rho=norm(A))
    n = size(A,1)
    x = zeros(n)
    for j = 1:N/2                      # 数値積分
        thetaj = (2 * j - 1) * pi / N
        zj = gamma + rho * exp(im * thetaj)
        wj = exp(im * thetaj) / N
        y = (zj * I(n) - A) \ b        # 線形方程式の求解
        x = x + rho * exp(zj) * wj * y
    end

    return x = 2 * real(x)
end