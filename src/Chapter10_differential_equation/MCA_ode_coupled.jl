using LinearAlgebra
function MCA_ode_coupled(m,k,x0,t)
    n = size(m,1)
    M = diagm(m)
    K = diagm(0 => (k[1:n]+k[2:n+1]), 1=>-k[2:n], -1=>-k[2:n])
    
    omega2, U_tilde = eigen(K,M)            # 一般化固有値問題の求解
    
    c = U_tilde' * (M * x0)
    omega = sqrt.(omega2)
    
    # x(t)の計算（工夫するとループ無しで1行で書くことができる）
    xt = U_tilde*(cos.(omega.*t') .* c)
    
    return xt
end