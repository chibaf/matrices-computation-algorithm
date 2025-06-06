using LinearAlgebra, SparseArrays
function MCA_pde_poisson(L,N,f)
    dx = L / (N+1)
    A1 = spdiagm(0 => -2*ones(N), -1=>ones(N-1), 1=>ones(N-1)) / dx^2
    A2 = kron(I(N),A1) + kron(A1,I(N)) # クロネッカー積による行列生成
    u = A2 \ f                         # 線形方程式の求解
    
    return u = reshape(u,N,N)
end