using LinearAlgebra
function MCA_pde_schrodinger(L,N,V,hm)
    dx = L / (N+1)
    A1 = diagm(0 => -2*ones(N), -1=>ones(N-1), 1=>ones(N-1)) / dx^2

    E,U = eigen(-hm*A1+diagm(V))  # 固有値問題の求解
    
    return E, U
end