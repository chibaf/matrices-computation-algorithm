using LinearAlgebra, SparseArrays
function MCA_pde_wave(L,c,u,v,h,tmax)
    N = size(u,1); dx = L / (N+1)
    kmax = Int(round(tmax/h))
    uvec = zeros(kmax,N); vvec = zeros(kmax,N); tvec = zeros(kmax)
    A = spdiagm(0 => -2*ones(N), -1=>ones(N-1), 1=>ones(N-1)) / dx^2

    # 陰的オイラー法
    Ah = h * c^2 * A
    F = lu(I(N) - h * Ah)         # 反復の外でLU分解を実施
    for k = 1:kmax
        t = k*h
        v = F \ (Ah * u + v)      # 線形方程式の求解（LU分解を繰り返し利用）
        u = u + h * v
        uvec[k,:] = u
        vvec[k,:] = v
        tvec[k] = t
    end

    return uvec, vvec, tvec
end