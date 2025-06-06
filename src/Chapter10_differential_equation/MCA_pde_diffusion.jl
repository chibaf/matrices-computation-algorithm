using LinearAlgebra, SparseArrays
function MCA_pde_diffusion(alpha,L,u0,ux0,uxL,h,tmax)
    N = size(u0,1); dx = L / (N+1);
    A = spdiagm(0 => -2*ones(N), -1=>ones(N-1), 1=>ones(N-1)) / dx^2;
    b = zeros(N); b[1] = ux0; b[N] = uxL; b = b / dx^2;
    f(t) = alpha * b;

    # 陰的オイラー法による常微分方程式の求解
    u, ~ = MCA_ode_implicit_euler(alpha*A,f,u0,h,tmax);

    return u
end