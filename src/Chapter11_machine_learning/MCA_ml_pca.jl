using LinearAlgebra
function MCA_ml_pca(X,ell)
    n = size(X,1)
    Xt = X .- sum(X, dims=1) / n
    ~, ~, V = svd(Xt)        # 特異値分解

    return (B = V[:,1:ell], Y = X * V[:,1:ell])
end