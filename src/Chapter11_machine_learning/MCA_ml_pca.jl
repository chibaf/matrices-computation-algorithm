using LinearAlgebra
function MCA_ml_pca(X,l)
    n = size(X,1)
    xbar = sum(X, dims=2) / n
    Xt = X .- xbar

    ~, ~, V = svd(Xt)        # 特異値分解
    Y = X * V[:,1:l]

    return V[:,1:l], Y
end