using LinearAlgebra, Clustering
function MCA_ml_spectral_clustering(X,k; ell=1,sigma=0.1)
    n, m = size(X)

    g = vec(sum(X.^2, dims=2))
    G = g*ones(n)' + ones(n)*g' - 2*X*X'
    W = exp.(-G / (2*sigma^2))
    W = W - diagm(diag(W))             # 完全グラフで隣接行列を生成
    
    D = diagm(vec(sum(W, dims=2)))
    L = D - W
    
    ~, U = eigen(L,D)                  # 一般化固有値問題の求解
    y = kmeans(U[:,2:ell+1]',k)        # データ表現 U のクラスタリング

    return y.assignments
end