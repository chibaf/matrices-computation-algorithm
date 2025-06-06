using Images, LinearAlgebra
function MCA_imgcomp_svd(img,k)
    A = Matrix{Float64}(img)           # 倍精度の行列に変換（各要素は[0,1]）
    U, S, V = svd(A)                   # 特異値分解に基づく低ランク近似
    Uk = U[:,1:k]; Sk = S[1:k]; Vk = V[:,1:k]
    Uk = Matrix{Float16}(Uk)           # 半精度の行列に変換
    Wk = Matrix{Float16}(Vk*diagm(Sk)) # 半精度の行列に変換
    return (Uk=Uk, Wk=Wk)
end
function MCA_imgcomp_svd_reconstruct(svdk)
    Ak = svdk.Uk * svdk.Wk'
    Ak = min.(1,max.(0,Ak))            # 各要素を[0,1]に揃える
    imgk = colorview(Gray,Ak)          # 近似画像の復元

    return imgk
end