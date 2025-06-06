using Images, LinearAlgebra
function MCA_imgcomp_dft(img,fs)
    A = Matrix{Float32}(img)           # 単精度の行列に変換（各要素は[0,1]）
    F = MCA_dft2(A,1)                  # 2次元離散フーリエ変換
    Ft = MCA_imgcomp_dft_filter(F,fs)  # フィルタリング
    return Ft = Matrix{ComplexF32}(Ft)
end
function MCA_imgcomp_dft_reconstruct(Ft,m,n,fs)
    F = Matrix{ComplexF32}(zeros(m,n))
    F[1:fs,1:fs] = Ft[1:fs,1:fs]
    F[1:fs,n-fs+1:n] = Ft[1:fs,fs+1:2*fs]
    F[m-fs+1:m,1:fs] = Ft[fs+1:2*fs,1:fs]
    F[m-fs+1:m,n-fs+1:n] = Ft[fs+1:2*fs,fs+1:2*fs]
    At = MCA_dft2(F,-1) / (m * n)
    imgt = colorview(Gray,real(At))
    return imgt
end
function MCA_dft2(A,pm)
    m, n = size(A)
    a = pm*2*im*pi/m
    W = exp.(a.*((0:m-1).*(0:m-1)'))
    F = W * A
    a = pm*2*im*pi/n
    W = exp.(a.*((0:n-1).*(0:n-1)'))
    F = F * W'
    return F
end
function MCA_imgcomp_dft_filter(F,fs)
    m, n = size(F);
    Ft = zeros(2*fs,2*fs)
    F11 = F[1:fs,1:fs]
    F12 = F[1:fs,n-fsize+1:n]
    F21 = F[m-fs+1:m,1:fs]
    F22 = F[m-fs+1:m,n-fs+1:n]
    Ft = [F11 F12; F21 F22]
    return Ft
end