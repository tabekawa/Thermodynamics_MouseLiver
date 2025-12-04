function [H, F, inv_nsigma, sq, H_unm, F_unm, sq_unm] = create_obj_func(lnc_exp,sigma_c,G_exp_t,sigma_G_t,lb,ub,time_interval,tp,unm_mes)

%% 二次変数項の係数行列: H
% 濃度または生成自由エネルギー変数の二次項の係数行列: H_y
sigmas = [sigma_c;sigma_G_t];
H_y = 2 ./ (sigmas.^2);
H_y(sigmas==0) = 0;
H_y(1:length(lnc_exp)) = H_y(1:length(lnc_exp)) ./ (length(lnc_exp(lnc_exp~=100000000))-1*tp);
H_y(length(lnc_exp)+1:end) = H_y(length(lnc_exp)+1:end) ./ length(G_exp_t(G_exp_t~=100000000));
% 系統誤差変数の二次項の係数行列: H_a
H_a = zeros(length(lnc_exp)/tp+length(G_exp_t),1);
for i = 1:length(lnc_exp)/tp
    for j = 1:tp
        H_a(i) = H_a(i) + H_y(i+length(lnc_exp)/tp*(j-1));
    end
end
H_a(length(lnc_exp)/tp+1:end) = H_y(length(lnc_exp)+1:end);
% 濃度または生成自由エネルギー変数と系統誤差変数の二次項の係数行列: H_side
H_side = [];
for i = 1:tp
    H_side = [H_side, diag(H_y(1+length(lnc_exp)/tp*(i-1):length(lnc_exp)/tp*i))];
end
H_side = blkdiag(H_side, diag(H_y(length(lnc_exp)+1:end)));

H_y = diag(H_y);
H_a = diag(H_a);
H = [H_y, H_side'; H_side, H_a];

%% 一次変数項の係数ベクトル: F
% 濃度または生成自由エネルギー変数の一次項の係数ベクトル: F_y
F_y = -1 .* [lnc_exp;G_exp_t] .* diag(H_y);
% 系統誤差変数の一次項の係数ベクトル: F_a
F_a = zeros(length(lnc_exp)/tp+length(G_exp_t),1);
for i = 1:length(lnc_exp)/tp
    for j = 1:tp
        F_a(i) = F_a(i) + F_y(i+length(lnc_exp)/tp*(j-1));
    end
end
F_a(length(lnc_exp)/tp+1:end) = F_y(length(lnc_exp)+1:end);

F = [F_y; F_a];

%% 系統誤差絶対値変数の一次項の係数ベクトル: inv_nsigma
% 濃度の時点圧縮前での、系統誤差絶対値変数の一次項の係数ベクトル: inv_nsigma_y
inv_nsigma_y = diag(H_y) .* sigmas ./ 2;

inv_nsigma = zeros(length(lnc_exp)/tp+length(G_exp_t),1);
for i = 1:length(lnc_exp)/tp
    for j = 1:tp
        inv_nsigma(i) = inv_nsigma(i) + inv_nsigma_y(i+length(lnc_exp)/tp*(j-1));
    end
end
inv_nsigma(length(lnc_exp)/tp+1:end) = inv_nsigma_y(length(lnc_exp)+1:end);

%% 定数項: double
sq = [lnc_exp;G_exp_t]' * (H_y./2) * [lnc_exp;G_exp_t];

%% 未計測濃度の時間変化および中間濃度正則化項の二次項の係数行列: H_unm
% 同時点濃度変数（未計測時点）の二次項の係数ベクトル: H_unm_sametp
H_unm_sametp = zeros(length(lnc_exp),1);
H_unm_sametp(lnc_exp==100000000) = 1;
for i = 1:tp
    n_fluct(i) = 1 / time_interval(i) + 1 / time_interval(i+1);
end
n_fluct([1,7,8,14]) = n_fluct([1,7,8,14]) - 1;
for i = 1:tp
    H_unm_sametp(1+length(lnc_exp)/tp*(i-1):length(lnc_exp)/tp*i) = H_unm_sametp(1+length(lnc_exp)/tp*(i-1):length(lnc_exp)/tp*i) .* n_fluct(i);
end
% H_unm_sametpに同時点濃度変数（未計測時点に隣接した計測時点）の二次項の係数を追加
% ob/obの最初時点またはWTの最終時点が該当する場合は注意、1を引く。今回のデータでは該当なし。
for i = unm_mes
    if i > length(lnc_exp)/tp
        if lnc_exp(i-length(lnc_exp)/tp) == 100000000
            H_unm_sametp(i) = H_unm_sametp(i) + 1 / time_interval(fix(i/length(lnc_exp)*tp)+1);
        end
    end
    if i < length(lnc_exp)*(tp-1)/tp+1
        if lnc_exp(i+length(lnc_exp)/tp) == 100000000
            H_unm_sametp(i) = H_unm_sametp(i) + 1 / time_interval(fix(i/length(lnc_exp)*tp)+2);
        end
    end
end
% 隣接時点濃度変数の二次項の係数ベクトル: H_unm_side
H_unm_side = zeros(length(lnc_exp)*(tp-1)/tp,1);
H_unm_side(lnc_exp(1:length(lnc_exp)*(tp-1)/tp)==100000000 | lnc_exp(length(lnc_exp)/tp+1:end)==100000000) = 1;
n_fluct_side = 1 ./ time_interval(2:end-1);
n_fluct_side(7) = 0;
for i = 1:tp-1
    H_unm_side(1+length(lnc_exp)/tp*(i-1):length(lnc_exp)/tp*i) = -1 .* H_unm_side(1+length(lnc_exp)/tp*(i-1):length(lnc_exp)/tp*i) .* n_fluct_side(i);
end
% 中間濃度正則化項における濃度変数の二次項の係数ベクトル: H_med
H_med = zeros(length(lnc_exp),1);
H_med(lnc_exp==100000000) = 10^(-3) * 1;
%    H_unm = diag(H_med) .* 2;

H_unm = diag(H_unm_sametp);
for i = 1:length(H_unm_side)
    H_unm(length(lnc_exp)/tp+i,i) = H_unm_side(i);
    H_unm(i,length(lnc_exp)/tp+i) = H_unm_side(i);
end
H_unm = H_unm .* 2 .* 10^4;
H_unm = (H_unm + diag(H_med)) .* 2;

%% 中間濃度正則化項における濃度変数の一次項の係数ベクトル: F_unm
med_mets = mean([lb,ub],2);
F_unm = -2 .* repmat(med_mets,tp,1) .* H_med;

%% 中間濃度正則化項における定数項: double_unm
sq_unm = repmat(med_mets,tp,1)' * diag(H_med) * repmat(med_mets,tp,1);
end