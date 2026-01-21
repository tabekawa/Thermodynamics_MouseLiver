function [kcat, Km] = sample_kinparam_brenda(lnKeq,Km_r,kcat_mean_r,kcat_sd_r,Km_mean_r,Km_sd_r,numSamples)

rng('shuffle')
l = 1;
for i = 1:length(Km_r)
    if Km_r(i)<0
        HR_sign(l) = -1;
        l = l + 1;
    elseif Km_r(i)>0
        HR_sign(l) = 1;
        l = l + 1;
    end
end
Km_mean_r = Km_mean_r(Km_mean_r~=0);
Km_sd_r = Km_sd_r(Km_sd_r~=0);
try
    lb = [kcat_mean_r-2.*kcat_sd_r;Km_mean_r-2.*Km_sd_r];
    ub = [kcat_mean_r+2.*kcat_sd_r;Km_mean_r+2.*Km_sd_r];
    HR = [1, -1, HR_sign];
    P.A = [HR;diag(ones(2+nnz(Km_r~=0),1));diag(-ones(2+nnz(Km_r~=0),1))];
    P.b = [log10(exp(lnKeq));ub;-1.*lb];
    param = polySampler(P,length(lb)^2,numSamples);
catch
    disp('bounds are difined within 3sigma')
    lb = [kcat_mean_r-3.*kcat_sd_r;Km_mean_r-3.*Km_sd_r];
    ub = [kcat_mean_r+3.*kcat_sd_r;Km_mean_r+3.*Km_sd_r];
    HR = [1, -1, HR_sign];
    P.A = [HR;diag(ones(2+nnz(Km_r~=0),1));diag(-ones(2+nnz(Km_r~=0),1))];
    P.b = [log10(exp(lnKeq));ub;-1.*lb];
    param = polySampler(P,length(lb)^2,numSamples);
end

kcat = 10.^param(1:2,:);
Km_r = repmat(Km_r,1,1,numSamples);
for i = 1:size(Km_r,3)
    Km_i = squeeze(Km_r(:,:,i));
    Km_i(Km_i<0) = -10.^param(2+find(HR_sign==-1),i);
    Km_i(Km_i>0) = 10.^param(2+find(HR_sign==1),i);
    Km_r(:,:,i) = Km_i;
end
Km = Km_r;
end