function [By2000,results1,results2,redchi2,fval] = cal_redchi2(Ac,Af,CC,ratio2,H,H_unm,F,F_y,F_unm,xs,ds,err,inv_nsigma,lb,ub,ub_ob,smpl,sigma,tp,T,lambdas,m,sq,sq_unm,no_sample)

rng('shuffle')
[BS2000,sd2000] = perturb_m(Ac,Af,xs,smpl,sigma,F,no_sample,tp);

for j = 1:size(BS2000,2)
    [F_bs,sq_bs] = coeff_gen_m(H,F,Ac,BS2000(:,j),smpl,tp,[16,41+15]);
    [By2000(:,j),fval(j),RSS,Ps,resid,preEPE(j)] = GLEAM_opt_direcFixed(H,H_unm,F_bs,F_unm,Ac,Af,CC,ratio2,lb,ub,ub_ob,m,sq_bs,sq_unm,inv_nsigma,smpl,sigma,lambdas,T,tp,ds,F,sq);
end

notopt = checkbootstrap(fval);
By2000(:,:,notopt) = [];
preEPE(notopt) = [];
By2000 = By2000(:,:,1:2000);
preEPE = preEPE(1:2000);
EPE = sum(preEPE)/length(preEPE);
df = length(find(F_y~=0))/2./sd2000.*(EPE - err);
redchi2 = length(find(F_y~=0)) .* EPE ./ df;
By2000 = squeeze(By2000);

end