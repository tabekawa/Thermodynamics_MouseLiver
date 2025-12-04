function [X,F,RSS,Ps,resid,err] = opt_get_sol(H_meas,H_unm,f_meas,f_unm,Ac,Af,C,ratio,lb,ub,ub_ob,sq,sq_unm,inv_nsigma,inv_nsigma_unm,smpl,sigma,lambda,T,tp,d,Acs,Cs,Ss,f_exp,sq_exp)

rG_lim = repmat(-0.01,size(Ac,1)*tp,1);
for r = 0:tp-1
    rG_lim(26+size(Ac,1)*r) = -0.01 - 2.578730581020227*log(10^(d(26+size(Ac,1)*r)*(7.2-7.4)));
    rG_lim(29+24+size(Ac,1)*r) = -0.01 - 2.578730581020227*log(10^(d(29+24+size(Ac,1)*r)*(7.2-7.4)));
end

[Z_meas, sol_meas] = opt_measured(H_meas,f_meas,Ac,Af,C,ratio,lb,ub,ub_ob,m,sq,inv_nsigma,lambda,T,tp,Acs,Cs,Ss,rG_lim,d);
if Z_meas < 10
    % 未計測代謝物推定
    f_conc = f(1:size(Ac,2)*tp);
    lb_meas_fixed = [repmat(lb,14,1);sol_meas(length(f_conc)+1:length(f_conc)+size(Af,2)+length(inv_nsigma))];
    ub_meas_fixed = [repmat(ub,14/2,1);repmat(ub_ob,14/2,1);sol_meas(length(f_conc)+1:length(f_conc)+size(Af,2)+length(inv_nsigma))];
    lb_meas_fixed(f_conc~=0) = sol_meas(f_conc~=0);
    ub_meas_fixed(f_conc~=0) = sol_meas(f_conc~=0);
    [Z_unm, sol_unm] = opt_unmeasured(H_unm,f_unm,Ac,Af,C,ratio,lb_meas_fixed,ub_meas_fixed,m,sq_unm,inv_nsigma_unm,lambda,T,tp,Acs,Cs,Ss,rG_lim,d);
    Z = Z_meas + Z_unm*10^(-10);
elseif Z_meas == 10
    Z = 10;
end

X = sol_unm;
F = Z;
err = X(1:size(f_exp,1))' * (H_meas./2) * X(1:size(f_exp,1)) + f_exp' * X(1:size(f_exp,1)) + sq_exp;
fc = f_exp(1:size(Ac,2)*tp);
ff = f_exp(size(Ac,2)*tp+1:size(Ac,2)*tp+size(Af,2));
fa = f_exp(size(Ac,2)*tp+size(Af,2)+1:size(Ac,2)*tp+size(Af,2)+size(inv_nsigma,1));
f_exp = f_exp(1:size(Ac,2)*tp+size(Af,2));
fa_effec = [];
for i = 0:tp-1
    single_fc = fc(1+size(Ac,2)*i:size(Ac,2)+size(Ac,2)*i);
    fa_effec = [fa_effec;find(single_fc~=0)];
end
fa_effec = [fa_effec;size(Ac,2)+find(ff~=0)];
ya = X(f_exp~=0) + X(size(f_exp,1)+fa_effec) - smpl(f_exp~=0);
as = X(size(f_exp,1)+find(fa~=0));
RSS = (ya./sigma(f_exp~=0))' * (ya./sigma(f_exp~=0));
Ps = length(as(-10^(-6)>=as | as>=10^(-6)));
resid = ya./sigma(f_exp~=0);

end