function demo_GLEAM_mouse()

load('MouseData.mat')

% creat the objective function for GLEAM
[H, F, inv_nsigma, sq, H_unm, F_unm, sq_unm] = create_obj_func(lnc_exp,sigma_c,G_exp_t,sigma_G_t,lb_met,ub_met,time_interval,14,unm_meas);

% estimate thermodynamically consistent metabolite concentrations and the standard Gibbs free energy of formation by GLEAM
[x,d,err] = GLEAM_opt(H,H_unm,F,F_unm,Ac,Af,CC,ratio,lb_met,ub_met,ub_met_ob,m_steady,sq,sq_unm,irr_pos,irr_neg,inv_nsigma,[lnc_exp;G_exp_t],[sigma_c;sigma_G_t],10^lambda,evector_cov,14,10,500);

% calculate the Gibbs free energy change of reaction
for i = 1:14
    drGs(:,i) = [Ac,Af*evector_cov] * [x(1+size(Ac,2)*(i-1):size(Ac,2)*i);x(1+size(Ac,2)*14:size(Ac,2)*14+size(Af,2))];
end
drGs(26,:) = drGs(26,:) + 2.578730581020227*log(10^(7.2-7.4));
drGs(29+24,:) = drGs(29+24,:) + 2.578730581020227*log(10^(7.2-7.4));

% calculate reduced chi-square using parametric bootstrap sampling
[By,redchi2] = cal_redchi2(Ac,Af,CC,ratio2,H,H_unm,F,F_unm,x,d,err,inv_nsigma,lb_metfree,ub_metfree,ub_metfree_obglucon,[lnc_exp;G_exp_t],[sigma_c;sigma_G_t],14,evector_cov,10^lambdas,m,sq,sq_unm,2000);

% calculate confidence intervals
CI = conf_int(By,0.05);
[ByDrG,Byratio] = By_DrG_ratio(Ac,Af,evector_cov,Cpair,By,14);
CI_DrG = conf_int(ByDrG,0.05);
CI_ratio = conf_int(Byratio,0.05);

save('MouseData_after_GLEAM.mat')

end