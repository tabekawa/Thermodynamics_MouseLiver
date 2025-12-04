function demo_GLEAM_mouse()

load('MouseData.mat')

[H, F, inv_nsigma, sq, H_unm, F_unm, sq_unm] = create_obj_func(lnc_exp,sigma_c,G_exp_t,sigma_G_t,lb_metfree,ub_metfree,time_interval,14,unm_mes);

[x,fval,ef,d,RSS,Ps,resid,output2,result2,Px2,Py2,err] = GLEAM_opt(H,H_unm,F,F_unm,Ac,Af,CC,ratio2,lb_metfree,ub_metfree,ub_metfree_obglucon,m,sq,sq_unm,irr_pos,irr_neg,inv_nsigma,[lnc_exp;G_exp_t],[sigma_c;sigma_G_t],10^lambdas,evector_cov,14,10,500);

for i = 1:14
    drGs(:,i) = [Ac,Af*evector_cov] * [x(1+size(Ac,2)*(i-1):size(Ac,2)*i);x(1+size(Ac,2)*14:size(Ac,2)*14+size(Af,2))];
end
drGs(26,:) = drGs(26,:) + 2.578730581020227*log(10^(7.2-7.4));
drGs(29+24,:) = drGs(29+24,:) + 2.578730581020227*log(10^(7.2-7.4));

[By2000,results1,results2,redchi2,fval] = cal_redchi2(Ac,Af,CC,ratio2,H,H_unm,F,F(1:size(Ac,2)*14+size(Af,2)),F_unm,x,d,err,inv_nsigma,lb_metfree,ub_metfree,ub_metfree_obglucon,[lnc_exp;G_exp_t],[sigma_c;sigma_G_t],14,evector_cov,10^lambdas,m,sq,sq_unm,2000);

CI = conf_int(By2000,0.05);
[ByGr,Byratio] = ByGrratio(Ac,Af,evector_cov,Cratio,By2000,14);
CI_DrG = conf_int(ByGr,0.05);
CI_ratio = conf_int(Byratio,0.05);

save('MouseData_after_GLEAM.mat')

end