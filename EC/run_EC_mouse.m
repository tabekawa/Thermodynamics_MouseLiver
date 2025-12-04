function run_EC_mouse()

load('MouseData_after_GLEAM.mat');
Cs = repmat("Cratio",1,14);
Cs = join(Cs,",");
Cs = ["blkdiag(",Cs,")"];
Cs = join(Cs,"");
Cs = eval(Cs);
x_ratio = Cs * x(1:size(Ac,2)*14);
clearvars Cs

x_ec_glyc_l = min_enzyme_cost(Ac,Af,Cratio,x_ratio,lb_metfree,ub_metfree,x,glyc_rxn_l(2:end),evector_cov,14,RT, ...
    [Molecular_Weight_norm(glyc_rxn_l(2:6)),2.*Molecular_Weight_norm(glyc_rxn_l(7:end))],Km_ori,d,zeros(size(Ac,2),1),100,fixmet_glygluc_l,0);
x_ec_glucon_l = min_enzyme_cost(Ac,Af,Cratio,x_ratio,lb_metfree,ub_metfree,x,glucon_rxn_l(1:end-1),evector_cov,14,RT, ...
    [2.*Molecular_Weight_norm(glucon_rxn_l(1:5)),Molecular_Weight_norm(glucon_rxn_l(6:end-1))],Km_ori,d,zeros(size(Ac,2),1),100,fixmet_glygluc_l,6);
x_ec_TCA0_l = min_enzyme_cost(Ac,Af,Cratio,x_ratio,lb_metfree,ub_metfree,x,TCA_rxn_l,evector_cov,14,RT, ...
    Molecular_Weight_norm(TCA_rxn_l),Km_ori,d,zeros(size(Ac,2),1),100,fixmet_TCA_l,0);
x_ec_TCA16_l = min_enzyme_cost(Ac,Af,Cratio,x_ratio,lb_metfree,ub_metfree,x,TCA_rxn_l,evector_cov,14,RT, ...
    Molecular_Weight_norm(TCA_rxn_l),Km_ori,d,zeros(size(Ac,2),1),100,fixmet_TCA_l,6);

[EC_tot_glyc_l(:,2), E_glyc_l(:,2), rev_glyc_ec_l(:,2), sat_s_glyc_ec_l(:,2), sat_p_glyc_ec_l(:,2),~,sat_unm_glyc_ec_l(:,2)] = enzyme_cost(Ac,Af,x_ec_glyc_l,x,glyc_rxn_l(2:end),evector_cov,RT,[Molecular_Weight_norm(glyc_rxn_l(2:6)),2.*Molecular_Weight_norm(glyc_rxn_l(7:end))],Km_ori,Km_micofactors,[],d,14,0);
[EC_tot_glucon_l(:,2), E_glucon_l(:,2), rev_glucon_ec_l(:,2), sat_s_glucon_ec_l(:,2), sat_p_glucon_ec_l(:,2),~,sat_unm_glucon_ec_l(:,2)] = enzyme_cost(Ac,Af,x_ec_glucon_l,x,glucon_rxn_l(1:end-1),evector_cov,RT,[2.*Molecular_Weight_norm(glucon_rxn_l(1:5)),Molecular_Weight_norm(glucon_rxn_l(6:end-1))],Km_ori,Km_micofactors,[],d,14,6);
[EC_tot_TCA0_l(:,2), E_TCA0_l(:,2), rev_TCA0_ec_l(:,2), sat_s_TCA0_ec_l(:,2), sat_p_TCA0_ec_l(:,2),~,sat_unm_TCA0_ec_l(:,2)] = enzyme_cost(Ac,Af,x_ec_TCA0_l,x,TCA_rxn_l,evector_cov,RT,Molecular_Weight_norm(TCA_rxn_l),Km_ori,Km_micofactors,[],d,14,0);
[EC_tot_TCA16_l(:,2), E_TCA16_l(:,2), rev_TCA16_ec_l(:,2), sat_s_TCA16_ec_l(:,2), sat_p_TCA16_ec_l(:,2),~,sat_unm_TCA16_ec_l(:,2)] = enzyme_cost(Ac,Af,x_ec_TCA16_l,x,TCA_rxn_l,evector_cov,RT,Molecular_Weight_norm(TCA_rxn_l),Km_ori,Km_micofactors,[],d,14,6);

[EC_tot_glyc_l(:,1), E_glyc_l(:,1), rev_glyc_ec_l(:,1), sat_s_glyc_ec_l(:,1), sat_p_glyc_ec_l(:,1),~,sat_unm_glyc_ec_l(:,1)] = enzyme_cost_invivo(Ac,Af,x,glyc_rxn_l(2:end),evector_cov,RT,[Molecular_Weight_norm(glyc_rxn_l(2:6)),2.*Molecular_Weight_norm(glyc_rxn_l(7:end))],Km_ori,Km_micofactors,[],d,14,0);
[EC_tot_glucon_l(:,1), E_glucon_l(:,1), rev_glucon_ec_l(:,1), sat_s_glucon_ec_l(:,1), sat_p_glucon_ec_l(:,1),~,sat_unm_glucon_ec_l(:,1)] = enzyme_cost_invivo(Ac,Af,x,glucon_rxn_l(1:end-1),evector_cov,RT,[2.*Molecular_Weight_norm(glucon_rxn_l(1:5)),Molecular_Weight_norm(glucon_rxn_l(6:end-1))],Km_ori,Km_micofactors,[],d,14,6);
[EC_tot_TCA0_l(:,1), E_TCA0_l(:,1), rev_TCA0_ec_l(:,1), sat_s_TCA0_ec_l(:,1), sat_p_TCA0_ec_l(:,1),~,sat_unm_TCA0_ec_l(:,1)] = enzyme_cost_invivo(Ac,Af,x,TCA_rxn_l,evector_cov,RT,Molecular_Weight_norm(TCA_rxn_l),Km_ori,Km_micofactors,[],d,14,0);
[EC_tot_TCA16_l(:,1), E_TCA16_l(:,1), rev_TCA16_ec_l(:,1), sat_s_TCA16_ec_l(:,1), sat_p_TCA16_ec_l(:,1),~,sat_unm_TCA16_ec_l(:,1)] = enzyme_cost_invivo(Ac,Af,x,TCA_rxn_l,evector_cov,RT,Molecular_Weight_norm(TCA_rxn_l),Km_ori,Km_micofactors,[],d,14,6);

RMSE_ECMmets_glycglucon = sqrt(mean(( log10(exp(x([mets_glyc_meas_l+size(Ac,2)*0;mets_glucon_meas_l+size(Ac,2)*6]))) - log10(exp([x_ec_glyc_l(mets_glyc_meas_l);x_ec_glucon_l(mets_glucon_meas_l)])) ).^2));
RMSE_ECMmets_TCA = sqrt(mean(( log10(exp(x([mets_TCA_meas_l+size(Ac,2)*0;mets_TCA_meas_l+size(Ac,2)*6]))) - log10(exp([x_ec_TCA0_l(mets_TCA_meas_l);x_ec_TCA16_l(mets_TCA_meas_l)])) ).^2));

[meg_glyc(2), x_meg_glyc(:,2)] = cal_meg(x_ec_glyc_l,x,Ac,Af,evector_cov,CC,ratio2,lb_metfree,ub_metfree,glyc_rxn_l(2:end),rxn_liver_rev,mets_glyc_l(2:end),d,0,RT);
[meg_glucon(2), x_meg_glucon(:,2)] = cal_meg(x_ec_glucon_l,x,Ac,Af,evector_cov,CC,ratio2,lb_metfree,ub_metfree,glucon_rxn_l(1:end-1),rxn_liver_rev,mets_glucon_l(2:end),d,6,RT);
[meg_TCA0(2), x_meg_TCA0(:,2)] = cal_meg(x_ec_TCA0_l,x,Ac,Af,evector_cov,CC,ratio2,lb_metfree,ub_metfree,TCA_rxn_l,TCA_rxn_l,mets_TCA_l,d,0,RT);
[meg_TCA16(2), x_meg_TCA16(:,2)] = cal_meg(x_ec_TCA16_l,x,Ac,Af,evector_cov,CC,ratio2,lb_metfree,ub_metfree,TCA_rxn_l,TCA_rxn_l,mets_TCA_l,d,6,RT);
[meg_glyc(1), x_meg_glyc(:,1)] = cal_meg_invivo(x,Ac,Af,evector_cov,CC,ratio2,lb_metfree,ub_metfree,glyc_rxn_l(2:end),rxn_liver_rev,mets_glyc_l(2:end),d,0,RT);
[meg_glucon(1), x_meg_glucon(:,1)] = cal_meg_invivo(x,Ac,Af,evector_cov,CC,ratio2,lb_metfree,ub_metfree,glucon_rxn_l(1:end-1),rxn_liver_rev,mets_glucon_l(2:end),d,6,RT);
[meg_TCA0(1), x_meg_TCA0(:,1)] = cal_meg_invivo(x,Ac,Af,evector_cov,CC,ratio2,lb_metfree,ub_metfree,TCA_rxn_l,TCA_rxn_l,mets_TCA_l,d,0,RT);
[meg_TCA16(1), x_meg_TCA16(:,1)] = cal_meg_invivo(x,Ac,Af,evector_cov,CC,ratio2,lb_metfree,ub_metfree,TCA_rxn_l,TCA_rxn_l,mets_TCA_l,d,6,RT);

meg_glgl(2) = exp(sqrt( (log(meg_glyc(2))^2*size(mets_glyc_l(2:end),1) + log(meg_glucon(2))^2*size(mets_glucon_l(2:end),1)) / (size(mets_glyc_l(2:end),1)+size(mets_glucon_l(2:end),1)) ));
meg_glgl(1) = exp(sqrt( (log(meg_glyc(1))^2*size(mets_glyc_l(2:end),1) + log(meg_glucon(1))^2*size(mets_glucon_l(2:end),1)) / (size(mets_glyc_l(2:end),1)+size(mets_glucon_l(2:end),1)) ));
meg_TCA(2) = exp(sqrt( (log(meg_TCA0(2))^2*size(mets_TCA_l,1) + log(meg_TCA16(2))^2*size(mets_TCA_l,1)) / (size(mets_TCA_l,1)*2) ));
meg_TCA(1) = exp(sqrt( (log(meg_TCA0(1))^2*size(mets_TCA_l,1) + log(meg_TCA16(1))^2*size(mets_TCA_l,1)) / (size(mets_TCA_l,1)*2) ));


for i = 1:2000
    By_x_ec_glyc_l(:,i) = min_enzyme_cost(Ac,Af,Cratio,Byratio(:,i),lb_metfree,ub_metfree,By2000(:,i),glyc_rxn_l(2:end),evector_cov,14,RT, ...
        [Molecular_Weight_norm(glyc_rxn_l(2:6)),2.*Molecular_Weight_norm(glyc_rxn_l(7:end))],Km_ori,d,zeros(size(Ac,2),1),100,fixmet_glygluc_l,0);
    By_x_ec_glucon_l(:,i) = min_enzyme_cost(Ac,Af,Cratio,Byratio(:,i),lb_metfree,ub_metfree,By2000(:,i),glucon_rxn_l(1:end-1),evector_cov,14,RT, ...
        [2.*Molecular_Weight_norm(glucon_rxn_l(1:5)),Molecular_Weight_norm(glucon_rxn_l(6:end-1))],Km_ori,d,zeros(size(Ac,2),1),100,fixmet_glygluc_l,6);
    By_x_ec_TCA0_l(:,i) = min_enzyme_cost(Ac,Af,Cratio,Byratio(:,i),lb_metfree,ub_metfree,By2000(:,i),TCA_rxn_l,evector_cov,14,RT, ...
        Molecular_Weight_norm(TCA_rxn_l),Km_ori,d,zeros(size(Ac,2),1),100,fixmet_TCA_l,0);
    By_x_ec_TCA16_l(:,i) = min_enzyme_cost(Ac,Af,Cratio,Byratio(:,i),lb_metfree,ub_metfree,By2000(:,i),TCA_rxn_l,evector_cov,14,RT, ...
        Molecular_Weight_norm(TCA_rxn_l),Km_ori,d,zeros(size(Ac,2),1),100,fixmet_TCA_l,6);

    [By_EC_tot_glyc_l(:,2,i), By_E_glyc_l(:,2,i), By_rev_glyc_ec_l(:,2,i), By_sat_s_glyc_ec_l(:,2,i), By_sat_p_glyc_ec_l(:,2,i),~,By_sat_unm_glyc_ec_l(:,2,i)] = enzyme_cost(Ac,Af,By_x_ec_glyc_l(:,i),By2000(:,i),glyc_rxn_l(2:end),evector_cov,RT,[Molecular_Weight_norm(glyc_rxn_l(2:6)),2.*Molecular_Weight_norm(glyc_rxn_l(7:end))],Km_ori,Km_micofactors,kcat_By_norm(:,:,i,1),d,14,0);
    [By_EC_tot_glucon_l(:,2,i), By_E_glucon_l(:,2,i), By_rev_glucon_ec_l(:,2,i), By_sat_s_glucon_ec_l(:,2,i), By_sat_p_glucon_ec_l(:,2,i),~,By_sat_unm_glucon_ec_l(:,2,i)] = enzyme_cost(Ac,Af,By_x_ec_glucon_l(:,i),By2000(:,i),glucon_rxn_l(1:end-1),evector_cov,RT,[2.*Molecular_Weight_norm(glucon_rxn_l(1:5)),Molecular_Weight_norm(glucon_rxn_l(6:end-1))],Km_ori,Km_micofactors,kcat_By_norm(:,:,i,1),d,14,6);
    [By_EC_tot_TCA0_l(:,2,i), By_E_TCA0_l(:,2,i), By_rev_TCA0_ec_l(:,2,i), By_sat_s_TCA0_ec_l(:,2,i), By_sat_p_TCA0_ec_l(:,2,i),~,By_sat_unm_TCA0_ec_l(:,2,i)] = enzyme_cost(Ac,Af,By_x_ec_TCA0_l(:,i),By2000(:,i),TCA_rxn_l,evector_cov,RT,Molecular_Weight_norm(TCA_rxn_l),Km_ori,Km_micofactors,kcat_By_norm(:,:,i,1),d,14,0);
    [By_EC_tot_TCA16_l(:,2,i), By_E_TCA16_l(:,2,i), By_rev_TCA16_ec_l(:,2,i), By_sat_s_TCA16_ec_l(:,2,i), By_sat_p_TCA16_ec_l(:,2,i),~,By_sat_unm_TCA16_ec_l(:,2,i)] = enzyme_cost(Ac,Af,By_x_ec_TCA16_l(:,i),By2000(:,i),TCA_rxn_l,evector_cov,RT,Molecular_Weight_norm(TCA_rxn_l),Km_ori,Km_micofactors,kcat_By_norm(:,:,i,1),d,14,6);

    [By_EC_tot_glyc_l(:,1,i), By_E_glyc_l(:,1,i), By_rev_glyc_ec_l(:,1,i), By_sat_s_glyc_ec_l(:,1,i), By_sat_p_glyc_ec_l(:,1,i),~,By_sat_unm_glyc_ec_l(:,1,i)] = enzyme_cost_invivo(Ac,Af,By2000(:,i),glyc_rxn_l(2:end),evector_cov,RT,[Molecular_Weight_norm(glyc_rxn_l(2:6)),2.*Molecular_Weight_norm(glyc_rxn_l(7:end))],Km_ori,Km_micofactors,kcat_By_norm(:,:,i,1),d,14,0);
    [By_EC_tot_glucon_l(:,1,i), By_E_glucon_l(:,1,i), By_rev_glucon_ec_l(:,1,i), By_sat_s_glucon_ec_l(:,1,i), By_sat_p_glucon_ec_l(:,1,i),~,By_sat_unm_glucon_ec_l(:,1,i)] = enzyme_cost_invivo(Ac,Af,By2000(:,i),glucon_rxn_l(1:end-1),evector_cov,RT,[2.*Molecular_Weight_norm(glucon_rxn_l(1:5)),Molecular_Weight_norm(glucon_rxn_l(6:end-1))],Km_ori,Km_micofactors,kcat_By_norm(:,:,i,1),d,14,6);
    [By_EC_tot_TCA0_l(:,1,i), By_E_TCA0_l(:,1,i), By_rev_TCA0_ec_l(:,1,i), By_sat_s_TCA0_ec_l(:,1,i), By_sat_p_TCA0_ec_l(:,1,i),~,By_sat_unm_TCA0_ec_l(:,1,i)] = enzyme_cost_invivo(Ac,Af,By2000(:,i),TCA_rxn_l,evector_cov,RT,Molecular_Weight_norm(TCA_rxn_l),Km_ori,Km_micofactors,kcat_By_norm(:,:,i,1),d,14,0);
    [By_EC_tot_TCA16_l(:,1,i), By_E_TCA16_l(:,1,i), By_rev_TCA16_ec_l(:,1,i), By_sat_s_TCA16_ec_l(:,1,i), By_sat_p_TCA16_ec_l(:,1,i),~,By_sat_unm_TCA16_ec_l(:,1,i)] = enzyme_cost_invivo(Ac,Af,By2000(:,i),TCA_rxn_l,evector_cov,RT,Molecular_Weight_norm(TCA_rxn_l),Km_ori,Km_micofactors,kcat_By_norm(:,:,i,1),d,14,6);

    [~, By_x_meg_glyc(:,2,i)] = cal_meg(By_x_ec_glyc_l(:,i),By2000(:,i),Ac,Af,evector_cov,CC,ratio2,lb_metfree,ub_metfree,intersect(glyc_rxn_l(2:end),rxn_liver_rev),mets_glyc_l(2:end),d,0,RT);
    [~, By_x_meg_glucon(:,2,i)] = cal_meg(By_x_ec_glucon_l(:,i),By2000(:,i),Ac,Af,evector_cov,CC,ratio2,lb_metfree,ub_metfree,intersect(glucon_rxn_l(1:end-1),rxn_liver_rev),mets_glucon_l(2:end),d,6,RT);
    [~, By_x_meg_TCA0(:,2,i)] = cal_meg(By_x_ec_TCA0_l(:,i),By2000(:,i),Ac,Af,evector_cov,CC,ratio2,lb_metfree,ub_metfree,TCA_rxn_l,mets_TCA_l,d,0,RT);
    [~, By_x_meg_TCA16(:,2,i)] = cal_meg(By_x_ec_TCA16_l(:,i),By2000(:,i),Ac,Af,evector_cov,CC,ratio2,lb_metfree,ub_metfree,TCA_rxn_l,mets_TCA_l,d,6,RT);
    [~, By_x_meg_glyc(:,1,i)] = cal_meg_invivo(By2000(:,i),Ac,Af,evector_cov,CC,ratio2,lb_metfree,ub_metfree,intersect(glyc_rxn_l(2:end),rxn_liver_rev),mets_glyc_l(2:end),d,0,RT);
    [~, By_x_meg_glucon(:,1,i)] = cal_meg_invivo(By2000(:,i),Ac,Af,evector_cov,CC,ratio2,lb_metfree,ub_metfree,intersect(glucon_rxn_l(1:end-1),rxn_liver_rev),mets_glucon_l(2:end),d,6,RT);
    [~, By_x_meg_TCA0(:,1,i)] = cal_meg_invivo(By2000(:,i),Ac,Af,evector_cov,CC,ratio2,lb_metfree,ub_metfree,TCA_rxn_l,mets_TCA_l,d,0,RT);
    [~, By_x_meg_TCA16(:,1,i)] = cal_meg_invivo(By2000(:,i),Ac,Af,evector_cov,CC,ratio2,lb_metfree,ub_metfree,TCA_rxn_l,mets_TCA_l,d,6,RT);
end

CI_EC_tot_glyc_l_ECM = conf_int(By_EC_tot_glyc_l(:,2,:));
CI_EC_tot_glucon_l_ECM = conf_int(By_EC_tot_glucon_l(:,2,:));
CI_EC_tot_TCA0_l_ECM = conf_int(By_EC_tot_TCA0_l(:,2,:));
CI_EC_tot_TCA16_l_ECM = conf_int(By_EC_tot_TCA16_l(:,2,:));

CI_EC_tot_glyc_l_liver = conf_int(By_EC_tot_glyc_l(:,1,:));
CI_EC_tot_glucon_l_liver = conf_int(By_EC_tot_glucon_l(:,1,:));
CI_EC_tot_TCA0_l_liver = conf_int(By_EC_tot_TCA0_l(:,1,:));
CI_EC_tot_TCA16_l_liver = conf_int(By_EC_tot_TCA16_l(:,1,:));

CI_E_glyc_l_ECM = conf_int(By_E_glyc_l(:,2,:));
CI_E_glucon_l_ECM = conf_int(By_E_glucon_l(:,2,:));
CI_E_TCA0_l_ECM = conf_int(By_E_TCA0_l(:,2,:));
CI_E_TCA16_l_ECM = conf_int(By_E_TCA16_l(:,2,:));

CI_E_glyc_l_liver = conf_int(By_E_glyc_l(:,1,:));
CI_E_glucon_l_liver = conf_int(By_E_glucon_l(:,1,:));
CI_E_TCA0_l_liver = conf_int(By_E_TCA0_l(:,1,:));
CI_E_TCA16_l_liver = conf_int(By_E_TCA16_l(:,1,:));

CI_rev_glyc_ec_l_ECM = conf_int(By_rev_glyc_ec_l(:,2,:));
CI_rev_glucon_ec_l_ECM = conf_int(By_rev_glucon_ec_l(:,2,:));
CI_rev_TCA0_ec_l_ECM = conf_int(By_rev_TCA0_ec_l(:,2,:));
CI_rev_TCA16_ec_l_ECM = conf_int(By_rev_TCA16_ec_l(:,2,:));

CI_rev_glyc_ec_l_liver = conf_int(By_rev_glyc_ec_l(:,1,:));
CI_rev_glucon_ec_l_liver = conf_int(By_rev_glucon_ec_l(:,1,:));
CI_rev_TCA0_ec_l_liver = conf_int(By_rev_TCA0_ec_l(:,1,:));
CI_rev_TCA16_ec_l_liver = conf_int(By_rev_TCA16_ec_l(:,1,:));

CI_sat_s_glyc_ec_l_ECM = conf_int(By_sat_s_glyc_ec_l(:,2,:));
CI_sat_s_glucon_ec_l_ECM = conf_int(By_sat_s_glucon_ec_l(:,2,:));
CI_sat_s_TCA0_ec_l_ECM = conf_int(By_sat_s_TCA0_ec_l(:,2,:));
CI_sat_s_TCA16_ec_l_ECM = conf_int(By_sat_s_TCA16_ec_l(:,2,:));

CI_sat_s_glyc_ec_l_liver = conf_int(By_sat_s_glyc_ec_l(:,1,:));
CI_sat_s_glucon_ec_l_liver = conf_int(By_sat_s_glucon_ec_l(:,1,:));
CI_sat_s_TCA0_ec_l_liver = conf_int(By_sat_s_TCA0_ec_l(:,1,:));
CI_sat_s_TCA16_ec_l_liver = conf_int(By_sat_s_TCA16_ec_l(:,1,:));

CI_x_meg_glyc_ECM = conf_int(By_x_meg_glyc(:,2,:));
CI_x_meg_glucon_ECM = conf_int(By_x_meg_glucon(:,2,:));
CI_x_meg_TCA0_ECM = conf_int(By_x_meg_TCA0(:,2,:));
CI_x_meg_TCA16_ECM = conf_int(By_x_meg_TCA16(:,2,:));

CI_x_meg_glyc_liver = conf_int(By_x_meg_glyc(:,1,:));
CI_x_meg_glucon_liver = conf_int(By_x_meg_glucon(:,1,:));
CI_x_meg_TCA0_liver = conf_int(By_x_meg_TCA0(:,1,:));
CI_x_meg_TCA16_liver = conf_int(By_x_meg_TCA16(:,1,:));

end