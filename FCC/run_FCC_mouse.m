function run_FCC_mouse()

load('MouseData_after_GLEAM.mat')
for i = 1:2000
    lnKeq = -( ( Af * evector_cov * By2000(size(Ac,2)*14+1:size(Ac,2)*14+size(Af,2),i) ) + ( Af * repmat(RT*log(10^(-3)),size(Af,2),1) ) ) ./ RT;
    [kcat_By(:,:,i,:), Km_By(:,:,i,:)] = gen_kinpararm(rxnnames,Km_ori,lnKeq,kcat_mean,kcat_sd,Km_mean,Km_sd,100);
end
kcat_By_norm = kcat_By ./ min(kcat_By(kcat_By~=0));
Km_By = Km_By ./ 1000;
Km_By(:,cofactors_l,:,:) = 1.*sign(Km_By(:,cofactors_l,:,:));

for i = 1:2000
    for j = 1:100
        [By_rev_glyc_fcc_l(:,i,j), By_sat_s_glyc_fcc_l(:,i,j), By_sat_p_glyc_fcc_l(:,i,j),By_sat_unm_glyc_fcc_l(:,i,j)] = cal_etas(Ac,Af,By2000(:,i),glyc_rxn_l(2:end),evector_cov,RT,Km_By(:,:,i,j),d,14,0);
        [By_rev_glucon_fcc_l(:,i,j), By_sat_s_glucon_fcc_l(:,i,j), By_sat_p_glucon_fcc_l(:,i,j),By_sat_unm_glucon_fcc_l(:,i,j)] = cal_etas(Ac,Af,By2000(:,i),glucon_rxn_l(1:end-1),evector_cov,RT,Km_By(:,:,i,j),d,14,6);
        [By_rev_TCA0_fcc_l(:,i,j), By_sat_s_TCA0_fcc_l(:,i,j), By_sat_p_TCA0_fcc_l(:,i,j),By_sat_unm_TCA0_fcc_l(:,i,j)] = cal_etas(Ac,Af,By2000(:,i),TCA_rxn_l,evector_cov,RT,Km_By(:,:,i,j),d,14,0);
        [By_rev_TCA16_fcc_l(:,i,j), By_sat_s_TCA16_fcc_l(:,i,j), By_sat_p_TCA16_fcc_l(:,i,j),By_sat_unm_TCA16_fcc_l(:,i,j)] = cal_etas(Ac,Af,By2000(:,i),TCA_rxn_l,evector_cov,RT,Km_By(:,:,i,j),d,14,6);

        [By_FCC_glyc_l(:,i,j), ~, By_elas_s_glyc_l(:,i,j), By_elas_p_glyc_l(:,i,j)] = calFluxCtrlCoef(By_rev_glyc_fcc_l(:,i,j),By_sat_s_glyc_fcc_l(:,i,j),By_sat_p_glyc_fcc_l(:,i,j),rxnnames_ec_glyc,1);
        [By_FCC_glucon_l(:,i,j), ~, By_elas_s_glucon_l(:,i,j), By_elas_p_glucon_l(:,i,j)] = calFluxCtrlCoef(By_rev_glucon_fcc_l(:,i,j),By_sat_s_glucon_fcc_l(:,i,j),By_sat_p_glucon_fcc_l(:,i,j),rxnnames_ec_glucon,2);
        [By_FCC_TCA0_l(:,i,j), ~, By_elas_s_TCA0_l(:,i,j), By_elas_p_TCA0_l(:,i,j)] = calFluxCtrlCoef(By_rev_TCA0_fcc_l(:,i,j),By_sat_s_TCA0_fcc_l(:,i,j),By_sat_p_TCA0_fcc_l(:,i,j),rxnnames_TCA,3);
        [By_FCC_TCA16_l(:,i,j), ~, By_elas_s_TCA16_l(:,i,j), By_elas_p_TCA16_l(:,i,j)] = calFluxCtrlCoef(By_rev_TCA16_fcc_l(:,i,j),By_sat_s_TCA16_fcc_l(:,i,j),By_sat_p_TCA16_fcc_l(:,i,j),rxnnames_TCA,3);
    end
end

end