function [kcat_set, Km_set] = gen_kinpararm(rxn,Km_ori,lnKeq,kcat_mean,kcat_sd,Km_mean,Km_sd,numSamples)

for i = 1:length(rxn)
    if i == 1 || i == 26
        kcat_set(i,:,:) = zeros(1,2,numSamples);
        Km_set(i,:,:) = zeros(1,size(Km_ori,2),numSamples);
    elseif i > 29
        break
    else
        Km_r = Km_ori(i,:);
        kcat_mean_r = kcat_mean(i,:)';
        kcat_sd_r = kcat_sd(i,:)';
        Km_mean_r = Km_mean(i,:)';
        Km_sd_r = Km_sd(i,:)';
        [kcat, Km] = sample_kinparam_brenda(lnKeq(i),Km_r,kcat_mean_r,kcat_sd_r,Km_mean_r,Km_sd_r,numSamples);
        kcat_set(i,:,:) = kcat;
        Km_set(i,:,:) = Km;
    end
end

end