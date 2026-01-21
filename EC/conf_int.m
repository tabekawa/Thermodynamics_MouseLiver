function [CI] = conf_int(By,alpha)
By = squeeze(By);

if size(By,2) == 1
    By = reshape(By, 1, []);
end

for i = 1:size(By,1)
    ordinal = sort(By(i,:));
    minI = alpha/2*size(By,2);
    maxI = (1-alpha/2)*size(By,2);
    CI(i,:) = [ordinal(minI),ordinal(maxI)];
end

end