function [BS,sd] = perturb_m(Ac,Af,opt,smpl,sigma,f,b,tp)

rng('shuffle')
fc = f(1:size(Ac,2)*tp);
ff = f(size(Ac,2)*tp+1:size(Ac,2)*tp+size(Af,2));
fa = f(size(Ac,2)*tp+size(Af,2)+1:size(Ac,2)*(tp+1)+size(Af,2)*2);
f = f(1:size(Ac,2)*tp+size(Af,2));
fa_effec = [];
for i = 0:tp-1
    single_fc = fc(1+size(Ac,2)*i:size(Ac,2)+size(Ac,2)*i);
    fa_effec = [fa_effec;find(single_fc~=0)];
end
fa_effec = [fa_effec;size(Ac,2)+find(ff~=0)];
ya = opt(f~=0) + opt(size(f,1)+fa_effec) - smpl(f~=0);
yas = ya./sigma(f~=0);
effn = length(f(f~=0));
RSmean = sum(yas)/effn;
RSS = sum((yas - RSmean).^2);
sd = sqrt(RSS/(effn-1));
BS = opt(1:size(f,1)) + sigma.*(randn(size(f,1),b).*sqrt(sd));

end