function [ByGr,Byratio] = By_DrG_ratio(Ac,Af,T,C,By,tp)

Acs = repmat("Ac",1,tp);
Acs = join(Acs,",");
Acs = ["blkdiag(",Acs,")"];
Acs = join(Acs,"");
Acs = eval(Acs);
Cs = repmat("C",1,tp);
Cs = join(Cs,",");
Cs = ["blkdiag(",Cs,")"];
Cs = join(Cs,"");
Cs = eval(Cs);

ByGr = [Acs,repmat(Af*T,tp,1)] * By(1:size(Ac,2)*tp+size(Af,2),:);
for r = 0:tp-1
    ByGr(26+size(Ac,1)*r,:) = ByGr(26+size(Ac,1)*r,:) + 2.578730581020227*log(10^(7.2-7.4));
    ByGr(29+24+size(Ac,1)*r,:) = ByGr(29+24+size(Ac,1)*r,:) + 2.578730581020227*log(10^(7.2-7.4));
end
Byratio = Cs * By(1:size(Ac,2)*tp,:);

end