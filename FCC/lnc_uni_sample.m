function lnc_us = lnc_uni_sample(x,Ac,Af,CC,ratio,lb,ub,d,T,RT,tp,fixmet,h)

d = d(1+h*size(Ac,1):(h+1)*size(Ac,1));
lb(lb>log(10^-7)) = log(10^-7);
lb([16,56]) = 0;
ub(ub<log(10^-2)) = log(10^-2);
lb(fixmet) = x(fixmet+size(Ac,2)*h);
ub(fixmet) = x(fixmet+size(Ac,2)*h);
rG_lim = repmat(-0.01,size(Ac,1),1);
rG_lim(26) = -0.01 - RT*log(10^(d(26)*(7.2-7.4)));
rG_lim(29+24) = -0.01 - RT*log(10^(d(29+24)*(7.2-7.4)));
rG_lim = rG_lim - Af*T.*d * x(size(Ac,2)*tp+1:size(Ac,2)*tp+size(Af,2));
P.A = [Ac.*d;diag(ones(size(Ac,2),1));diag(repmat(-1,size(Ac,2),1))];
P.b = [rG_lim;ub;-1.*lb];
lnc_us = polySampler(P,size(Ac,2)^2,20);

end