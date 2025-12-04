function [meg, c] = cal_meg_invivo(x,Ac,Af,T,CC,ratio,lb,ub,rxn,rev,mets,d,h,RT)

revrxn = intersect(rev,rxn);
tp=14;
rG_lim_rev = zeros(size(Ac,1),1);
rG_lim_rev = rG_lim_rev - Af*T * x(size(Ac,2)*tp+1:size(Ac,2)*tp+size(Af,2));
rG_lim_rev = rG_lim_rev(revrxn);

d = d(1+h*size(Ac,1):(h+1)*size(Ac,1));
irrev_idx = setdiff(rxn, revrxn); % aに含まれないインデックスを取得
rG_lim_irrev = repmat(-0.01,size(Ac,1),1);
rG_lim_irrev(26) = rG_lim_irrev(26) - RT*log(10^(d(26)*(7.2-7.4)));
rG_lim_irrev(29+24) = rG_lim_irrev(29+24) - RT*log(10^(d(29+24)*(7.2-7.4)));
rG_lim_irrev = rG_lim_irrev - Af*T.*d * x(size(Ac,2)*tp+1:size(Ac,2)*tp+size(Af,2));
rG_lim_irrev = rG_lim_irrev(irrev_idx);

model.Q = sparse(diag(ones(size(mets,1),1)));
model.obj = -2 .* x(mets+h*size(Ac,2));
model.A = sparse([Ac(revrxn,mets);-1.*Ac(revrxn,mets);Ac(irrev_idx,mets).*d(irrev_idx)]);
model.rhs = [rG_lim_rev;-1.*rG_lim_rev;rG_lim_irrev];
model.lb = lb(mets);
model.ub = ub(mets);
model.objcon = sum(x(mets+h*size(Ac,2)).^2);

param.BarHomogeneous = 1;
param.OptimalityTol = 10^(-9);
param.FeasibilityTol = 10^(-9);
param.NumericFocus = 3;

RESULT = gurobi(model,param);
if strcmp(RESULT.status, 'OPTIMAL') || strcmp(RESULT.status, 'SUBOPTIMAL')
    meg = exp(sqrt(RESULT.objval / size(mets,1)));
    c = x(1+h*size(Ac,2):(1+h)*size(Ac,2));
    c(mets) = RESULT.x;
    c = exp(c);
elseif strcmp(RESULT.status, 'INFEASIBLE') || strcmp(RESULT.status, 'INF_OR_UNBD')
    meg = 10^8;
    c = repmat(10^8,size(Ac,2),1);
end

end