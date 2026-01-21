function [Z, X] = opt_measured(H,f,Ac,Af,C,ratio,lb,ub,ub_ob,m,sq,inv_nsigma,lambda,T,tp,Acs,Cs,Ss,rG_lim,d)

model.Q = sparse([H,zeros(size(H,1),size(Ac,1)*tp+size(inv_nsigma,1));zeros(size(Ac,1)*tp+size(inv_nsigma,1),size(H,2)+size(Ac,1)*tp+size(inv_nsigma,1))])./2;
model.obj = [f;zeros(size(Ac,1)*tp,1);(inv_nsigma.*lambda)];
model.A = sparse([[Acs,repmat(Af*T,tp,1),zeros(size(Ac,1)*tp,size(inv_nsigma,1)*2+size(Ac,1)*tp)].*d'; ...
    Cs,zeros(size(C,1)*tp,size(H,1)-size(C,2)*tp+size(Ac,1)*tp+size(inv_nsigma,1));...
    zeros(size(inv_nsigma,1),size(H,1)-size(inv_nsigma,1)),[diag(ones(size(Ac,2),1)),zeros(size(Ac,2),size(inv_nsigma,1)-size(Ac,2));zeros(size(inv_nsigma,1)-size(Ac,2),size(Ac,2)),T].*(-1),zeros(size(inv_nsigma,1),size(Ac,1)*tp),diag(repmat(-1,size(inv_nsigma,1),1)); ...
    zeros(size(inv_nsigma,1),size(H,1)-size(inv_nsigma,1)),[diag(ones(size(Ac,2),1)),zeros(size(Ac,2),size(inv_nsigma,1)-size(Ac,2));zeros(size(inv_nsigma,1)-size(Ac,2),size(Ac,2)),T],zeros(size(inv_nsigma,1),size(Ac,1)*tp),diag(repmat(-1,size(inv_nsigma,1),1)); ...
    zeros(size(m,1)*tp,size(H,2)),(Ss.*d')',zeros(size(m,1)*tp,size(inv_nsigma,1)); ...
    [zeros(size(m,1)*tp,size(H,2)),(Ss.*d')',zeros(size(m,1)*tp,size(inv_nsigma,1))].*(-1)]);
model.rhs = [rG_lim;repmat(ratio,tp,1);zeros(size(inv_nsigma,1)*2+size(m,1)*tp*2,1)];
model.lb = [repmat(lb,tp,1);repmat(-10^4,size(Af,2),1);repmat(-10^4,size(inv_nsigma,1),1);repmat(10^2,size(Ac,1)*tp,1);zeros(size(inv_nsigma,1),1)];
model.ub = [repmat(ub,tp/2,1);repmat(ub_ob,tp/2,1);repmat(10^4,size(Af,2),1);repmat(10^4,size(inv_nsigma,1),1);repmat(10^3,size(Ac,1)*tp,1);repmat(10^4,size(inv_nsigma,1),1)];
model.objcon = sq;

param.BarHomogeneous = 1;
param.OptimalityTol = 10^(-9);
param.FeasibilityTol = 10^(-9);
param.NumericFocus = 3;

RESULT = gurobi(model,param);
if strcmp(RESULT.status, 'OPTIMAL') || strcmp(RESULT.status, 'SUBOPTIMAL')
    Z = RESULT.objval;
    X = RESULT.x;
elseif strcmp(RESULT.status, 'INFEASIBLE') || strcmp(RESULT.status, 'INF_OR_UNBD')
    Z = 10;
    X = zeros(length(f)+size(Ac,1)*tp+length(inv_nsigma), 1);
else
    Z = 10;
    X = zeros(length(f)+size(Ac,1)*tp+length(inv_nsigma), 1);
    disp(RESULT.status)
end

end