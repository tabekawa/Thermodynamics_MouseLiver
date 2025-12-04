function [Z_unm, X_unm] = opt_unmeasured(H_unm,f_unm,Ac,Af,C,ratio,lb,ub,m,sq_unm,inv_nsigma_unm,lambda,T,tp,Acs,Cs,Ss,rG_lim,d)

model.Q = sparse([H_unm,zeros(size(H_unm,1),size(Ac,1)*tp+size(inv_nsigma_unm,1));zeros(size(Ac,1)*tp+size(inv_nsigma_unm,1),size(H_unm,2)+size(Ac,1)*tp+size(inv_nsigma_unm,1))])./2;
model.obj = [f_unm;zeros(size(Ac,1)*tp,1);(inv_nsigma_unm.*lambda)];
model.A = sparse([[Acs,repmat(Af*T,tp,1),zeros(size(Ac,1)*tp,size(inv_nsigma_unm,1)*2+size(Ac,1)*tp)].*d'; ...
    Cs,zeros(size(C,1)*tp,size(H_unm,1)-size(C,2)*tp+size(Ac,1)*tp+size(inv_nsigma_unm,1));...
    zeros(size(inv_nsigma_unm,1),size(H_unm,1)-size(inv_nsigma_unm,1)),[diag(ones(size(Ac,2),1)),zeros(size(Ac,2),size(inv_nsigma_unm,1)-size(Ac,2));zeros(size(inv_nsigma_unm,1)-size(Ac,2),size(Ac,2)),T].*(-1),zeros(size(inv_nsigma_unm,1),size(Ac,1)*tp),diag(repmat(-1,size(inv_nsigma_unm,1),1)); ...
    zeros(size(inv_nsigma_unm,1),size(H_unm,1)-size(inv_nsigma_unm,1)),[diag(ones(size(Ac,2),1)),zeros(size(Ac,2),size(inv_nsigma_unm,1)-size(Ac,2));zeros(size(inv_nsigma_unm,1)-size(Ac,2),size(Ac,2)),T],zeros(size(inv_nsigma_unm,1),size(Ac,1)*tp),diag(repmat(-1,size(inv_nsigma_unm,1),1)); ...
    zeros(size(m,1)*tp,size(H_unm,2)),(Ss.*d')',zeros(size(m,1)*tp,size(inv_nsigma_unm,1)); ...
    [zeros(size(m,1)*tp,size(H_unm,2)),(Ss.*d')',zeros(size(m,1)*tp,size(inv_nsigma_unm,1))].*(-1)]);
model.rhs = [rG_lim;repmat(ratio,tp,1);zeros(size(inv_nsigma_unm,1)*2+size(m,1)*tp*2,1)];
model.lb = [lb;repmat(10^2,size(Ac,1)*tp,1);zeros(size(inv_nsigma_unm,1),1)];
model.ub = [ub;repmat(10^3,size(Ac,1)*tp,1);repmat(10^4,size(inv_nsigma_unm,1),1)];
model.objcon = sq_unm;

param.BarHomogeneous = 1;
param.OptimalityTol = 10^(-9);
param.FeasibilityTol = 10^(-9);
param.NumericFocus = 3;

RESULT = gurobi(model,param);
if strcmp(RESULT.status, 'OPTIMAL') || strcmp(RESULT.status, 'SUBOPTIMAL')
    Z_unm = RESULT.objval;
    if Z_unm > 10^6
        Z_unm = 10^6;
    end
    X_unm = RESULT.x;
elseif strcmp(RESULT.status, 'INFEASIBLE') || strcmp(RESULT.status, 'INF_OR_UNBD')
    Z_unm = 10^11;
    X_unm = zeros(length(f_unm)+size(Ac,1)*tp+length(inv_nsigma_unm), 1);
else
    Z_unm = 10^11;
    X_unm = zeros(length(f_unm)+size(Ac,1)*tp+length(inv_nsigma_unm), 1);
    disp(RESULT.status)
end

end