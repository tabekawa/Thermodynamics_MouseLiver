function [X,F,RSS,Ps,resid,result,err] = opt_gqp_t_m(H,H_unm,f,Ac,Af,C,ratio,lb,ub,sq,inv_nsigma,smpl,sigma,lambda,T,tp,d,f_exp,sq_exp)

    H_pre = H;
    H_unm = H_unm .* 0;
    H(1:size(Ac,2)*tp,1:size(Ac,2)*tp) = H(1:size(Ac,2)*tp,1:size(Ac,2)*tp) + H_unm;

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
    
    rG_lim = repmat(-0.01,size(Ac,1)*tp,1);
    for r = 0:tp-1
        rG_lim(26+size(Ac,1)*r) = -0.01 - 2.578730581020227*log(10^(d(26+size(Ac,1)*r)*(7.2-7.4)));
        rG_lim(29+24+size(Ac,1)*r) = -0.01 - 2.578730581020227*log(10^(d(29+24+size(Ac,1)*r)*(7.2-7.4)));
    end

    model.Q = sparse([H,zeros(size(H,1),size(inv_nsigma,1));zeros(size(inv_nsigma,1),size(H,2)+size(inv_nsigma,1))])./2;
    model.obj = [f;(inv_nsigma.*lambda)];
    model.A = sparse([[Acs,repmat(Af*T,tp,1),zeros(size(Ac,1)*tp,size(inv_nsigma,1)*2)].*d; ...
        Cs,zeros(size(C,1)*tp,size(H,1)-size(C,2)*tp+size(inv_nsigma,1));...
        zeros(size(inv_nsigma,1),size(H,1)-size(inv_nsigma,1)),[diag(ones(size(Ac,2),1)),zeros(size(Ac,2),size(inv_nsigma,1)-size(Ac,2));zeros(size(inv_nsigma,1)-size(Ac,2),size(Ac,2)),T].*(-1),diag(repmat(-1,size(inv_nsigma,1),1)); ...
        zeros(size(inv_nsigma,1),size(H,1)-size(inv_nsigma,1)),[diag(ones(size(Ac,2),1)),zeros(size(Ac,2),size(inv_nsigma,1)-size(Ac,2));zeros(size(inv_nsigma,1)-size(Ac,2),size(Ac,2)),T],diag(repmat(-1,size(inv_nsigma,1),1))]);
    model.rhs = [rG_lim;repmat(ratio,tp,1);zeros(size(inv_nsigma,1)*2,1)];     
    model.lb = [lb;zeros(size(inv_nsigma,1),1)];
    model.ub = [ub;repmat(10^4,size(inv_nsigma,1),1)];
    model.objcon = sq;

    param.BarHomogeneous = 1;
    param.OptimalityTol = 10^(-9);
    param.FeasibilityTol = 10^(-9);
    param.NumericFocus = 3;

    RESULT = gurobi(model,param);

    result = RESULT.status;
    if strcmp(RESULT.status, 'INFEASIBLE') || strcmp(RESULT.status, 'INF_OR_UNBD') || strcmp(RESULT.status, 'NUMERIC')
        X = [lb;zeros(size(inv_nsigma,1),1)];
        F = 0;
        err = 0;
        RSS = 0;
        Ps = 0;
        resid = 0;
    else
        X = RESULT.x;
        F = RESULT.objval;
        err = X(1:size(f,1))' * (H_pre./2) * X(1:size(f,1)) + f_exp' * X(1:size(f,1)) + sq_exp;
        fc = f(1:size(Ac,2)*tp);
        ff = f(size(Ac,2)*tp+1:size(Ac,2)*tp+size(Af,2));
        fa = f(size(Ac,2)*tp+size(Af,2)+1:size(Ac,2)*tp+size(Af,2)+size(inv_nsigma,1));
        f = f(1:size(Ac,2)*tp+size(Af,2));
        fa_effec = [];
        for i = 0:tp-1
            single_fc = fc(1+size(Ac,2)*i:size(Ac,2)+size(Ac,2)*i);
            fa_effec = [fa_effec;find(single_fc~=0)];
        end
        fa_effec = [fa_effec;size(Ac,2)+find(ff~=0)];
        ya = X(f~=0) + X(size(f,1)+fa_effec) - smpl(f~=0);
        as = X(size(f,1)+find(fa~=0));
        RSS = (ya./sigma(f~=0))' * (ya./sigma(f~=0));
        Ps = length(as(-10^(-6)>=as | as>=10^(-6)));
        resid = ya./sigma(f~=0);
    end
    
end