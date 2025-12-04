function [X,F,E,O] = gaqp(H,H_unm,f,f_unm,Ac,Af,C,ratio,lb,ub,ub_ob,m,sq,sq_unm,inv_nsigma,inv_nsigma_unm,lambda,T,tp,Acs,Cs,Ss,lb_d,ub_d,glyc_const,popsize)

options = optimoptions('ga', 'PlotFcn', {@gaplotonlybest,@gaplotscores}, 'FunctionTolerance', 10^(-9), 'PopulationSize', popsize);
glyc = [1,3,5,6,7,8,9];
rG_lim = repmat(-0.01,size(Ac,1)*tp,1);

% lower-level problem for continuous variables
    function Z = qp_continuous(d)
        d = d(size(Ac,1)*tp+1:size(Ac,1)*tp*2);
        d = round(d);
        d(d==0) = -1;

        % 解糖糖新生を揃える
        for q = 0:tp-1
            d(glyc+size(Ac,1)*q) = repmat(d(5+size(Ac,1)*q),length(glyc),1);
        end

        % MCTの電位差考慮
        for r = 0:tp-1
            rG_lim(26+size(Ac,1)*r) = -0.01 - 2.578730581020227*log(10^(d(26+size(Ac,1)*r)*(7.2-7.4)));
            rG_lim(29+24+size(Ac,1)*r) = -0.01 - 2.578730581020227*log(10^(d(29+24+size(Ac,1)*r)*(7.2-7.4)));
        end
        % 計測代謝物推定
        [Z_meas, sol_meas] = opt_measured(H,f,Ac,Af,C,ratio,lb,ub,ub_ob,m,sq,inv_nsigma,lambda,T,tp,Acs,Cs,Ss,rG_lim,d);
        if Z_meas < 10
            % 未計測代謝物推定
            f_conc = f(1:size(Ac,2)*tp);
            lb_meas_fixed = [repmat(lb,14,1);sol_meas(length(f_conc)+1:length(f_conc)+size(Af,2)+length(inv_nsigma))];
            ub_meas_fixed = [repmat(ub,14/2,1);repmat(ub_ob,14/2,1);sol_meas(length(f_conc)+1:length(f_conc)+size(Af,2)+length(inv_nsigma))];
            lb_meas_fixed(f_conc~=0) = sol_meas(f_conc~=0);
            ub_meas_fixed(f_conc~=0) = sol_meas(f_conc~=0);
            Z_unm = opt_unmeasured(H_unm,f_unm,Ac,Af,C,ratio,lb_meas_fixed,ub_meas_fixed,m,sq_unm,inv_nsigma_unm,lambda,T,tp,Acs,Cs,Ss,rG_lim,d);
            Z = Z_meas + Z_unm*10^(-10);
        elseif Z_meas == 10
            Z = 10;
        end
    end

    function state = gaplotonlybest(options, state, flag)
        xlabel('Generation')
        ylabel('Best Score')
        hold on
        best = min(state.Score);
        plot(state.Generation,best,'b.','Tag','gaplotonlybest');
    end

[X,F,E,O] = ga(@qp_continuous, size(Ac,1)*tp*2, ...
    [diag(repmat(-1,size(Ac,1)*tp,1)),diag(repmat(10^3,size(Ac,1)*tp,1)); ...
    [diag(repmat(-1,size(Ac,1)*tp,1)),diag(repmat(10^3,size(Ac,1)*tp,1))].*(-1); ...
    zeros(size(glyc_const,1),size(Ac,1)*tp),glyc_const], ...
    [repmat(-10^2+10^3,size(Ac,1)*tp,1);repmat(-10^2,size(Ac,1)*tp,1);zeros(size(glyc_const,1),1)], ...
    [(Ss)',zeros(size(m,1)*tp,size(Ac,1)*tp)], zeros(size(m,1)*tp,1), ...
    [repmat(-10^3,size(Ac,1)*tp,1);lb_d], [repmat(10^3,size(Ac,1)*tp,1);ub_d], [], size(Ac,1)*tp+1:size(Ac,1)*tp*2, options);
X = X(size(Ac,1)*tp+1:size(Ac,1)*tp*2);
for i = 0:tp-1
    X(glyc+size(Ac,1)*i) = repmat(X(5+size(Ac,1)*i),length(glyc),1);
end
X = round(X);
end