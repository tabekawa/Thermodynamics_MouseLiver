function x = min_enzyme_cost(Ac,Af,CC,ratio,lb,ub,opt,rxn,T,tp,RT,MW,Km,d,x0,no_init,fixmet,h)

rng('shuffle')
d = d(1+h*size(Ac,1):(h+1)*size(Ac,1));
A = [Ac,Af*T] .* d;
lb(lb>log(10^-7)) = log(10^-7);
lb([16,56]) = 0;
ub(ub<log(10^-2)) = log(10^-2);
ub([16,56]) = 0;
lb(fixmet) = opt(fixmet+size(Ac,2)*h);
ub(fixmet) = opt(fixmet+size(Ac,2)*h);
ratio = ratio(1+size(CC,1)*h:size(CC,1)*(h+1));
rG_lim = repmat(-0.01,size(A,1),1);
rG_lim(26) = rG_lim(26) - RT*log(10^(d(26)*(7.2-7.4)));
rG_lim(29+24) = rG_lim(29+24) - RT*log(10^(d(29+24)*(7.2-7.4)));
rG_lim = rG_lim - Af*T.*d * opt(size(Ac,2)*tp+1:size(Ac,2)*tp+size(Af,2));
Km = Km .* d(1:size(Km,1));
Km = Km(rxn,:);
options = optimoptions('fmincon','OptimalityTolerance',10^(-7),'ConstraintTolerance',10^(-9),'MaxFunctionEvaluations',1000000,'MaxIterations',100000,'Display','final');

    function f = fun(xx)
        conc = exp(xx(1:size(Ac,2)));
        DrG = A(rxn,:) * [xx;opt(size(Ac,2)*tp+1:size(Ac,2)*tp+size(Af,2))];
        l = 1;
        for j = 1:size(Km,1)
            for i = 1:size(Km,2)
                if Km(j,i) < 0
                    if Km(j,i) ~= -1
                        cKm(l) = conc(i)/(-Km(j,i));
                        l = l+1;
                    end
                elseif Km(j,i) > 0
                    if Km(j,i) ~= 1
                        cKm(l) = conc(i)/Km(j,i);
                        l = l+1;
                    end
                end
            end
        end
        sat_unm = geomean(cKm);
        for j = 1:size(Km,1)
            sub_q = 1;
            prod_q = 1;
            for i = 1:size(Km,2)
                if Km(j,i) < 0
                    if Km(j,i) ~= -1
                        sub_q = sub_q * conc(i)/(-Km(j,i));
                    elseif Km(j,i) == -1
                        sub_q = sub_q * sat_unm;
                    end
                elseif Km(j,i) > 0
                    if Km(j,i) ~= 1
                        prod_q = prod_q * conc(i)/Km(j,i);
                    elseif Km(j,i) == 1
                        prod_q = prod_q * sat_unm;
                    end
                end
            end
            rev(j) = 1-exp(DrG(j)/RT);
            sat_s(j) = sub_q / (1 + sub_q + prod_q);
            sat_p(j) = prod_q / (1 + sub_q + prod_q);
            E(j) = MW(j) * 1/rev(j) * 1/sat_s(j);
        end
        f = sum(E);
    end

problem = createOptimProblem('fmincon','objective',@fun, ...
    'Aineq',Ac.*d,'bineq',rG_lim, ...
    'Aeq',CC,'beq',ratio,'lb',lb,'ub',ub,'options',options,'x0',x0);
ms = MultiStart;
ms.UseParallel = 'always';
[x,fval,exitflag,output,solutions] = run(ms,problem,no_init);
end