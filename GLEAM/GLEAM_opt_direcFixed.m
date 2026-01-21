function [X,F,RSS,Ps,resid,err] = GLEAM_opt_direcFixed(H,H_unm,f,f_unm,Ac,Af,C,ratio,lb,ub,ub_ob,m,sq,sq_unm,inv_nsigma,smpl,sigma,lambda,T,tp,d,f_exp,sq_exp)

rng('shuffle')
H_meas = H;
H = zeros(size(H,1),size(H,2));
H(1:size(Ac,2)*tp,1:size(Ac,2)*tp) = H(1:size(Ac,2)*tp,1:size(Ac,2)*tp) + H_unm;
H_unm = H;
f_meas = f;
f = zeros(size(f_exp,1),1);
f(1:size(Ac,2)*tp) = f(1:size(Ac,2)*tp) + f_unm;
f_unm = f;
inv_nsigma_unm = zeros(size(inv_nsigma,1),1);

S = Ac;
S(S>0) = 1;
S(S<0) = -1;
S = S(:,m);

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
Ss = repmat("S",1,tp);
Ss = join(Ss,",");
Ss = ["blkdiag(",Ss,")"];
Ss = join(Ss,"");
Ss = eval(Ss);

[X,F,RSS,Ps,resid,err] = opt_get_sol(H_meas,H_unm,f_meas,f_unm,Ac,Af,C,ratio,lb,ub,ub_ob,sq,sq_unm,inv_nsigma,inv_nsigma_unm,smpl,sigma,lambda,T,tp,d,Acs,Cs,Ss,f_exp,sq_exp);

end