function [X,F,E,D,RSS,Ps,resid,O,Px,Py,err] = GLEAM_opt(H,H_unm,f,f_unm,Ac,Af,C,ratio,lb,ub,ub_ob,m,sq,sq_unm,irr_pos,irr_neg,inv_nsigma,smpl,sigma,lambda,T,tp,n,popsize)

rng('shuffle')
D = zeros(size(Ac,1)*tp,n);
F = repmat(10^8,n,1);
E = repmat(10^8,n,1);
Py = cell(n,1);
Px = zeros(n,1);
H_meas = H;
H = zeros(size(H,1),size(H,2));
H(1:size(Ac,2)*tp,1:size(Ac,2)*tp) = H(1:size(Ac,2)*tp,1:size(Ac,2)*tp) + H_unm;
H_unm = H;
f_meas = f;
f = zeros(size(f,1),1);
f(1:size(Ac,2)*tp) = f(1:size(Ac,2)*tp) + f_unm;
f_unm = f;
inv_nsigma_unm = zeros(size(inv_nsigma,1),1);

S = Ac;
S(S>0) = 1;
S(S<0) = -1;
S = S(:,m);
lb_d = zeros(size(Ac,1)*tp,1);
lb_d(irr_pos) = 1;
ub_d = ones(size(Ac,1)*tp,1);
ub_d(irr_neg) = 0;

% 糖新生→解糖疎外
glyc_const = zeros(tp-2,size(Ac,1)*tp);
for i = 1:tp/2-1
    glyc_const(i,size(Ac,1)*(i-1)+5) = -1;
    glyc_const(i,size(Ac,1)*i+5) = 1;
end
for i = tp/2:tp-2
    glyc_const(i,size(Ac,1)*i+5) = -1;
    glyc_const(i,size(Ac,1)*(i+1)+5) = 1;
end

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

for i = 1:n
    [x,fval,exitflag,output] = gaqp(H_meas,H_unm,f_meas,f_unm,Ac,Af,C,ratio,lb,ub,ub_ob,m,sq,sq_unm,inv_nsigma,inv_nsigma_unm,lambda,T,tp,Acs,Cs,Ss,lb_d,ub_d,glyc_const,popsize);
    D(:,i) = x';
    F(i) = fval;
    E(i) = exitflag;
    O(i)= output;
    h = findobj(gcf,'Tag','gaplotonlybest');
    BestF = repmat(10^8,1,length(h)-1);
    for k = 2:length(h)
        BestF(k-1) = h(k).YData;
    end
    Py{i} = BestF;
    Px(i) = length(h)-2;
end
[~,I] = min(F);
D = D(:,I);
D(D==0) = -1;
[X,F,RSS,Ps,resid,err] = opt_get_sol(H_meas,H_unm,f_meas,f_unm,Ac,Af,C,ratio,lb,ub,ub_ob,sq,sq_unm,inv_nsigma,inv_nsigma_unm,smpl,sigma,lambda,T,tp,D,Acs,Cs,Ss,f_meas,sq);
figure(2)
xlabel('Generation')
ylabel('Best Score')
hold on
for j = 1:n
    plot(Px(j):-1:0,Py{j},'b-');
end

end