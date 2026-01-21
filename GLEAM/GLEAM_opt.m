function [X,D,err] = GLEAM_opt(H,H_unm,f,f_unm,Ac,Af,CC,ratio,lb,ub,ub_ob,m_steady,sq,sq_unm,irr_pos,irr_neg,inv_nsigma,smpl,sigma,lambda,T,tp,n,popsize)
% GLEAM_opt - estimate thermodynamically consistent metabolite concentrations and the standard Gibbs free energy of formation by GLEAM
%
% Inputs:
%   H : coefficient matrix of the quadratic term of residual sum of squares constituting the objective function for GLEAM
%   H_unm : coefficient matrix of the quadratic term for the unique estimation of unmeasured metabolite concentrations constituting the objective function for GLEAM
%   f : coefficient vetor of the linear term of residual sum of squares constituting the objective function for GLEAM
%   f_unm : coefficient matrix of the linear term for the unique estimation of unmeasured metabolite concentrations constituting the objective function for GLEAM
%   Ac : The linear constraint matrix regarding thermodynamic law for metabolite concentrations
%   Af : The linear constraint matrix regarding thermodynamic law for the standard Gibbs free energy of formation
%   CC : The linear constraint matrix regarding cofactor ratios for metabolite concentrations
%   ratio : bounds vector for cofactor ratios
%   lb : lower bounds vector of metabolite concentrations
%   ub : upper bounds vector of metabolite concentrations
%   ub_ob : upper bounds vector of metabolite concentrations in ob/ob mouse
%   m_steady : indices vector of metabolite concentrations assumed to be at steady state
%   sq : constant offset for residual sum of squares constituting the objective function for GLEAM
%   sq_unm : constant offset for the unique estimation of unmeasured metabolite concentrations constituting the objective function for GLEAM
%   irr_pos : indices vector of irreversible reactions assumed to proceed in forward direction
%   irr_neg : indices vector of irreversible reactions assumed to proceed in backward direction
%   inv_nsigma : coefficient vetor of the regularization term of the systematic error constituting the objective function for GLEAM
%   smpl : mean value vector of metabolite concentrations and the standard Gibbs free energy of formation transformed by transposed matrix of the eigenvectors of the variance-covariance matrix
%   sigma : standard variation vector of metabolite concentrations and the standard Gibbs free energy of formation transformed by transposed matrix of the eigenvectors of the variance-covariance matrix
%   lambda : regularization parameter for the systematic error
%   T : matrix of the eigenvectors of the standard Gibbs free energy variance-covariance matrix
%   tp : sum of the number of time points
%   n : the number of GLEAM iterations
%   popsize : population size for th egenetic algorithm
%
% Outputs:
%   X : optimal solution of GLEAM other than binary variables
%   D : optimal solution of binary variables of GLEAM
%   err : residual sum of squares in GLEAM

%% prepare inputs for GLEAM
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
S = S(:,m_steady);
lb_d = zeros(size(Ac,1)*tp,1);
lb_d(irr_pos) = 1;
ub_d = ones(size(Ac,1)*tp,1);
ub_d(irr_neg) = 0;

% prohibit the shift from gluconeogenesis to glycolysis
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
CCs = repmat("CC",1,tp);
CCs = join(CCs,",");
CCs = ["blkdiag(",CCs,")"];
CCs = join(CCs,"");
CCs = eval(CCs);
Ss = repmat("S",1,tp);
Ss = join(Ss,",");
Ss = ["blkdiag(",Ss,")"];
Ss = join(Ss,"");
Ss = eval(Ss);

%% estimate metabolite concentrations and the standard Gibbs free energy of formation by GLEAM
for i = 1:n
    [x,fval,exitflag,output] = gaqp(H_meas,H_unm,f_meas,f_unm,Ac,Af,CC,ratio,lb,ub,ub_ob,m_steady,sq,sq_unm,inv_nsigma,inv_nsigma_unm,lambda,T,tp,Acs,CCs,Ss,lb_d,ub_d,glyc_const,popsize);
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
[X,F,RSS,Ps,resid,err] = opt_get_sol(H_meas,H_unm,f_meas,f_unm,Ac,Af,CC,ratio,lb,ub,ub_ob,sq,sq_unm,inv_nsigma,inv_nsigma_unm,smpl,sigma,lambda,T,tp,D,Acs,CCs,Ss,f_meas,sq);
figure(2)
xlabel('Generation')
ylabel('Best Score')
hold on
for j = 1:n
    plot(Px(j):-1:0,Py{j},'b-');
end

end