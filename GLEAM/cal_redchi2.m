function [By,redchi2] = cal_redchi2(Ac,Af,CC,ratio,H,H_unm,F,F_unm,x,d,err,inv_nsigma,lb,ub,ub_ob,smpl,sigma,tp,T,lambda,m_steady,sq,sq_unm,no_sample)
% cal_redchi2 - calculate reduced chi-square using parametric bootstrap sampling
%
% Inputs:
%   Ac : The linear constraint matrix regarding thermodynamic law for metabolite concentrations
%   Af : The linear constraint matrix regarding thermodynamic law for the standard Gibbs free energy of formation
%   CC : The linear constraint matrix regarding cofactor ratios for metabolite concentrations
%   ratio : bounds vector for cofactor ratios
%   H : coefficient matrix of the quadratic term of residual sum of squares constituting the objective function for GLEAM
%   H_unm : coefficient matrix of the quadratic term for the unique estimation of unmeasured metabolite concentrations constituting the objective function for GLEAM
%   F : coefficient vetor of the linear term of residual sum of squares constituting the objective function for GLEAM
%   F_unm : coefficient matrix of the linear term for the unique estimation of unmeasured metabolite concentrations constituting the objective function for GLEAM
%   x : optimal solution of GLEAM other than binary variables
%   d : optimal solution of binary variables of GLEAM
%   err : residual sum of squares in GLEAM
%   inv_nsigma : coefficient vetor of the regularization term of the systematic error constituting the objective function for GLEAM
%   lb : lower bounds vector of metabolite concentrations
%   ub : upper bounds vector of metabolite concentrations
%   ub_ob : upper bounds vector of metabolite concentrations in ob/ob mouse
%   smpl : mean value vector of metabolite concentrations and the standard Gibbs free energy of formation transformed by transposed matrix of the eigenvectors of the variance-covariance matrix
%   sigma : standard variation vector of metabolite concentrations and the standard Gibbs free energy of formation transformed by transposed matrix of the eigenvectors of the variance-covariance matrix
%   tp : sum of the number of time points
%   T : matrix of the eigenvectors of the standard Gibbs free energy variance-covariance matrix
%   lambda : regularization parameter for the systematic error
%   m_steady : indices vector of metabolite concentrations assumed to be at steady state
%   sq : constant offset for residual sum of squares constituting the objective function for GLEAM
%   sq_unm : constant offset for the unique estimation of unmeasured metabolite concentrations constituting the objective function for GLEAM
%   no_sample : the number of parametric bootsrap samples
%
% Outputs:
%   By : optimal solution of GLEAM using parametric bootstrap samples as inputs
%   redchi2 : reduced chi-square at the regularization parameter

%% parametric bootstrap sampling for metabolite concentrations and the standard Gibbs free energy
rng('shuffle')
[BS,Bsd] = perturb(Ac,Af,x,smpl,sigma,F,no_sample*5/4,tp);

%% apply GLEAM to the parametric bootstrap samples
F_y = F(1:size(Ac,2)*tp+size(Af,2));
for j = 1:size(BS,2)
    [F_bs,sq_bs] = gen_coeff(H,F,Ac,BS(:,j),smpl,tp,[16,41+15]);
    [By(:,j),fval(j),RSS,Ps,resid,preEPE(j)] = GLEAM_opt_direcFixed(H,H_unm,F_bs,F_unm,Ac,Af,CC,ratio,lb,ub,ub_ob,m_steady,sq_bs,sq_unm,inv_nsigma,smpl,sigma,lambda,T,tp,d,F,sq);
end

%% calculate reduced chi-square
notopt = checkbootstrap(fval);
By(:,:,notopt) = [];
preEPE(notopt) = [];
By = By(:,:,1:no_sample);
preEPE = preEPE(1:no_sample);
EPE = sum(preEPE)/length(preEPE);
df = length(find(F_y~=0))/2./Bsd.*(EPE - err);
redchi2 = length(find(F_y~=0)) .* EPE ./ df;
By = squeeze(By);

end