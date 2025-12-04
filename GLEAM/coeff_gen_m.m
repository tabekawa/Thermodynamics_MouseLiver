function [F,SQ] = coeff_gen_m(H,f,Ac,bssmpl,exp,tp,H2O)

H = H(1:size(exp,1),1:size(exp,1));
fy = f(1:size(exp,1));
fy = fy./exp.*bssmpl;
for i = 0:tp-1
    fy(H2O+size(Ac,2)*i) = 0;
end
fa = [zeros(size(Ac,2),1);fy(size(Ac,2)*tp+1:size(exp,1))];
for j = 0:tp-1
    for i = 1:size(Ac,2)
        fa(i) =  fa(i) + fy(size(Ac,2)*j+i);
    end
end
F = [fy;fa];
SQ = bssmpl' * (H./2) * bssmpl;

end