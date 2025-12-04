function [EC_tot, E, rev, sat_s, sat_p, kcat_d, sat_unm, G, cKm] = enzyme_cost_invivo(Ac,Af,x,rxn,T,RT,MW,Km,Km_micofactors,kcat,d,tp,h)
d = d(1+h*size(Ac,1):(h+1)*size(Ac,1));
A = [Ac,Af*T] .* d;
conc = exp(x(1+h*size(Ac,2):(h+1)*size(Ac,2)));
DrG = A(rxn,:) * [x(1+h*size(Ac,2):(h+1)*size(Ac,2));x(size(Ac,2)*tp+1:size(Ac,2)*tp+size(Af,2))];
d = d(1:size(Km_micofactors,1));
Km_micofactors = Km_micofactors .* d;
Km_micofactors = Km_micofactors(rxn,:);
kcat_d = 1;
l = 1;
for j = 1:size(Km_micofactors,1)
    for i = 1:size(Km_micofactors,2)
        if Km_micofactors(j,i) < 0
            if Km_micofactors(j,i) ~= -1
                cKm(l) = conc(i)/(-Km_micofactors(j,i));
                l = l+1;
            end
        elseif Km_micofactors(j,i) > 0
            if Km_micofactors(j,i) ~= 1
                cKm(l) = conc(i)/Km_micofactors(j,i);
                l = l+1;
            end
        end
    end
end
sat_unm = geomean(cKm);
for j = 1:size(Km_micofactors,1)
    sub_q = 1;
    prod_q = 1;
    for i = 1:size(Km_micofactors,2)
        if Km_micofactors(j,i) < 0
            if Km_micofactors(j,i) ~= -1
                sub_q = sub_q * conc(i)/(-Km_micofactors(j,i));
            elseif Km_micofactors(j,i) == -1
                sub_q = sub_q * sat_unm;
            end
        elseif Km_micofactors(j,i) > 0
            if Km_micofactors(j,i) ~= 1
                prod_q = prod_q * conc(i)/Km_micofactors(j,i);
            elseif Km_micofactors(j,i) == 1
                prod_q = prod_q * sat_unm;
            end
        end
    end
    rev(j) = 1-exp(DrG(j)/RT);
    sat_s(j) = sub_q / (1 + sub_q + prod_q);
    sat_p(j) = prod_q / (1 + sub_q + prod_q);
    E(j) = MW(j) * 1/rev(j) * 1/sat_s(j);
end
EC_tot = sum(E);
end