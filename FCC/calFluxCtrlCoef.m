function [FCC, varFCC, elas_s, elas_p] = calFluxCtrlCoef(rev,sat_s,sat_p,rxn,pathway)
elas_s = rev.^(-1) - sat_s;
elas_p = 1 - rev.^(-1) - sat_p;
if pathway == 1
    FCC(1) = prod(elas_s(2:end));
    for i = 2:5
        FCC(i) = prod(-elas_p(1:i-1)) * prod(elas_s(i+1:end));
    end
    FCC(6) = prod(-elas_p(1:4)) * (elas_s(5)-elas_p(5)) * prod(elas_s(7:end));
    for i = 7:length(rxn)-1
        FCC(i) = prod(-elas_p(1:4)) * (elas_s(5)-elas_p(5)) * prod(-elas_p(6:i-1)) * prod(elas_s(i+1:end));
    end
    FCC(length(rxn)) = prod(-elas_p(1:4)) * (elas_s(5)-elas_p(5)) * prod(-elas_p(6:length(rxn)-1));
    FCC = FCC./sum(FCC);
    FCC = FCC';
end

if pathway == 2
    FCC(length(rxn)) = prod(-elas_p(1:end-1));
    for i = 6:length(rxn)-1
        FCC(i) = prod(-elas_p(1:i-1)) * prod(elas_s(i+1:end));
    end
    FCC(5) = prod(-elas_p(1:4)) * (elas_s(6)-elas_p(6)) * prod(elas_s(7:end));
    for i = 2:4
        FCC(i) = prod(-elas_p(1:i-1)) * prod(elas_s(i+1:5)) * (elas_s(6)-elas_p(6)) * prod(elas_s(7:end));
    end
    FCC(1) = prod(elas_s(2:5)) * (elas_s(6)-elas_p(6)) * prod(elas_s(7:end));
    FCC = FCC./sum(FCC);
    FCC = FCC';
end

if pathway == 3
    FCC(1) = prod(elas_s(2:end));
    for i = 2:length(rxn)-1
        FCC(i) = prod(-elas_p(1:i-1)) * prod(elas_s(i+1:end));
    end
    FCC(length(rxn)) = prod(-elas_p(1:length(rxn)-1));
    FCC = FCC./sum(FCC);
    FCC = FCC';
end

varFCC = var(FCC);
end