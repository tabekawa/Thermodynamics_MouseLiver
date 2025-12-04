function notopt = checkbootstrap(Z)
    notopt = [];
        for j = 1:lenght(Z)
            if Z(j) >= 10
            notopt = [notopt;j];
            disp(Z(j))
            end
        end
end