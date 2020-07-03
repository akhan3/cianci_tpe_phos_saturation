function str = PStr(P)
    if P >= 1
        str = sprintf('%g W', P);
    elseif P <= 1000e-3 && P >= .1e-3 
        str = sprintf('%g mW', P*1e3);
    elseif P < .1e-3
        str = sprintf('%g uW', P*1e6);
    else 
        keyboard
        error('P out of range');
    end
