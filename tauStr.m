function str = tauStr(tau)
    if tau < 100e-9
        str = sprintf('%g ns', tau*1e9);
    else
        str = sprintf('%g us', tau*1e6);
    end
