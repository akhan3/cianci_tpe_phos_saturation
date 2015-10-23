function [N1_ss] = cianci_pulseTrain(P, gamma, tpa, excitationType, doPlotting)

%% Excitation source
lambda = 780e-9; % 800 nm
f = 80e6;       % Hz
fwhm = 100e-15; % s

%% Gaussian beam PSF
w0 = .35e-6; % normalized to beam waist
r = 0;    z = 0;
[R,Z] = meshgrid(r,z);
z_rayleigh = pi * w0^2 / lambda;
w = w0 * sqrt(1 + (Z./z_rayleigh).^2);
S_gaussian = (2/pi)./(w.^2) .* exp(-2 * (R./w).^2);
Sr = S_gaussian;
assert(numel(Sr) == 1);


%% Loop until SS is reached
iC = 0; % Initial condition
while (true) 
    iC = iC+1;
    if iC == 1
        N1_initial = 1e-10;
        x = [0 0]; y = x; % prepare two element buffer for derivative
    end
    

    if doPlotting
        [t_ss_exc(iC), N1_ss_exc(iC), t(:,iC), N1(:,iC), pulse(:,iC)] = cianci_model(P, lambda, f, fwhm, Sr(1,1), tpa, gamma, N1_initial, excitationType);
        t(:,iC) = t(:,iC) + (iC-1)/f;
        N1_initial = N1(end, iC);
        N1_ss = N1_ss_exc(iC);
    else
        [t_ss_exc, N1_ss_exc, t, N1, ~] = cianci_model(P, lambda, f, fwhm, Sr(1,1), tpa, gamma, N1_initial, excitationType);
        t = t + (iC-1)/f;
        N1_initial = N1(end);
        N1_ss = N1_ss_exc;
    end


    %% check the asymptote
    x = [x(2), t_ss_exc(end) + (iC-1)/f];
    y = [y(2), N1_ss_exc(end)];
    if iC == 2
        dydx_1 = diff(y) / diff(x);
    end
    if iC >= 2
        dydx = diff(y) / diff(x);
        if dydx/dydx_1 < exp(-5) % must have passed 5 time constants (within <1% the SS value)
            fprintf('#%5d# [%.2f] (P = %s, N1ss = %f)', iC, log(dydx/dydx_1), PStr(P), N1_ss);
            break
        end
        if dydx_1 == 0
            fprintf('#%5d# [ZERO] (P = %s, N1ss = %.2f)', iC, PStr(P), N1_ss);
            break;
        end
    end
end % while


if N1_ss < 0 || N1_ss > 0.5 
    fprintf('\tWARNING: Abnormal N1_ss\t');
    keyboard;
end

%% Plotting
if doPlotting
    ph = plot(  ... % t(:), 0.6*pulse(:),'-r',...
                t(:), N1(:),'-b',...
                t_ss_exc(end) + (iC-1)/f, N1_ss,'o');
    col = get(ph(end),'color');            
    set(ph(end), 'markerfacecolor', col);
    set(ph(1), 'color', col);
    grid on
    ylim([.0 .5])
    set(gca, 'xtickLabel',    get(gca, 'xtick') / 1e-6);
    xlabel('Time (us)');
    str = sprintf('%G GM, %s, %s', tpa/1e-58, tauStr(1/gamma), PStr(P)); 
    title(str);
    myplot
    drawnow;
end
