function [N1_ss] = cianci_pulseTrain(P, gamma, excitationType, quietStatus)

% clear
% clc
% clf; drawnow
% close all
% figure

charCount = 0;


%% physical constants
h = 6.63e-34; % J.s
c = 3e8; % m/s

%% Fluorophore
tpa = 100 * 1e-58;  % GM = 1e-58 m^4 / (photon/s)
% gamma = 1 / 1e-6;

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


%% time variables
% timeConstant = 1/(2 * tpa * Sr^2 *(P / (h*c/lambda)).^2 + gamma);
% numCycles = 1 + ceil(5*timeConstant / (1/f))
% numCycles = 20


% tic;
% N1 = zeros([size(R), length(P)]);
% for iP = 1:length(P)
%     disp( P(iP) )

%% Pre-allocate variables
% N1_ss = zeros(length(z), length(r), length(P));
% N1_ss_exc = zeros(length(z), length(r), length(iC), 
% t(:,iz,ir,iC), 
% N1(:,iz,ir,iC), 
% gp(:,iC)]

%% Initial conditions
iC = 1;

%% Loop
while (true)
%     fprintf(1, repmat('\b',1,charCount)); % delete line before
%     charCount = fprintf('\t iC=%d', iC); 

    for ir = 1:length(r)
        for iz = 1:length(z)
            if iC == 1
                N1_initial = 1e-10;
            else
                N1_initial = N1(end,iz,ir,iC-1);
            end
%             [~, t(:,iz,ir,iC), N1(:,iz,ir,iC), gp(:,iC)] = cianci_model(P, lambda, f, alfa, Sr(iz,ir), tpa, gamma, N1_initial, true);
            [t_ss_exc(iz,ir,iC), N1_ss_exc(iz,ir,iC), t(:,iz,ir,iC), N1(:,iz,ir,iC), gp(:,iC)] = cianci_model(P, lambda, f, fwhm, Sr(iz,ir), tpa, gamma, N1_initial, excitationType, true);
            t(:,iz,ir,iC) = t(:,iz,ir,iC) + (iC-1)/f; 
            if (~quietStatus && false)
                fprintf('N1_ss_exc(%d,%d,%d) = %G\n',             ir,iz,iC, N1_ss_exc(iz,ir,iC)); 
    %             fprintf('N1   (%d,%d,%d,%d) = %G\n',   size(t,1),ir,iz,iC, N1(end,iz,ir,iC)); 
            end
        end
    end
    
    if ~quietStatus && false
%         cla
        ph = plot(  squeeze(t(:)), .6*squeeze(gp(:)),'o-r',...
                        squeeze(t(:)), squeeze(N1(:)),'o-b',...
                        squeeze(t_ss_exc(end)) + (iC-1)/f, squeeze(N1_ss_exc(end)),'sg'); 
%         set(gca,'xscale', 'log');
        grid on
%         set(ph, 'linewidth', 2);
%         ylim([-.05 .65])
        str = sprintf('P = %.2f mW', P*1e3); title(str);
%         myplot;
        drawnow;
%         keyboard
    end

    N1_ss = N1_ss_exc(z==0,r==0,iC);
    if sum(isnan(N1_ss ) | isinf(N1_ss ) | N1_ss < 0)
        fprintf('\tSS_exc value is negative!!!!');
%         keyboard;
    end


    %% check the asymptote
    x(iC) = squeeze(t_ss_exc(:,:,iC)) + (iC-1)/f;
    y(iC) = squeeze(N1_ss_exc(:,:,iC));
    if iC >= 2
        dydx = diff(y)./diff(x);
        if dydx(end)/dydx(1) < exp(-5)  % == 0
            fprintf('#%d (%G, %G)', iC, P, N1_ss);
%             plot(x(1:end-1), 0.5*dydx./dydx(1), '-');
            break
        end
%         if  N1(end) <= 0 
%             fprintf('\tDecay asymptote goes negative! ');
%             break
%         end
        if dydx(1) == 0
            fprintf('#%d [zero-slope] (%G, %G)', iC, P, N1_ss);
            break;
        end
    end

%     if abs(1-N1_ss/N1_ss_prev) < 0.000136487992059  % 1e-3
%         break
%     end

%     if iC > 10
%         [Ass,timeConst] = fitRisingExp(x', y');
%     end
%     elseif isnan(N1_ss)
%         fprintf('ERROR: Power too high or bad parameters\n')
%         break
%     end
    
    iC = iC+1;
    
%     if iC == 10
%         break;
%     end
    
end %% while



if sum(isnan(N1_ss ) | isinf(N1_ss ) | N1_ss < 0)
    fprintf('\tSS is negative!!!!');
%     N1_ss = 0;
    %         keyboard;
end

% fprintf('\n');

fprintf(1, repmat('\b',1,charCount)); % delete line before
%     charCount = fprintf('\tiC = %d', iC);

% return


%% Plotting
% cla
% ph = plot(  squeeze(t(:)), .6*squeeze(gp(:)),'-r',...
% ph = plot(  ...  %squeeze(t(:)), .5*squeeze(gp(:)),'-',...
%                 squeeze(t(:)), squeeze(N1(:)),'-b'); grid on

ph = plot(  ... %squeeze(t(:)), .6*squeeze(gp(:)),'-r',...
                squeeze(t(:)), squeeze(N1(:)),'-b',...
                t_ss_exc(:,:,end) + (iC-1)/f, N1_ss_exc(:,:,end),'o');
col = get(ph(end),'color');            
set(ph(end), 'markerfacecolor', col);
set(ph(1), 'color', col);
grid on

% set(gca,'yscale','log')            
% ylim([-.05 .65])
% ylim([.0 .5])
% set(ph, 'linewidth', 2);
% set(gca, 'xtickLabel',    get(gca, 'xtick') / 1e-6);
xlabel('Time (s)');
    tauNS = (1/gamma) / 1e-9;
    if tauNS < 100
        tauStr = sprintf('%.1f ns', tauNS);
    else
        tauStr = sprintf('%.1f us', tauNS/1e3);
    end
str = sprintf('$P_{avg}$ = %.2f mW (%s) , %s', P*1e3, excitationType, tauStr); 
title(str,'interpreter','latex');
axis square 
myplot
drawnow;
% pause(.5)

t0 = fwhm*10;
tCondensed = sort(unique([t0, [t0-fwhm*3 : fwhm/12 : t0+fwhm*3]]));

x = squeeze(t(1,:,:,:)) + tCondensed(end); 
y = squeeze(N1_ss_exc(:));
W = tpa * ( Sr * (P / (h*c/lambda)) )^2;
timeConst = 1 / (2*W + gamma);
% [P   W    (2*W + gamma)  5*timeConst    t(end,end,end,end)    N1(end,end,end,end)      W/(2*W+gamma)    N1_ss_exc(1)]'
% save('xy',   'x','y');

return   















% fprintf('%d\n', iC);
% [N1_ss_prev, N1_ss, abs(1-N1_ss/N1_ss_prev),   abs(1-N1_ss/N1_ss_prev) < 0.000136487992059]
N1_00 = squeeze(N1(:,z==0,r==0,:));

% toc;




x = squeeze(t(:));
y = squeeze(N1(:));
yy = max(y)-y;

% title([Ntime, N1_00])
% keyboard 
return

end




function [Ass,timeConst] = fitRisingExp(x, y)
    [cfun,gof,fitAlgo] = fit(x,y, 'exp2');
    A = cfun.a;
    ta = -1/cfun.b;
    B = cfun.c;
    tb = -1/cfun.d;

    if A > 0 && B < 0
        Ass = A;
        timeConst = tb;
    elseif A < 0 && B > 0
        Ass = B;
        timeConst = ta;
    else
        error('Curve fitting failed');
    end

    % yfitted = Ass*(1-exp(x/timeConst));
    % hold off;
    % plot(cfun, x,y);
    % hold on;
    % plot(x,yfitted());
    % hold off;

end
