clear all;
close all;
clc;
% clf; drawnow
 
figure('windowStyle','docked');
ax1 = subplot(121); hold on; myplot; axis square;
ax2 = subplot(122); hold on; myplot; axis square;
drawnow;

numTauPoints = 5;
numPowerPoints = 12; % must be a mnimum of 5 for 'rat22' curve fitting model
excitationType = 'Sech2Pulse'; 
% excitationType = 'GaussianPulse'; 
% excitationType = 'CW'; 



TAU = logspace(log10(1e-9), log10(100e-6), 11)';
% TAU = [TAU; logspace(log10(1e-9), log10(10e-9), 5)'];
% TAU = [TAU; logspace(log10(10e-9), log10(10e-6), 5)'];
% TAU = [TAU; 50e-6; 100e-6; 500e-6; 1000e-6];


TAU = sort(unique(TAU));
% TAU = sort(TAU,'descend');


for kkk = 1:length(TAU)
    gamma = 1./TAU(kkk);
    tauNS = TAU(kkk) / 1e-9;
    if tauNS < 100
        tauStr = sprintf('%.1f ns', tauNS);
    else
        tauStr = sprintf('%.1f us', tauNS/1e3);
    end
    fprintf('%d/%d: TAU = %s:\t',   kkk, length(TAU), tauStr);

    P = 1.36e-6/sqrt(2.03e-8 + 1/gamma) * logspace(-1,1,numPowerPoints)';

    %% Sanity check
    r = 0;    z = 0;
    assert(length(r) == 1);
    assert(length(z) == 1);
    quietStatus = true;



    %% Pre-allocate variables
    N1_ss = zeros(size(P));
    % P = sort(P, 'descend');
    if (~quietStatus && false) || true
        figure('windowStyle','docked');
        subplot(121); hold on;
    end

    %% Loop
    for iP = 1:length(P)
        fprintf('\n\t[%2d/%2d]\t',   iP,length(P));
%         fprintf('\n\t[%2d] %.1e W\t',   iP,P(iP));
        for ir = 1:length(r)
            for iz = 1:length(z)
                if ~quietStatus
                    subplot(121); hold on
                end
                N1_ss(iP) = cianci_pulseTrain(P(iP), gamma, excitationType, quietStatus);
%                 fprintf('%d/%d: N1_ss(P=%g) = %f\n',   iP, length(P), P(iP), N1_ss(iz,ir,iP)); 

    % load xy
    % return            
                %% plotting in loop
                if ~quietStatus 
                    subplot(122)
                    hold on;
                    loglog(P(1:iP), squeeze(N1_ss(:,:,1:iP)), 'ob', 'markerfacecolor','w');
                    grid on; 
                    xlim([min(P) max(P)])
                    ylim([min(squeeze(N1_ss(:,:,1:iP))) 1])
                    set(gca,'YTick', sort(unique([.25, .5, get(gca,'YTick')])) ); 
                    str = sprintf('$\\tau$ = %s', tauStr);
                    title(str, 'Interpreter','latex')
                    axis square
                    myplot
                    drawnow
                end
            end
        end
    %     pause(1) 

    end


    % return




    %% Find saturation threshold
    x1 = P;
    y1 = squeeze(N1_ss);
    
    x1 = x1(~isnan(y1));
    y1 = y1(~isnan(y1));
    
%     if length(y) ~= length(unique(y))
%         jj = min(find(diff(y)==0));
%     else
%         jj = length(y);
%     end
%     x = x(1:jj);
%     y = y(1:jj);
    
    %% curve fitting with second-degree:             1/2
       %                                   y =  --------------
       %                                         1 + (x0/x)^2
    fop = fitoptions('rat22', 'Lower',[.5 0 0, 0 0], ...
                              'Upper',[.5 0 0, 0 inf], ...
                              'StartPoint', [0 0 0, 0 median(x1)^2]);
    [cfun,gof,fitAlgo] = fit(x1,y1, 'rat22',fop);
    Psat_fit_2 = sqrt(cfun.q2);

    %% curve fitting with general-exponent:            1/2
       %                                    y =  --------------
       %                                          1 + (x0/x)^m
    [xData, yData] = prepareCurveData( P, N1_ss );
    ft = fittype( '0.5/(1 + (x0/x)^m)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [0 0];
    opts.StartPoint = [2 0.616044676146639];
    opts.Upper = [10 12];
    [fr, gf] = fit( xData, yData, ft, opts );

    
    fprintf('\n');
%     fprintf('\tPsat = %G W\n', Psat_fit_2);
    fprintf('\t0.5/(1+(%.1e/P)^2), gof.rsquare=%g\n', Psat_fit_2, gof.rsquare);
%     fprintf('\t0.5/(1+(%f/P)^%f), gof.rsquare=%f\n', fr.x0, fr.m, gf.rsquare);

    %% make array
    PPsat_fit_2(kkk,1) = Psat_fit_2;
    PPsat_fit_m(kkk,1) = fr.x0;
    PPsat_mmm_m(kkk,1) = fr.m;

%     PP{kkk} = P;
%     NN1_ss{kkk} = squeeze(N1_ss);


    %% plotting N1-P saturation curve
%     figure('windowStyle','docked');
%     axes(ax2); hold on;
    subplot(122)
    myplot
    ph1 = loglog(P, squeeze(N1_ss), 'o');
    ph1 = loglog(P, squeeze(N1_ss), 'o');
    col = get(ph1,'color');            
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    grid on; 
    hold on;
    limX = get(gca,'xlim');
%     x = logspace(log10(min(limX)), log10(max(limX)), 100);
    x = logspace(log10(min(P)), log10(max(P)),100);
    y = feval(cfun,x);
    ph2 = plot(x,y,'-');;
    delete(ph1)
%     ph1 = loglog(P, squeeze(N1_ss), 'o');
    ph1 = loglog(x1, y1, 'o');
    col = get(ph1,'color');            
    set(ph1, 'markerfacecolor', col);
    set(ph2,'color',col);
    plot(Psat_fit_2, .25,'sk', 'markerfacecolor', 'k');
    plot(fr.x0, .25,'sr', 'markerfacecolor', 'r');
    ph2 = plot(x,feval(fr,x),'--');
    hold off;
%     lh = legend({'$\frac{1/2}{1+\left({P_{sat}}/{P}\right)^2}$ fitting'});
%     set(lh,'location','southeast');
%     set(lh,'Interpreter','latex','fontsize',18);
    xlabel('Average power (W)')
    ylabel('Excited state population @ (\rho=0,z=0)')
    axis square

    tauNS = (1/gamma) / 1e-9;
    if tauNS < 100
        str = sprintf('$\\tau$ = %.1f ns', tauNS);
    else
        str = sprintf('$\\tau$ = %.1f ${\\mu}$s', tauNS/1e3);
    end
    title(str, 'Interpreter','latex')

    myplot
    drawnow

    
    PPP(:,kkk) = P;
    NN1(:,kkk) = N1_ss;
    if sum(y1<0)
%         keyboard
    end
    
%     return
    
end % for kkk

% return


%% Curve fitting for lowpass filter (Bode plot):             B
   %                                             y =  ----------------
   %                                                   sqrt(1 + x/x0)
[xData, yData] = prepareCurveData(TAU, PPsat_fit_2);
ft = fittype( 'B/(1 + x/x0)^m', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 .5 0];
opts.StartPoint = [1 .5 median(TAU)];
opts.Upper = [inf .5 inf];
[fr, gf] = fit( xData, yData, ft, opts);
B = fr.B;
tau0 = fr.x0;
A = B*sqrt(tau0);

K = 2.2456e36; % [photon/J/m^2] % K = Sr/(2*h*c/lambda*f*(fwhm/1.763));


fprintf('B/sqrt(1+T/T0) = %.2e/sqrt(1+T/%.2e) = %.2e/sqrt(%.2e + T), gof.rsquare=%g\n', ...
            B, tau0, A, tau0, gf.rsquare);


%% new figure    
axes(ax1); hold off; cla;
clf

%% Plotting tau-Psat
cla; hold off;
xx = logspace(log10(min(TAU)), log10(max(TAU)), 100)';
loglog(xx, feval(fr,xx),'r-');
hold on;
loglog(xx, A./sqrt(xx),'--k');
loglog(xx, B*ones(size(xx)),'--k');
ph = loglog(TAU,PPsat_fit_2,'ob');
% ph2 = loglog(TAU,PPsat_fit_m,'sr');
hold off;
set(ph,'MarkerFaceColor','w');
% set(ph2,'MarkerFaceColor','w');
xlabel('Probe Lifetime, \tau (s)')
ylabel('Saturation onset, P_{sat} (W)')
str = sprintf('$\\frac{%.2E}{\\sqrt{1+\\tau/%.2E}} = \\frac{%.2E}{\\sqrt{%.2E+\\tau}} $', B, tau0, A, tau0 );
lh = legend({str});
set(lh,'Interpreter','latex','location','southwest','fontsize',18,'box','off');
% axis square

myplot

% close(1:10);

