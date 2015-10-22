clear all;
close all;
clc;
% clf; drawnow
 
figure('windowStyle','docked');
ax1 = subplot(121); hold on; myplot; axis square;
ax2 = subplot(122); hold on; myplot; axis square;
drawnow;

numTauPoints = 6;
numPowerPoints = 8; % must be a mnimum of 5 for 'rat22' curve fitting model
% excitationType = 'Sech2Pulse'; 
% excitationType = 'GaussianPulse'; 
excitationType = 'CW'; 

%% Fluorophore
tpa = 1 * 1e-58;  % GM = 1e-58 m^4 / (photon/s)


TAU = logspace(log10(.1e-6), log10(100e-6), numTauPoints)';
% TAU = [TAU; logspace(log10(1e-9), log10(10e-9), 5)'];
% TAU = [TAU; logspace(log10(10e-9), log10(10e-6), 5)'];
% TAU = [TAU; 50e-6; 100e-6; 500e-6; 1000e-6];

TAU = 1e-6 * [.001, .01, .1, 1, 4, 10, 60, 100]' ;


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

    % P = 1.35e-6/sqrt(2.04e-8 + 1/gamma) * logspace(-1,1,numPowerPoints)';
    P = 1e2 * (1.35e-34/sqrt(tpa)) / sqrt(2.04e-8 + 1/gamma) * logspace(-1,1,numPowerPoints)';

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
                N1_ss(iP) = cianci_pulseTrain(P(iP), gamma, tpa, excitationType, quietStatus);
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
    fprintf('\n');


    % return




    %% Find saturation threshold
%     x1 = x1(~isnan(y1));
%     y1 = y1(~isnan(y1));
    [xData, yData] = prepareCurveData( P, squeeze(N1_ss) );

    %% curve fitting with second-degree:             1/2
       %                                   y =  --------------
       %                                         1 + (x0/x)^2
    fop = fitoptions('rat22', 'Lower',[.5 0 0, 0 0], ...
                              'Upper',[.5 0 0, 0 inf], ...
                              'StartPoint', [0 0 0, 0 median(xData)^2]);
    fop.TolFun = 1e-10;
    fop.TolX = 1e-10;
    fop.Display = 'final';
    [cfun,gof,fitAlgo] = fit(xData,yData, 'rat22',fop);
    Psat_fit_2 = sqrt(cfun.q2);

%     %% curve fitting with general-exponent:            1/2
%        %                                    y =  --------------
%        %                                          1 + (x0/x)^m
%     [xData, yData] = prepareCurveData( P, N1_ss );
%     ft = fittype( '0.5/(1 + (x0/x)^m)', 'independent', 'x', 'dependent', 'y' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'final';
%     opts.TolFun = 1e-10;
%     opts.TolX = 1e-10;
%     opts.Lower = [0 0];
%     opts.StartPoint = [2 0.616044676146639];
%     opts.Upper = [10 12];
%     [fr,gf,fitAlgo] = fit( xData, yData, ft, opts );

    
%     fprintf('\tPsat = %G W\n', Psat_fit_2);
    fprintf('\t0.5/(1+(%.1e/P)^2), gof.rsquare=%g\n', Psat_fit_2, gof.rsquare);
%     fprintf('\t0.5/(1+(%f/P)^%f), gof.rsquare=%f\n', fr.x0, fr.m, gf.rsquare);

    %% make array
    PPsat_fit_2(kkk,1) = Psat_fit_2;
%     PPsat_fit_m(kkk,1) = fr.x0;
%     PPsat_mmm_m(kkk,1) = fr.m;

%     PP{kkk} = P;
%     NN1_ss{kkk} = squeeze(N1_ss);


    %% plotting N1-P saturation curve
%     figure('windowStyle','docked');
%     axes(ax2); hold on;
    subplot(122)
    myplot
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    grid on; 
    hold on;
    xx = logspace(log10(min(P)), log10(max(P)),100);
    yy = feval(cfun,xx);
    ph1 = plot(xx,yy,'-r');
    ph2 = loglog(xData, yData, 'o');
    col = get(ph1,'color');            
    set(ph2, 'markerfacecolor', col);
    set([ph1,ph2],'color',col);
    plot(Psat_fit_2, .25,'sk', 'markerfacecolor', 'k');
%     plot(fr.x0, .25,'sr', 'markerfacecolor', 'r');
%     ph2 = plot(x,feval(fr,x),'--');
    hold off;
%     lh = legend({'$\frac{1/2}{1+\left({P_{sat}}/{P}\right)^2}$ fitting'});
%     set(lh,'location','southeast');
%     set(lh,'Interpreter','latex','fontsize',18);
    xlabel('Average power (W)')
    ylabel('Excited state population @ (\rho=0,z=0)')
    text(Psat_fit_2*1.2,.25/1.2, sprintf('%.2E W',Psat_fit_2));
    axis square

    tauNS = (1/gamma) / 1e-9;
    if tauNS < 100
        str = sprintf('TPA = %G GM, \\tau = %.1f ns', tpa/1e-58, tauNS);
    else
        str = sprintf('TPA = %G GM, \\tau = %.1f us', tpa/1e-58, tauNS/1e3);
    end
    title(str);

    myplot
    drawnow

    
    PPP(:,kkk) = P;
    NN1(:,kkk) = N1_ss;
    if sum(yData<0)
%         keyboard
    end
    
%     return
    
end % for kkk

% return

%% Seed the main curve fitting method with the power model % y = a*x^b
[xData, yData] = prepareCurveData(TAU, PPsat_fit_2);
ft = fittype( 'power1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.TolFun = 1e-10;
opts.TolX = 1e-10;
opts.Display = 'final';
opts.StartPoint = [1e-6 -.5];
[fr, gf] = fit( xData, yData, ft, opts );
A_start = fr.a;

%% Curve fitting for lowpass filter (Bode plot):             B
   %                                             y =  ----------------
   %                                                   sqrt(1 + x/x0)
[xData, yData] = prepareCurveData(TAU, PPsat_fit_2);
ft = fittype( 'A/sqrt(x0 + x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.TolFun = 1e-10;
opts.TolX = 1e-10;
opts.Display = 'final';
opts.Lower = [0 0];
opts.StartPoint = [A_start median(TAU)];
opts.Upper = [inf inf];
[fr,gf,fitAlgo] = fit( xData, yData, ft, opts);
A = fr.A;
tau0 = fr.x0;
% B = A/sqrt(tau0);

K = 2.2456e36; % [photon/J/m^2] % K = Sr/(2*h*c/lambda*f*(fwhm/1.763));


% fprintf('B/sqrt(1+T/T0) = %.2e/sqrt(1+T/%.2e) = %.2e/sqrt(%.2e + T), gof.rsquare=%g\n', ...
%             B, tau0, A, tau0, gf.rsquare);
fprintf('A/sqrt(T0+T) = %g/sqrt(%g + T), gof.rsquare=%g\n', ...
            A, tau0, gf.rsquare);


%% new figure    
axes(ax1); hold off; cla;
clf

%% Plotting tau-Psat
cla; hold off;
xx = logspace(log10(min(TAU)), log10(max(TAU)), 100)';
loglog(xx, feval(fr,xx),'r-');
hold on;
% loglog(xx, A./sqrt(xx),'--k');
% loglog(xx, B*ones(size(xx)),'--k');
ph = loglog(TAU,PPsat_fit_2,'ob');
% ph2 = loglog(TAU,PPsat_fit_m,'sr');
hold off;
set(ph,'MarkerFaceColor','b');
% set(ph2,'MarkerFaceColor','w');
xlabel('Probe Lifetime, \tau (s)')
ylabel('Saturation onset, P_{sat} (W)')
str = sprintf('$\\frac{A}{\\sqrt{\\tau_0+\\tau}} = \\frac{%G}{\\sqrt{%G+\\tau}}$', A, tau0 );
lh = legend({str});
set(lh,'Interpreter','latex','location','southwest','fontsize',18,'box','off');
str = sprintf('TPA = %G GM, \\tau = [%g,%g] s', tpa/1e-58, min(TAU), max(TAU));
title(str)
% axis square
%         set(gca,'xscale','linear');
%         set(gca,'yscale','linear');
myplot

% close(1:10);



