% clc
clear
% close all

verbose = false;

% TPA_GM = [1 10 30 100 200]';
TPA_GM = logspace(log10(1),log10(200),4)' ;

for jj = 1:length(TPA_GM)
    

figure('windowStyle','docked');
ax1 = subplot(121); hold on; myplot; axis square;
ax2 = subplot(122); hold on; myplot; axis square;
drawnow;

numTauPoints = 6;
numPowerPoints = 8; % must be a mnimum of 5 for 'rat22' curve fitting model
excitationType = 'Sech2'; 
% excitationType = 'Gaussian'; 
% excitationType = 'CW'; 

%% Fluorophore
tpaGM = TPA_GM(jj);
tpa = tpaGM * 1e-58;  % GM = 1e-58 m^4 / (photon/s)
clear tpaGM;

if strcmp(excitationType, 'CW')
    TAU = 1e-6 * [      .01, .1, 1, 4, 10, 60, 100]' ; % CW
else
    TAU = 1e-6 * [.001, .01, .1, 1, 4, 10, 60, 100]' ; 
end


TAU = sort(unique(TAU));
% TAU = sort(TAU,'descend');


for kkk = 1:length(TAU)
    gamma = 1./TAU(kkk);
    fprintf('%d/%d: TPA = %g GM:\t%d/%d: TAU = %s:\t', jj, length(TPA_GM), TPA_GM(jj),  kkk, length(TAU), tauStr(TAU(kkk)));

    if strcmp(excitationType, 'CW')
        P = 3.5336e-32 / sqrt(tpa) / sqrt(1/gamma) * logspace(-1,1,numPowerPoints)'; % CW
    else
        P = (1.35e-34/sqrt(tpa)) / sqrt(2.04e-8 + 1/gamma) * logspace(-1,1,numPowerPoints)';
    end
    
    % P = sort(P,'descend');


    %% Pre-allocate variables
    N1_ss = zeros(size(P));

    %% Create figure window for time traces
    if verbose
        figure('windowStyle','docked');
        myplot;
        hold on;
    end
    
    %% Sweep power values
    for iP = 1:length(P)
        fprintf('\n\t[%d/%d]\t',   iP,length(P));
        N1_ss(iP) = cianci_pulseTrain(P(iP), gamma, tpa, excitationType, verbose);
    end
    fprintf('\n');




    %% Find saturation threshold
    [xData, yData] = prepareCurveData(P, N1_ss);

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
    Ps = sqrt(cfun.q2);
    fprintf('\t0.5/(1+(%g/P)^2), gof.rsquare=%g\n', Ps, gof.rsquare);

    %% make array
    Psat(kkk,1) = Ps;
    Psat_gof(kkk,1) = gof.rsquare;
    

    %% plotting N1-P saturation curve
    axes(ax2); hold on;
    subplot(122)
    myplot
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    grid on; 
    hold on;
    xx = logspace(log10(min(P)), log10(max(P)),100);
    yy = feval(cfun,xx);
    ph1 = plot(xx,yy,'-');
    ph2 = loglog(xData, yData, 'o');
    col = get(ph1,'color');            
    set(ph2, 'markerfacecolor', col);
    set([ph1,ph2],'color',col);
    plot(Ps, .25,'sk', 'markerfacecolor', 'k');
    hold off;
    %     str = sprintf('$\\frac{1/2}{1+\\left({%G}/{P}\\right)^2}$ fitting', Psat)
    %     lh = legend(str);
    %     text(Ps*1.2,.25/1.2, sprintf('%.2E W',Ps));
    %     set(lh,'location','southeast','Interpreter','latex','fontsize',16);
    set(gca,'xlim',[1e-6 1e-0]);
    if strcmp(excitationType, 'CW')
        set(gca,'xlim',[1e-3 1e3]);
        set(gca,'xtick',[1e-3, 1e0, 1e3]);
    else
        set(gca,'ylim',[1e-3 1e-0]);
    end
    hold on
    plot(get(gca,'xlim'), .25*ones(size(xlim)), '--k');
    plot(get(gca,'xlim'), .50*ones(size(xlim)), '--k');
    hold off
    xlabel('Average power, <P> (W)')
    ylabel('Excited state population, N_1');   % @ (\rho=0,z=0)')
    title(sprintf('%s, \\sigma_{TPA} = %G GM', excitationType, tpa/1e-58));
    axis square
    myplot
    drawnow

    
    PPP(:,kkk,jj) = P;
    NN1(:,kkk,jj) = N1_ss;
    TTAU(kkk,jj) = TAU(kkk);
    
    
end % for kkk

% return

%% Seed the main curve fitting method with the power model % y = a*x^b
[xData, yData] = prepareCurveData(TAU, Psat);
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
[xData, yData] = prepareCurveData(TAU, Psat);
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
a = A^2*tpa;

% K = 2.2456e36; % [photon/J/m^2] % K = Sr/(2*h*c/lambda*f*(fwhm/1.763));


% fprintf('B/sqrt(1+T/T0) = %.2e/sqrt(1+T/%.2e) = %.2e/sqrt(%.2e + T), gof.rsquare=%g\n', ...
%             B, tau0, A, tau0, gf.rsquare);
fprintf('A/sqrt(T0+T) = %g/sqrt(%g + T), gof.rsquare=%g\n', ...
            A, tau0, gf.rsquare);


%% new figure    
axes(ax1); hold off; cla;
% clf

%% Plotting tau-Psat
cla; hold off;
xx = logspace(log10(min(TAU)), log10(max(TAU)), 100)';
loglog(xx, feval(fr,xx),'r-');
hold on;
% loglog(xx, A./sqrt(xx),'--k');
% loglog(xx, B*ones(size(xx)),'--k');
ph = loglog(TAU,Psat,'ob');
% ph2 = loglog(TAU,PPsat_fit_m,'sr');
hold off;
set(ph,'MarkerFaceColor','b');
% set(ph2,'MarkerFaceColor','w');
xlabel('Probe Lifetime, \tau (s)')
ylabel('Saturation onset, P_{sat} (W)')
% str = sprintf('$\\frac{%G}{\\sqrt{%G+\\tau}}$', A, tau0 );
str = sprintf('$\\sqrt{\\frac{%G/%G}{%G+\\tau}}$', a, tpa, tau0 );
lh = legend({str});
set(lh,'Interpreter','latex','location','southwest','fontsize',18,'box','off');
str = sprintf('%G GM, \\tau = [%s, %s]', tpa/1e-58, tauStr(min(TAU)), tauStr(max(TAU)));
title(str)
axis square
set(gca,'xlim',[1e-9 1e-4]);
if strcmp(excitationType, 'CW')
    set(gca,'ylim',[1e-2 1e2]);
else
    set(gca,'ylim',[1e-4 1e-1]);
end
myplot

% close(1:10);

%% printing data for saving
folderName = sprintf('FILES_%s', excitationType);
fileName = sprintf('Data %s %g GM %s - %s', excitationType, tpa/1e-58, tauStr(min(TAU)), tauStr(max(TAU)))
fh = fopen(fullfile(folderName, [fileName,'.txt']), 'w');

CW_fprintf(fh, '\n');
CW_fprintf(fh, '##################################################################\n')
CW_fprintf(fh, '\tTPA\t%g\t[GM]\t%g\t[m^4.s/photon]\n', tpa/1e-58, tpa);
CW_fprintf(fh, '\n')
CW_fprintf(fh, '\tTau [s]\tPsat [W]\tgof\n');
for kt = 1:length(TAU)
    CW_fprintf(fh, '\t%g\t%g\t%g\t\t%s\n', TAU(kt), Psat(kt), Psat_gof(kt), tauStr(TAU(kt)));  
end
CW_fprintf(fh, '\n');
CW_fprintf(fh, '\tA [W.sqrt(s)]\tTau_0 [s]\tgof\ta [J^2.m^4]\n');
CW_fprintf(fh, '\t%g\t%g\t%g\t%g\n', A, tau0, gf.rsquare, a);
CW_fprintf(fh, '##################################################################\n')
CW_fprintf(fh, '\n');

fclose(fh);

print(gcf, fullfile(folderName, fileName), '-dpdf');

save(fullfile(folderName, fileName));



end %% jj

