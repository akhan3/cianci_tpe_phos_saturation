clear all;
% close all;
% clc;
% clf; drawnow
 
figure('windowStyle','docked');
ax1 = subplot(121); hold on; myplot; axis square;
ax2 = subplot(122); hold on; myplot; axis square;
drawnow;

numTauPoints = 10;
numPowerPoints = 20; % must be a mnimum of 5 for 'rat22' curve fitting model
excitationType = 'Sech2Pulse'; 
% excitationType = 'GaussianPulse'; 
% excitationType = 'CW'; 


% P = logspace(log10(1e-6), log10(1), numPowerPoints)';
TAU = logspace(log10(.1e-9), log10(10e-6), numTauPoints)';
% TAU = [[1:5:100]*1e-9, 1e-6, 10e-6, 100e-6]';

% TAU = 100e-9;

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

    P = 1.26e-6*sqrt(gamma) * logspace(-2,2,numPowerPoints)';

    %% Sanity check
    r = 0;    z = 0;
    assert(length(r) == 1);
    assert(length(z) == 1);
    quietStatus = true;



    %% Pre-allocate variables
    N1_ss = zeros(length(z), length(r), length(P));
    % P = sort(P, 'descend');
    if (~quietStatus && false) || true
        figure('windowStyle','docked');
        subplot(121); hold on;
    end

    %% Loop
    for iP = 1:length(P)
    %     disp( P(iP) )
        for ir = 1:length(r)
            for iz = 1:length(z)
                if ~quietStatus
                    subplot(121); hold on
                end
                N1_ss(iz,ir,iP) = cianci_pulseTrain(P(iP), gamma, excitationType, quietStatus);
                fprintf('%d,',   iP);
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
%     if length(y) ~= length(unique(y))
%         jj = min(find(diff(y)==0));
%     else
%         jj = length(y);
%     end
%     x = x(1:jj);
%     y = y(1:jj);
    
    %% curve fitting
    fop = fitoptions('rat22', 'Lower',[.5 0 0, 0 0], ...
                              'Upper',[.5 0 0, 0 inf], ...
                              'StartPoint', [0 0 0, 0 median(x1)^2]);
    [cfun,gof,fitAlgo] = fit(x1,y1, 'rat22',fop);
    Psat_fit = sqrt(cfun.q2);
%     interp1(y1,x1,.25)
    D = [P, squeeze(N1_ss)];


    %% make array
    PPsat_fit(kkk) = Psat_fit;
    PP{kkk} = P;
    NN1_ss{kkk} = squeeze(N1_ss);


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
    ph1 = loglog(P, squeeze(N1_ss), 'o');
    col = get(ph1,'color');            
    set(ph1, 'markerfacecolor', col);
    set(ph2,'color',col);
    plot(Psat_fit, .25,'sk', 'markerfacecolor', 'k');
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

    fprintf('\n');
    
    PPP(:,kkk) = P;
    NN1(:,kkk) = N1_ss;
    if sum(y1<0)
%         keyboard
    end
    
end % for kkk



PPsat_fit = PPsat_fit';
% [TAU PPsat_fit]

%% Curve fitting log-log
X = log(TAU);
Y = log(PPsat_fit);
fop = fitoptions('poly1', 'Lower',[-.5 -inf], 'Upper',[-.5 inf]);
[cfun,gof,fitAlgo] = fit(X,Y, 'poly1',fop);
A = exp(cfun.p2);

K = 2.2456e36; % [photon/J/m^2] % K = Sr/(2*h*c/lambda*f*(fwhm/1.763));


%% Curve fitting rat02
% fop = fitoptions('rat02', 'Lower',[-inf 0], ...
%                           'Upper',[.5 0 0, 0 inf], ...
%                           'StartPoint', [0 0 0, 0 median(x)^2]);
% [cfun,gof,fitAlgo] = fit(TAU,PPsat_fit, 'rat02');
% B = cfun.p1;
% x0 = sqrt(cfun.q1);


%% new figure    
axes(ax1); hold off; cla;

%% Plotting tau-Psat
cla; hold off;
loglog(TAU,A./sqrt(TAU),'r--');
hold on;
ph = loglog(TAU,PPsat_fit,'ob');
% plot(TAU,feval(cfun,TAU),'-.g');
hold off;
set(ph(:),'MarkerFaceColor','b');
xlabel('Probe Lifetime (s)')
ylabel('Saturation thereshold (W)')
str = sprintf('${%g}/{\\sqrt{\\tau}}$', A);
lh = legend({str,'Data'});
set(lh,'Interpreter','latex');
axis square

myplot

% close(1:10);

return

%% Plotting all P-N1 traces together
figure('windowStyle','docked', 'name','Combined');
clf
mark = {'s','d','^','v','o'};
legendStr = {'Pulsed Phosphorescence', 'Pulsed Fluorescence', 'CW Phosphorescence', 'CW Fluorescence'};
hold on
for k = 1:length(PPsat_fit)
    ph(k) = loglog(PP{k}, NN1_ss{k}, 'o-k', 'markerfacecolor','w');
%     set(ph(k), 'marker', mark{k});

    if PPsat_fit(k) <= 1
        th = text(PPsat_fit(k), .25, sprintf('%.1f mW',PPsat_fit(k)/1e-3));
    else
        th = text(PPsat_fit(k), .25, sprintf('%.1f W',PPsat_fit(k)));
    end
    th.BackgroundColor = 'w';
    th.Position = [th.Extent(1) th.Extent(2)*.5];
%     th2 = text(PP{k}(1), .0001, legendStr{k});
%     th2.Rotation = 59;
%     th2.FontWeight = 'bold';
%     th2.BackgroundColor = 'w';
end
plot(PPsat_fit, .25*ones(size(PPsat_fit)), 'ob', 'markersize',5, 'markerfacecolor','b','linewidth',2)

set(gca,'XScale', 'log');
set(gca,'YScale', 'log');
grid on; 
% %     ylim([1e-9 .5]); %axis square
xlabel('Average power (W)')
ylabel('Excited state population')
% axis square
% xlim([min(PP(:)) max(PP(:))])
% ylim([min(squeeze(N1_ss(:,:,1:iP))) 1])
% set(gca,'YTick', sort(unique([.25, .5, get(gca,'YTick')])) ); 
% title(['{\tau} = ', num2str(1/gamma)])
% legend(legendStr, 'location','southeast')

% lines at .5 and .25
a1 = plot(get(gca,'xlim'),0.5*[1 1], 'r--');
a2 = plot(get(gca,'xlim'),.25*[1 1], 'r--');
uistack([a1 a2], 'bottom')
hold off

myplot


return
