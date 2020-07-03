clear all;
clc;
% clf; drawnow
 

numPowerPoints = 20;

for kkk = [1,2,3,4]

%% Sweep power
switch kkk
    case 1  % pulsed phos 
        excitationType = 'Sech2Pulse'; gamma = 1/1e-6; P = logspace(log10(.02e-3), log10(100e-3), numPowerPoints)'; 
    case 2  % pulsed fluor 
        excitationType = 'Sech2Pulse'; gamma = 1/1e-9; P = logspace(log10(.2e-3), log10(1000e-3), numPowerPoints)'; 
    case 3  % CW phos 
        excitationType = 'CW'; gamma = 1/1e-6; P = logspace(log10(3e-3), log10(10), numPowerPoints)'; 
    case 4  % CW fluor 
        excitationType = 'CW'; gamma = 1/1e-9; P = logspace(log10(100e-3), log10(300), numPowerPoints)';  %% upto 300W
    otherwise
        error('Invalid choice of kkk');
end
 
% P = sort(P,'descend');





%% Sanity check
r = 0;    z = 0;
assert(length(r) == 1);
assert(length(z) == 1);
quietStatus = true;


%% Pre-allocate variables
N1_ss = zeros(length(z), length(r), length(P));
% P = sort(P, 'descend');
figure('windowStyle','docked');

%% Loop
for iP = 1:length(P)
%     disp( P(iP) )
    for ir = 1:length(r)
        for iz = 1:length(z)
            subplot(121); hold on
            N1_ss(iz,ir,iP) = cianci_pulseTrain(P(iP), gamma, excitationType, quietStatus);
            fprintf('%d/%d: N1_ss(P=%g) = %f\n',   iP, length(P), P(iP), N1_ss(iz,ir,iP)); 

% load xy
% return            
            %% plotting in loop
            subplot(122)
            loglog(P(1:iP), squeeze(N1_ss(:,:,1:iP)), 'o-b', 'markerfacecolor','w');
            grid on; 
            xlim([min(P) max(P)])
            ylim([min(squeeze(N1_ss(:,:,1:iP))) 1])
            set(gca,'YTick', sort(unique([.25, .5, get(gca,'YTick')])) ); 
            title(['{\tau} = ', num2str(1/gamma)])
            axis square
            drawnow
        end
    end
%     pause(1) 

end

  
% return

%% plotting
% figure('windowStyle','docked')
subplot(122)
% ah1 = subplot(121);
loglog(P, squeeze(N1_ss), 'o-b', 'markerfacecolor','w');
grid on; 
%     ylim([1e-9 .5]); %axis square
xlabel('Average power (W)')
ylabel('Excited state population @ (\rho=0,z=0)')
axis square
xlim([min(P) max(P)])
ylim([min(squeeze(N1_ss(:,:,1:iP))) 1])
set(gca,'YTick', sort(unique([.25, .5, get(gca,'YTick')])) ); 
title(['{\tau} = ', num2str(1/gamma)])
%     Plot(gcf);



%% Find saturation threshold
x = P;
y = squeeze(N1_ss);
if length(y) ~= length(unique(y))
    jj = min(find(diff(y)==0));
else
    jj = length(y);
end
x = x(1:jj);
y = y(1:jj);
Psat = interp1(y, x, .25)
D = [P, squeeze(N1_ss)];


%% make array
PPsat(kkk) = Psat;
PP{kkk} = P;
NN1_ss{kkk} = squeeze(N1_ss);



end % for kkk


figure('windowStyle','docked', 'name','Combined');

%% Plotting all together
clf
mark = {'s','d','^','v','o'};
legendStr = {'Pulsed Phosphorescence', 'Pulsed Fluorescence', 'CW Phosphorescence', 'CW Fluorescence'};
hold on
for k = 1:length(PPsat)
    ph(k) = loglog(PP{k}, NN1_ss{k}, 'o-k', 'markerfacecolor','w');
    set(ph(k), 'marker', mark{k});

    if PPsat(k) <= 1
        th = text(PPsat(k), .25, sprintf('%.1f mW',PPsat(k)/1e-3));
    else
        th = text(PPsat(k), .25, sprintf('%.1f W',PPsat(k)));
    end
%     th.BackgroundColor = 'w';
    th.Position = [th.Extent(1) th.Extent(2)*.5];
%     th2 = text(PP{k}(1), .0001, legendStr{k});
%     th2.Rotation = 59;
%     th2.FontWeight = 'bold';
%     th2.BackgroundColor = 'w';
end
plot(PPsat, .25*ones(size(PPsat)), 'ob', 'markersize',5, 'markerfacecolor','b','linewidth',2)

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
legend(legendStr, 'location','southeast')

% lines at .5 and .25
a1 = plot(get(gca,'xlim'),0.5*[1 1], 'r--');
a2 = plot(get(gca,'xlim'),.25*[1 1], 'r--');
uistack([a1 a2], 'bottom')
hold off


% beautifyPlot

return

p = Plot(gcf);
% p.BoxDim = [12 6];
set(gca,'FontSize',12);
% set(gca,'FontWeight','bold');
p.MarkerSpacing = 1;
p.Colors = {'b'};
p.YTick = sort([logspace(-5,0,6), .25, .5]);
% p.Markers = [repelem({'none'}, 4), 's'];
p.Markers = mark;
p.LineStyle{end} = 'none';
p.Colors{end} = 'r';


for k = 1:length(PPsat)
    if PPsat(k) <= 1
        th = text(PPsat(k), .25, sprintf('%.1f mW',PPsat(k)/1e-3));
    else
        th = text(PPsat(k), .25, sprintf('%.1f W',PPsat(k)));
    end
%     th.BackgroundColor = 'w';
    th.FontSize = 12;
    th.Position = [th.Extent(1) th.Extent(2)*1];
%     th2 = text(PP{k}(1)*.4, 2e-5, legendStr{k});
%     th2.Rotation = 58;
%     th2.BackgroundColor = 'w';
%     th2.FontSize = 20;
end

p.Legend = legendStr;
p.LegendLoc = 'northoutside';
p.LegendBox = 'northoutside';
% legend(legendStr, 'location','southeast')


return
