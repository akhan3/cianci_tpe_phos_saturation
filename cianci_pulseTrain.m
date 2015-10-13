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
% P = 1E-3;
lambda = 780e-9; % 800 nm
f = 80e6;       % Hz
alfa = 100e-15; % s


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
N1_ss_prev = 0/0;


%% Loop
while (true)
    fprintf(1, repmat('\b',1,charCount)); % delete line before
    charCount = fprintf('\t iC = %d', iC); 

    for ir = 1:length(r)
        for iz = 1:length(z)
            if iC == 1
                N1_initial = 1e-10;
            else
                N1_initial = N1(end,iz,ir,iC-1);
            end
%             [~, t(:,iz,ir,iC), N1(:,iz,ir,iC), gp(:,iC)] = cianci_model(P, lambda, f, alfa, Sr(iz,ir), tpa, gamma, N1_initial, true);
            [t_ss_exc(iz,ir,iC), N1_ss_exc(iz,ir,iC), t(:,iz,ir,iC), N1(:,iz,ir,iC), gp(:,iC)] = cianci_model(P, lambda, f, alfa, Sr(iz,ir), tpa, gamma, N1_initial, excitationType, true);
            t(:,iz,ir,iC) = t(:,iz,ir,iC) + (iC-1)/f; 
            if (~quietStatus)
                fprintf('N1_ss_exc(%d,%d,%d) = %G\n',             ir,iz,iC, N1_ss_exc(iz,ir,iC)); 
    %             fprintf('N1   (%d,%d,%d,%d) = %G\n',   size(t,1),ir,iz,iC, N1(end,iz,ir,iC)); 
            end
        end
    end
    
    if ~quietStatus && false
%         cla
        ph = plot(  squeeze(t(:)), .6*squeeze(gp(:)),'-r',...
                        squeeze(t(:)), squeeze(N1(:)),'-b',...
                        squeeze(t_ss_exc(end)), squeeze(N1_ss_exc(end)),'ok'); 
        grid on
%         set(ph, 'linewidth', 2);
%         ylim([-.05 .65])
        str = sprintf('P = %.2f mW', P*1e3); title(str);
%         myplot;
        drawnow;
%         keyboard
    end
    
%     N1_ss = N1(end,z==0,r==0,iC);
    N1_ss = N1_ss_exc(z==0,r==0,iC);
    
    if sum(isnan(N1_ss ) | isinf(N1_ss ) | N1_ss < 0)
        fprintf('\n'); 
        keyboard;
    end
    
    if abs(1-N1_ss/N1_ss_prev) < 0.000136487992059  % 1e-3
        break
    elseif N1_ss == N1_ss_prev
        break
    elseif isnan(N1_ss)
        fprintf('ERROR: Power too high or bad parameters\n')
        break
    end
    N1_ss_prev = N1_ss;

    iC = iC+1;
    
%     if iC == 10
%         break;
%     end
    
end

% fprintf('\n');

fprintf(1, repmat('\b',1,charCount)); % delete line before
%     charCount = fprintf('\tiC = %d', iC);

% fprintf('%d\n', iC);
% [N1_ss_prev, N1_ss, abs(1-N1_ss/N1_ss_prev),   abs(1-N1_ss/N1_ss_prev) < 0.000136487992059]
N1_00 = squeeze(N1(:,z==0,r==0,:));

% toc;

%% Plotting
% cla
% ph = plot(  squeeze(t(:)), .6*squeeze(gp(:)),'-r',...
% ph = plot(  ...  %squeeze(t(:)), .5*squeeze(gp(:)),'-',...
%                 squeeze(t(:)), squeeze(N1(:)),'-b'); grid on

ph = plot(  squeeze(t(:)), .6*squeeze(gp(:)),'-r',...
                squeeze(t(:)), squeeze(N1(:)),'-b',...
                t_ss_exc(:,:,end) + (iC-1)/f, N1_ss_exc(:,:,end),'ob'); 
set(ph(2), 'markerfacecolor', 'b');
grid on

% set(gca,'yscale','log')            
ylim([-.05 .65])
% set(ph, 'linewidth', 2);
% set(gca, 'xtickLabel',    get(gca, 'xtick') / 1e-6);
xlabel('Time (s)');
str = sprintf('P = %.2f mW (%s)', P*1e3, excitationType); title(str);
drawnow;
% pause(.5)

t0 = alfa*10;
tCondensed = sort(unique([t0, [t0-alfa*3 : alfa/12 : t0+alfa*3]]));

x = squeeze(t(1,:,:,:)) + tCondensed(end); 
y = squeeze(N1_ss_exc(:));
W = tpa * ( Sr * (P / (h*c/lambda)) )^2;
timeConst = 1 / (2*W + gamma);
% [P   W    (2*W + gamma)  5*timeConst    t(end,end,end,end)    N1(end,end,end,end)      W/(2*W+gamma)    N1_ss_exc(1)]'
% save('xy',   'x','y');

return   



x = squeeze(t(:));
y = squeeze(N1(:));
yy = max(y)-y;

% title([Ntime, N1_00])
% keyboard 
return