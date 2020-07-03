clear
% clc
% close all


figure('windowStyle','docked');
ax1 = subplot(121); hold on; myplot; axis square;
ax2 = subplot(122); hold on; myplot; axis square;
drawnow;

% h = pauseButton;
% pause(0.01); % To create the button

verbosity = 100;

%% Excitation source
beamWaist = 0.35e-6;    % Gaussian beam waist [m]
lambda = 780e-9;        % 800 nm
f = 80e6;               % 80 MHz
fwhm = 100e-15;         % 100 fs


excitationType = 'Sech2';
% excitationType = 'Gaussian';
% excitationType = 'Rect';
% excitationType = 'CW';

% MODEL_TAU0 = 2.2586e-8;
MODEL_K2 = 1;%[0.991059853842970];

switch excitationType
    case 'CW'
        fwhm = 1/f; % s
        fact = 1;
        MODEL_TAU0 = 0;
    case 'Gaussian'
        %fact = sqrt(pi); % * (2*sqrt(log(2)));
        fact = 2*sqrt(log(2));
        MODEL_TAU0 = [1.87957753681724e-08];
    case 'Sech2'
        %fact = 2; % * (2*acosh(sqrt(2)));
        fact = 2*acosh(sqrt(2));
        MODEL_TAU0 = [1.87957753681724e-08];
    case 'Rect'
        fact = 1;
        MODEL_TAU0 = [1.87957753681724e-08];
end

TPA_GM = round(logspace(log10(1),log10(200), 4))' ;
TPA_GM = 200;

for idx_TPA = 1:length(TPA_GM)
    
    if verbosity >= 0
        fprintf('====================== %d/%d: TPA = %g GM ======================\t\n', idx_TPA, length(TPA_GM), TPA_GM(idx_TPA) );
    end
    
    %     figure('windowStyle','docked');
    %     ax1 = subplot(121); hold on; myplot; axis square;
    %     ax2 = subplot(122); hold on; myplot; axis square;
    %     drawnow;
    
    
    %% Fluorophore
    tpaGM = TPA_GM(idx_TPA);
    tpa = tpaGM * 1e-58;  % GM = 1e-58 m^4 / (photon/s)
    clear tpaGM;
    
    TAU = 1e-6 * [0.001, 0.01, 0.1, 1, 10, 100/100]' ;
    TAU = 1e-6 * [0.001, 0.01];
    TAU = sort(unique(TAU));
    %     TAU = sort(TAU,'descend');
    
    
    
    for idx_TAU = 1:length(TAU)
        tau = TAU(idx_TAU);
        if verbosity >= 1
            fprintf('\t%d/%d: TAU = %s:\n', idx_TAU, length(TAU), tauStr(TAU(idx_TAU)));
        end
        
        for attempt = 1:2
            if attempt == 1
                numPowerPoints = 6;
                pf = 1;
            elseif attempt == 2
                numPowerPoints = 18;
                pf = 1;
            end
            
            if strcmp(excitationType, 'CW')
                P = 3.5336e-32 / sqrt(tpa) / sqrt(tau) * logspace(-1, 1, numPowerPoints)'; % CW
            else
                P = (1.35e-34/sqrt(tpa)) / sqrt(2.04e-8 + tau) * logspace(-1, 1, numPowerPoints)';
            end
            
            % P = logspace(-4,-1,numPowerPoints)' ;
            % P = sort(P,'descend');
            
            
            %% Pre-allocate variables
            N1_ss = zeros(size(P));
            
            %% Create figure window for time traces
            if verbosity >= 5
                if attempt == 1
                    figure('windowStyle','docked');
                elseif attempt == 2
                    clf;
                end
                hold on;
                myplot;
            end
            
            %% Sweep power values
            for idx_P = 1:length(P)
                if verbosity >= 2
                    fprintf('\t[%d/%d]\t',   idx_P,length(P));
                end
                [N1_ss(idx_P), lastSlope] = cianci_pulseTrain(P(idx_P), lambda, f, fwhm, 1/tau, tpa, beamWaist, excitationType, verbosity);
                if verbosity >= 2; fprintf('\n'); end;
                if lastSlope == 0 && attempt == 1
                    break;
                end
            end
            
            if lastSlope == 0 && attempt == 1
                continue;
            end
            break
        end % for attempt
        
        
        %% Curve-fitting to find saturation threshold (Psat)
        [xData, yData] = prepareCurveData(P, N1_ss);
        if attempt == 1
            if verbosity >= 1; fprintf('\t'); end;
            ft = fittype( '0.5/(1+(x0/x)^2)', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.TolFun = 1e-10;
            opts.TolX = 1e-10;
            if verbosity >=2
                opts.Display = 'final';
            end
            opts.Lower = 0;
            opts.StartPoint =  median(xData);
            opts.Upper = inf;
            [fr, gof, fitAlgo] = fit( xData, yData, ft, opts );
            Ps = fr.x0;
        elseif attempt == 2
            % interpolation only
            g1 = find(unique(yData));
            P_ = xData(g1);
            N_ = N1_ss(g1);
            Ps_interp1 = interp1(N_,P_,.25);
            Ps = Ps_interp1;
            fr = 'interp1';
            gof.rsquare = nan;
            fitAlgo = 'interp1';
        end
        fit1_Psat(idx_TAU,idx_TPA) = Ps;
        fit1_fr  {idx_TAU,idx_TPA} = fr;
        fit1_gof {idx_TAU,idx_TPA} = gof;
        fit1_algo{idx_TAU,idx_TPA} = fitAlgo;
        
        %
        % %% Curve-fitting to find model TAU0 for in Psat fitting
        %
        %     N1_m = modelN1ss(P, tpa, tau, 0, f, fwhm, fact, lambda, beamWaist);
        %     [P N1_ss N1_m]
        %     funcFit1 = @(K1,xt0,x)    K1 * modelN1ss(x, tpa, tau, xt0, f, fwhm, fact, lambda, beamWaist);
        %     ft = fittype( funcFit1, 'independent', 'x', 'dependent', 'y' );
        %     opts.Lower = [0 0];
        %     opts.StartPoint = [1 TAU0_MODEL];
        %     opts.Upper = [inf inf];
        %
        % [xData, yData] = prepareCurveData(P, N1_ss);
        % if verbosity >= 1; fprintf('\t'); end;
        % ft = fittype( '0.5/(1+(x0/x)^2)', 'independent', 'x', 'dependent', 'y' );
        % opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        % opts.TolFun = 1e-10;
        % opts.TolX = 1e-10;
        % if verbosity >=2
        %     opts.Display = 'final';
        % end
        % opts.Lower = 0;
        % opts.StartPoint =  median(xData);
        % opts.Upper = inf;
        % [fr, gof, fitAlgo] = fit( xData, yData, ft, opts );
        % Ps = fr.x0;
        
        
        
        %% plotting N1-P saturation curve
        axes(ax2); hold on;
        subplot(122)
        myplot
        set(gca,'xscale','log');
        set(gca,'yscale','log');
        grid on;
        hold on;
        xx = logspace(log10(min(P)), log10(max(P)),100)' ;
        if attempt == 1
            yy = feval(fr,xx);
        elseif attempt == 2
            yy = interp1(xData,yData,xx);
        end
        ph1 = plot(xx,yy,'-');
        ph2 = loglog(xData, yData, 'o');
        col = get(ph1,'color');
        set(ph2, 'markerfacecolor', col);
        set([ph1,ph2],'color',col);
        plot(Ps, 0.25,'sk', 'markerfacecolor', 'w');
        
        % model overlay
        %         xPhi = power2FluxDensity(xx, lambda, beamWaist);
        %         xPhiSat = 1./sqrt(2*tpa) .* 1./sqrt(TAU0_MODEL + tau) .* sqrt(f*fwhm*fact);
        %         yN1 = (1/2) ./ (1 + (xPhiSat./xPhi).^2);
        %         xPs = fluxDensity2Power(xPhiSat, lambda, beamWaist);
        
        % xPs = MODEL_K2/MODEL_K2 * modelPsat(tpa, tau, MODEL_TAU0, f, fwhm, fact, lambda, beamWaist);
        xPs = modelPsat(tpa, tau, MODEL_TAU0, f, fwhm, fact, lambda, beamWaist);
        yN1 = (1/2) ./ (1 + (xPs./xx).^2);
        
        plot(xx, yN1, '--r');
        [yy, yN1] = prepareCurveData(yy, yN1);
        rmse = rms(1 - yy(:) ./ yN1(:));
        hold off;
        
        if verbosity >= 1
            fprintf('\t0.5/(1+(%g/P)^2), gof.rsquare=%g', fit1_Psat(idx_TAU,idx_TPA), fit1_gof{idx_TAU,idx_TPA}.rsquare);
            fprintf('\tPs_model=%g, RMSE_model = %.2f%%\n', xPs, 100*rmse);
        end
        
        
        % set(gca,'xlim',[1e-6 1e-0]);
        % if strcmp(excitationType, 'CW')
        %     set(gca,'xlim',[1e-3 1e3]);
        %     set(gca,'xtick',[1e-3, 1e0, 1e3]);
        % else
        %     set(gca,'ylim',[1e-3 1e-0]);
        % end
        % hold on
        % plot(get(gca,'xlim'), 0.25*ones(size(xlim)), '--k');
        % plot(get(gca,'xlim'), 0.50*ones(size(xlim)), '--k');
        % hold off
        
        xlabel('Average power, <P> (W)')
        ylabel('Excited state population, N_1');   % @ (\rho=0,z=0)')
        str = sprintf('%s, %G GM, [%s, %s]', excitationType, tpa/1e-58, tauStr(min(TAU)), tauStr(max(TAU)));
        title(str);
        
        axis square
        myplot
        drawnow
        
%         return
        
        PPP{idx_TAU,idx_TPA} = P;
        NN1{idx_TAU,idx_TPA} = N1_ss;
        TTAU(idx_TAU,idx_TPA) = tau;
        TTPA(idx_TAU,idx_TPA) = tpa;
    end % for idx_TAU
    
    clear Ps fr gof fitAlgo
    
    return

    %     %% Seed the main curve fitting method with the power model % y = a*x^b
    %     [xData, yData] = prepareCurveData(TAU, fit1_Psat);
    %     ft = fittype( 'power1' );
    %     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    %     opts.TolFun = 1e-10;
    %     opts.TolX = 1e-10;
    %     opts.Display = 'final';
    %     opts.Lower = [0 -.5];
    %     opts.StartPoint = [1 -.5];
    %     opts.Upper = [inf -.5];
    %     [fr, gf] = fit( xData, yData, ft, opts );
    %     A_start = fr.a;
    %
    %     A = A_start;
    %     tau0 = 0;
    %     a = A^2*tpa;
    
    %     PSAT_TPA(:,idx_TPA) = fit1_Psat;
    %     AA_TPA(1,idx_TPA) = A;
    %     TTAU0_TPA(1,idx_TPA) = fit2_tau0;
    
    %% Curve fitting for lowpass filter (Bode plot):             A
    %                                             y =  ----------------
    %                                                   sqrt(x + x0)
    [xData, yData] = prepareCurveData(TAU, fit1_Psat(:,idx_TPA));
    funcFit2 = @(K2,x0,x)    K2 * modelPsat(tpa, x, x0, f, fwhm, fact, lambda, beamWaist);
    ft = fittype( funcFit2, 'independent', 'x', 'dependent', 'y' );
    %         ft = fittype( 'A/sqrt(x+x0)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.TolFun = 1e-10;
    opts.TolX = 1e-10;
    opts.Display = 'final';
    opts.Lower = [0 0];
    opts.StartPoint = [1 1];
    if strcmp(excitationType, 'CW')
        opts.Upper = [inf 0];   % fix x0 at 0 for CW
    else
        opts.Upper = [inf inf];
    end
    [fr,gf,fitAlgo] = fit( xData, yData, ft, opts);
    
    fit2_K2(idx_TPA) = fr.K2;
    fit2_tau0(idx_TPA) = fr.x0;
    fit2_fr{idx_TPA} = fr;
    fit2_gof{idx_TPA} = gf;
    fit2_algo{idx_TPA} = fitAlgo;
    clear fr gf fitAlgo;
    
    fprintf('Psat = K2 / sqrt(2*tpa) / sqrt(tau0+tau) * sqrt(f*fwhm*fact): (K2=%g, tau0=%g), fit2_gof.rsquare=%g\n', ...
        fit2_K2(idx_TPA), fit2_tau0(idx_TPA), fit2_gof{idx_TPA}.rsquare);
    
    
    %% new figure
    axes(ax1); hold on; %hold off; cla;
    % clf
    
    %% Plotting tau-Psat
    % cla; hold off;
    xx2 = logspace(log10(min(TAU)), log10(max(TAU)), 100)';
    loglog(xx2, feval(fit2_fr{idx_TPA}, xx2),'r-');
    hold on;
    ph = loglog(TAU,fit1_Psat,'ob');
    % yModelPsat = modelPsat(tpa, xx2, fit2_tau0(idx_TPA), f, fwhm, fact, lambda, beamWaist);
    yModelPsat = modelPsat(tpa, xx2, MODEL_TAU0, f, fwhm, fact, lambda, beamWaist);
    plot(xx2, yModelPsat, '--g');
    hold off;
    set(ph,'MarkerFaceColor','b');
    xlabel('Probe Lifetime, \tau (s)')
    ylabel('Saturation onset, P_{sat} (W)')
    %         str = sprintf('$\\frac{%G/\\sqrt{%G}}{\\sqrt{\\tau+%G}}$', A*sqrt(tpa), tpa, fit2_tau0);
    str = sprintf('$(K_2=%G, {\\tau}_0=%G)$', fit2_K2(idx_TPA), fit2_tau0(idx_TPA));
    lh = legend({str});
    set(lh,'Interpreter','latex','location','southwest','fontsize',14,'box','off');
    str = sprintf('%s, %G GM', excitationType, tpa/1e-58);
    title(str)
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    axis square
    
    %     yModelPsat = modelPsat(tpa, xx2, MODEL_TAU0, f, fwhm, fact, lambda, beamWaist);
    rmse_(1) = rms(1 - yModelPsat ./ feval(fit2_fr{idx_TPA}, xx2));
    rmse_(2) = rms(1 - yModelPsat ./ (fit2_K2(idx_TPA)*yModelPsat));
    rmse_
    
    
    % set(gca,'xlim',[1e-9 1e-4]);
    % if strcmp(excitationType, 'CW')
    %     set(gca,'ylim',[1e-2 1e2]);
    % else
    %     set(gca,'ylim',[1e-4 1e-1]);
    
    % end
    myplot
    
    if 0
        %% printing data for saving
        folderName = sprintf('FILES_%s', excitationType);
        fileName = sprintf('Data %s %g GM %s - %s', excitationType, tpa/1e-58, tauStr(min(TAU)), tauStr(max(TAU)))
        if exist(folderName, 'file') ~= 7
            mkdir(folderName);
        end
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
        
        print(gcf, fullfile(folderName, [fileName,'.pdf']), '-dpdf');
        
        save(fullfile(folderName, [fileName,'.mat']));
    end
    
    %     return
    
end %% for idx_TPA
