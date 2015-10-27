clear
% clc
% close all

figure('windowStyle','docked');
ax1 = subplot(121); hold on; myplot; axis square;
ax2 = subplot(122); hold on; myplot; axis square;
drawnow;

% h = pauseButton;
% pause(0.01); % To create the button

verbosity = 1;

excitationType = 'Sech2'; 
% excitationType = 'Gaussian'; 
% excitationType = 'Rect'; 
% excitationType = 'CW'; 

%% Excitation source
beamWaist = 0.35e-6;    % Gaussian beam waist [m]
lambda = 780e-9;        % 800 nm
f = 80e6;               % 80 MHz
fwhm = 100e-15;         % 100 fs

switch excitationType
    case 'CW'
        fwhm = 1/f; % s
        fact = 1;
    case 'Gaussian'
        %fact = sqrt(pi); % * (2*sqrt(log(2)));
        fact = 2*sqrt(log(2));
    case 'Sech2'
        %fact = 2; % * (2*acosh(sqrt(2)));
        fact = 2*acosh(sqrt(2));
    case 'Rect'
        fact = 1;
end

% TPA_GM = 200;
TPA_GM = round(logspace(log10(1),log10(200), 4))' ;

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

    if strcmp(excitationType, 'CW')
        TAU = 1e-6 * [.001, .01, .1, 1, 10, 100]' ; 
%         TAU = 1e-6 * [.001, .01, .1, 1, 10]' ; 
    else
        % TAU = 1e-6 * [.001, .01, .1, 1, 10, 100]' ; 
        TAU = 1e-6 * [           .1, 1, 10, 100]' ;
        TAU = 1e-6 * [           .1, 1, 10]' ;
        TAU = 1e-6 * logspace(-1,1,9)' ;
    end
    
%     TAU = 1e-6 * [.1, 1, 10]';

    TAU = sort(unique(TAU));
%     TAU = sort(TAU,'descend');



    for idx_TAU = 1:length(TAU)
        tau = TAU(idx_TAU);
        if verbosity >= 1
            fprintf('%d/%d: TPA = %g GM:\t%d/%d: TAU = %s:\n', idx_TPA, length(TPA_GM), TPA_GM(idx_TPA),  idx_TAU, length(TAU), tauStr(TAU(idx_TAU)));
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
            opts.StartPoint = median(xData);
            opts.Upper = inf;
            [fr, gof] = fit( xData, yData, ft, opts );
            Ps_fit = fr.x0;

            Ps = Ps_fit;
            Psat_gof(idx_TAU,1) = gof.rsquare;

        elseif attempt == 2
            % interpolation only
            g1 = find(unique(yData));
            P_ = xData(g1);
            N_ = N1_ss(g1);
            Ps_interp1 = interp1(N_,P_,.25);
            Ps = Ps_interp1;
            Psat_gof(idx_TAU,1) = nan;
        end
        Psat(idx_TAU,1) = Ps;


        %% plotting N1-P saturation curve
        axes(ax2); hold on;
        subplot(122)
        myplot
        set(gca,'xscale','log');
        set(gca,'yscale','log');
        grid on; 
        hold on;
        xx = logspace(log10(min(P)), log10(max(P)),100);
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
        plot(Ps, .25,'sk', 'markerfacecolor', 'w');
        
        % model overlay
        TAU0_MODEL = 2.2586e-8;
        xPhi = power2FluxDensity(xx, lambda, beamWaist);
        xPhiSat = 1./sqrt(2*tpa) .* 1./sqrt(TAU0_MODEL + tau) .* sqrt(f*fwhm*fact);
        yN1 = (1/2) ./ (1 + (xPhiSat./xPhi).^2);
        plot(xx, yN1, '--g');
        [yy, yN1] = prepareCurveData(yy, yN1);
        rmse = rms(1 - yy(:) ./ yN1(:));
        hold off;

        if verbosity >= 1
            fprintf('\t0.5/(1+(%g/P)^2), gof.rsquare=%g\n', Ps, Psat_gof(idx_TAU,1));
            fprintf('\tRMSE = %.2f%%\n', 100*rmse);
        end


    % set(gca,'xlim',[1e-6 1e-0]);
    % if strcmp(excitationType, 'CW')
    %     set(gca,'xlim',[1e-3 1e3]);
    %     set(gca,'xtick',[1e-3, 1e0, 1e3]);
    % else
    %     set(gca,'ylim',[1e-3 1e-0]);
    % end
    % hold on
    % plot(get(gca,'xlim'), .25*ones(size(xlim)), '--k');
    % plot(get(gca,'xlim'), .50*ones(size(xlim)), '--k');
    % hold off

        xlabel('Average power, <P> (W)')
        ylabel('Excited state population, N_1');   % @ (\rho=0,z=0)')
        str = sprintf('%s, %G GM, [%s, %s]', excitationType, tpa/1e-58, tauStr(min(TAU)), tauStr(max(TAU)));
        title(str);

        axis square
        myplot
        drawnow


        PPP{idx_TAU,idx_TPA} = P;
        NN1{idx_TAU,idx_TPA} = N1_ss;
        TTAU(idx_TAU,idx_TPA) = TAU(idx_TAU);
        TTPA(idx_TAU,idx_TPA) = tpa;
    end % for idx_TAU


    %% Seed the main curve fitting method with the power model % y = a*x^b
    [xData, yData] = prepareCurveData(TAU, Psat);
    ft = fittype( 'power1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.TolFun = 1e-10;
    opts.TolX = 1e-10;
    opts.Display = 'final';
    opts.Lower = [0 -.5];
    opts.StartPoint = [1 -.5];
    opts.Upper = [inf -.5];
    [fr, gf] = fit( xData, yData, ft, opts );
    A_start = fr.a;

    A = A_start;
    tau0 = 0;
    a = A^2*tpa;

    if 1
        %% Curve fitting for lowpass filter (Bode plot):             A
           %                                             y =  ----------------
           %                                                   sqrt(x + x0)
        [xData, yData] = prepareCurveData(TAU, Psat);
        ft = fittype( 'A/sqrt(x+x0)', 'independent', 'x', 'dependent', 'y' );
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
    end


    fprintf('A/sqrt(T0+T) = %g/sqrt(%g + T), gof.rsquare=%g\n', ...
                A, tau0, gf.rsquare);

    %% Book-keeping
    PSAT_TPA(:,idx_TPA) = Psat;
    AA_TPA(1,idx_TPA) = A;
    TTAU0_TPA(1,idx_TPA) = tau0;
            

    %% new figure    
    axes(ax1); hold on; %hold off; cla;
    % clf

    %% Plotting tau-Psat
    % cla; hold off;
    xx2 = logspace(log10(min(TAU)), log10(max(TAU)), 100)';
    loglog(xx2, feval(fr,xx2),'r-');
    hold on;
    ph = loglog(TAU,Psat,'ob');
    yP = fluxDensity2Power( sqrt( 1 ./ (2 * tpa .* (TAU0_MODEL + xx2)) ) .* sqrt(f*fwhm*fact), lambda, beamWaist);
    plot(xx2, yP, '--g');
    hold off;
    set(ph,'MarkerFaceColor','b');
    xlabel('Probe Lifetime, \tau (s)')
    ylabel('Saturation onset, P_{sat} (W)')
    if tau0 == 0
        str = sprintf('$\\frac{\\sqrt{%G}}{\\sqrt{%G}\\sqrt{\\tau}}$', a, tpa);
    else
        str = sprintf('$\\frac{%G/\\sqrt{%G}}{\\sqrt{\\tau+%G}}$', A*sqrt(tpa), tpa, tau0);
    end
    lh = legend({str});
    set(lh,'Interpreter','latex','location','southwest','fontsize',18,'box','off');
    str = sprintf('%s, %G GM', excitationType, tpa/1e-58);
    title(str)
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    axis square

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

%% Conversion to Flux
PHISAT_TPA = power2FluxDensity(PSAT_TPA, lambda, beamWaist);
% xPhiSat = 1./sqrt(2*tpa/gamma) .* sqrt(f*fwhm*fact)
PPP2 = fluxDensity2Power( sqrt( 1 ./ (2 * TTPA .* TTAU) ) .* sqrt(f*fwhm*fact), lambda, beamWaist);
PPP2 ./ PSAT_TPA ;

% [f A]



