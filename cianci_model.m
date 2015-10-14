function [t_ss,N1_ss,t,N1,pulse] = cianci_model(P, lambda, f, fwhm, Sr, tpa, gamma, N1_0, excitationType, doPlotting)

    %% physical constants
    h = 6.63e-34; % J.s
    c = 3e8; % m/s

    %% simulation parameters
    timeInit = 0;
    t0 = timeInit + fwhm*10;
    %     t0 = 0;
    %     timeInit = t0 - alfa*10;
    timeFinal = t0 + 10*fwhm;
    timeFinal = 1/f;

%     numPoints = 30;
%     dt = (timeFinal-timeInit) / numPoints;

%     dt = (1/gamma) / 200;
    dt = (1/f) / 1000;

    t = [timeInit:dt:timeFinal];
    tCondensed = sort(unique([t0, [t0-fwhm*3 : fwhm/12 : t0+fwhm*3]]));
    t = sort(unique([t, tCondensed]))';
%     t = sort(unique([t, t0, [t0-alfa*3 : alfa/12 : t0+alfa*3]]));
    if ~doPlotting 
        t = [t(1), t(end)];
    end
%     length(t)


    %% Solution
    [N1,pulse] = solDiffEq(N1_0, t,t0,f,fwhm,P,lambda,Sr,tpa,gamma, excitationType);

%     N1_neg = N1(N1<0);
%     if ~isempty(N1_neg)
%         N1 = N1 - min(N1_neg);
%     end
    
% % % % % %     fprintf('\t# %d negative values %s\n', sum(N1<0), repmat(['#'],[1,sum(N1<0)]))
    t_ss = tCondensed(end);
    N1_ss = N1(t == tCondensed(end));


    %% plotting
    if doPlotting && true && false
%         figure('windowStyle','docked')
%         hold on
        plot(t,pulse, '-', t,N1, '-o');  
%         hold off
        grid on;
%         ylim([0 1]);
%         xlim([0 2.5e-7]);
        drawnow;
%         pause(.05)
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Rectangular pulse function
function y = U(x)
    y = zeros(size(x));
%     y(x < 0) = 0;
    y(x == 0) = 1/2;
    y(x > 0) = 1;
end


%% solution of differential equation 
function [N1,pulse] = solDiffEq(N1_0, time,t0,f,fwhm,P,lambda,Sr,tpa,gamma, excitationType)
    SE = true;  % stimulated emission
    assert(time(1) == 0);
    assert(SE == true);

    dN_0 = (1-2*N1_0);
    
    %% physical constants
    h = 6.63e-34; % J.s
    c = 3e8; % m/s

    switch excitationType 
        case 'CW'
            pulse = ones(size(time));
            phi = P/(h*c/lambda) * Sr * (pulse);
            W = tpa * phi.^2;
%             M = cumtrapz(time, (1+SE)*W+gamma); % numerical
            M = time .* ((1+SE)*W+gamma); % analytical
        case 'GaussianPulse'
            alfa = fwhm / (2*sqrt(log(2)));
            pulse = exp(-((time-t0)/alfa).^2);
            phi = P/(h*c/lambda) * Sr * 1/f * (1/sqrt(pi) * 1/alfa * pulse);
            W = tpa * phi.^2;
%             M = cumtrapz(time, (1+SE)*W+gamma); % numerical
            W0 = tpa * (P/(h*c/lambda) * Sr * 1/f).^2;
            M = gamma*time + (1+SE)*W0 * (1/sqrt(8*pi) * 1/alfa) * ( erf(sqrt(2)*(time-t0)/alfa) + erf(sqrt(2)*t0/alfa) ) ; % analytical
        case 'Sech2Pulse'
            alfa = fwhm / (2*acosh(sqrt(2)));
            pulse = sech((time-t0)/alfa).^2;
            phi = P/(h*c/lambda) * Sr * 1/f * (1/2 * 1/alfa * pulse);
            W = tpa * phi.^2;
%             M = cumtrapz(time, (1+SE)*W+gamma); % numerical
            W0 = tpa * (P/(h*c/lambda) * Sr * 1/f).^2;
%             YSA = 1/(6*alfa) + 1/(24*alfa) * (sech((time-t0)/alfa)).^3 .* (3*sinh((time-t0)/alfa) + sinh(3*(time-t0)/alfa));
            YSA = 1/(6*alfa) *                               (tanh((time-t0)/alfa) -                       tanh(-t0/alfa)) ...
                + 1/(12*alfa) * ((sech((time-t0)/alfa)).^2 .* tanh((time-t0)/alfa) - (sech(t0/alfa)).^2 .* tanh(-t0/alfa));
            M = gamma*time + (1+SE)*W0 * YSA; % analytical
        case 'RectPulse'
            alfa = fwhm;
            pulse = U(time-t0+alfa/2) - U(time-t0-alfa/2);
            phi = P/(h*c/lambda) * Sr * 1/f * (1/alfa * pulse);
            W = tpa * phi.^2;
            M = cumtrapz(time, (1+SE)*W+gamma); % numerical
    otherwise
        error('Invalid choice of excitationType');
    end

    %% Numerical integration
    dN = dN_0*exp(-M) +  exp(-M) .* cumtrapz( time, (gamma+(SE-1)*W).*exp(M) );

    %% Correct for infinity
    Q = find(isnan(dN) | isinf(dN));
    if ~isempty(Q)
        fprintf('\nCorrecting...\n') 
        for iq = 1:length(Q)
            k = Q(iq); 
            dN(k) = dN_0*exp(-M(k)) + trapz(time(1:k), (gamma+(SE-1)*W(1:k)) .* exp(M(1:k)-M(k)));
        end
    end
    
    %% convert dN to N1
    N1 = 1/2*(1-dN);    

    %% Correct over and under-saturation
    N1(N1<0) = 0;
    N1(N1>0.5) = 0.5;

%     if sum(isnan(N1) | isinf(N1) | N1<0 )
%         fprintf('\n') 
%         keyboard
%     end

    
end
