function [N1_ss,t,N1,gp] = cianci_model(P, lambda, f, alfa, Sr, tpa, gamma, N1_0, excitationType, doPlotting)

    %% physical constants
    h = 6.63e-34; % J.s
    c = 3e8; % m/s

    %% simulation parameters
    timeInit = 0;
    t0 = timeInit + alfa*10;
    %     t0 = 0;
    %     timeInit = t0 - alfa*10;
    timeFinal = t0 + 10*alfa;
    timeFinal = 1/f;

%     numPoints = 30;
%     dt = (timeFinal-timeInit) / numPoints;

%     dt = (1/gamma) / 200;
    dt = (1/f) / 1000;

    t = [timeInit:dt:timeFinal];
    tCondensed = sort(unique([t0, [t0-alfa*3 : alfa/12 : t0+alfa*3]]));
    t = sort(unique([t, tCondensed]))';
%     t = sort(unique([t, t0, [t0-alfa*3 : alfa/12 : t0+alfa*3]]));
    if ~doPlotting 
        t = [t(1), t(end)];
    end
%     length(t)


    %% Solution
    [N1,gp] = solDiffEq(N1_0, t,t0,f,alfa,P,lambda,Sr,tpa,gamma, excitationType);

%     N1_neg = N1(N1<0);
%     if ~isempty(N1_neg)
%         N1 = N1 - min(N1_neg);
%     end
    
% % % % % %     fprintf('\t# %d negative values %s\n', sum(N1<0), repmat(['#'],[1,sum(N1<0)]))
    N1_ss = N1(t == tCondensed(end));


    %% plotting
    if doPlotting && true && false
%         figure('windowStyle','docked')
%         hold on
        plot(t,gp, '-', t,N1, '-o');  
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

%% Gaussian pulsed function
function y = GaussianPulse(t,t0,alfa)
    m = sqrt(8 * log(2)) / alfa; 
    y = exp(-(1/2) * (m*(t-t0)).^2);
end


%% Rectangular pulse function
function y = U(x)
    y = zeros(size(x));
%     y(x < 0) = 0;
    y(x == 0) = 1/2;
    y(x > 0) = 1;
end


%% solution of differential equation 
function [N1,gp] = solDiffEq(N1_0, time,t0,f,alfa,P,lambda,Sr,tpa,gamma, excitationType)
%     ti = time(1);
    assert(time(1) == 0);

    dN_0 = (1-2*N1_0);
    
    %% physical constants
    h = 6.63e-34; % J.s
    c = 3e8; % m/s

    switch excitationType 
        case 'CW'
            beta = (P / (h*c/lambda)).^2;
            gp = ones(size(time));
            W = @(t) tpa * Sr^2 * beta .* ones(size(t));
            M = cumtrapz(time, 2*W(time)+gamma);
            % M = time .* (2*W(time)+gamma);
            dN = dN_0*exp(-M) +  gamma * exp(-M) .* cumtrapz( time, exp(M) );
            Q = find(isnan(dN) | isinf(dN));
            for iq = 1:length(Q)
                k = Q(iq); 
                dN(k) = dN_0*exp(-M(k)) + gamma * trapz(time(1:k), exp(M(1:k)-M(k)));
            end
        case 'GaussianPulse'
            m = sqrt(8 * log(2)) / alfa; 
%             beta = 1/2*(P / (h*c/lambda) * 1/f).^2 .* m ./ (2 * sqrt(pi));
            gp = GaussianPulse(time,t0,alfa);
            W0 = tpa * Sr^2 * (P/(h*c/lambda) * 1/(f*alfa)).^2  .* (m*alfa/(2*pi))^2;
            W = @(t) W0 .* GaussianPulse(t,t0,alfa).^2;
            M = cumtrapz(time, 2*W(time)+gamma);
            dN = dN_0*exp(-M) +  gamma * exp(-M) .* cumtrapz( time, exp(M) );
%             Ma = @(t)  gamma*t + 2*W0/m * ( erf(m*(t-t0)) + erf(m*t0) ) ;
%             dN = dN_0*exp(-Ma(time)) +  gamma * exp(-Ma(time)) .* cumtrapz( time, exp(Ma(time)) );
        case 'RectPulse'
            beta = 1/2* 2*(P / (h*c/lambda) * 1/(f*alfa)).^2;
            gp = U(time-t0+alfa/2) - U(time-t0-alfa/2);
            W = @(t) tpa * Sr^2 * beta .* gp.^2;
            M = cumtrapz(time, 2*W(time)+gamma);
            dN = dN_0*exp(-M) +  gamma * exp(-M) .* cumtrapz( time, exp(M) );
        case 'RectPulse_old'
            warning('Did you really choose rectangular pulse?')
            error('Did you really choose rectangular pulse? Comment this line')
            beta = 2 * (P / (h*c/lambda) * 1/(f*alfa)).^2;
            gp = U(time-t0+alfa/2) - U(time-t0-alfa/2);
            intF = @(t)  exp( gamma * t + (tpa * Sr^2 * beta) * ( (t-t0+alfa/2) .* (U(t-t0+alfa/2) - U(t-t0-alfa/2) ) + alfa     .*       U(t-t0-alfa/2)   )   );
    otherwise
        error('Invalid choice of excitationType');
    end

    N1 = 1/2*(1-dN);    % convert dN to N1
    
%     M = cumtrapz(time, 2*W(time)+gamma);
%     N1 = N1_0*exp(-M) +  exp(-M) .* cumtrapz( time, W(time) .* exp(M) );



    if isnan(N1) | isinf(N1) | N1<0
%         fprintf('\n') 
%         keyboard
    end

    
end
