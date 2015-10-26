function [t_ss,N1_ss,t,N1,pulse] = cianci_model(P, lambda, f, fwhm, Sr, tpa, gamma, N1_0, excitationType, verbosity)
    
    %% Define time vector
    timeFinal = 1/f;
    dt1 = (1/f) / 120;
    dt2 = (1/gamma) / 10;
    dt = min([dt1, dt2*inf]);
    t = sort(unique([[0:dt:timeFinal], timeFinal]));
    if strcmp(excitationType, 'CW') % correction for CW
        fwhm = 100e-15;
        t0 = 0 + fwhm*10;
        tCondensed = sort(unique([0, t0, 2*t0, t0+fwhm*10]));
    else
        t0 = 0 + fwhm*10;
        tCondensed = sort(unique([0, t0, 2*t0, [t0-fwhm*2 : fwhm/12  : t0+fwhm*2], t0+fwhm*10]));
    end
    t = sort(unique([t, tCondensed]))';

    %% Solution
    [N1,pulse] = solDiffEq(N1_0, t,t0,f,fwhm,P,lambda,Sr,tpa,gamma, excitationType, verbosity);
    t_ss = tCondensed(end);
    N1_ss = N1(t == t_ss);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Rectangular pulse function
function y = U(x)
    y = zeros(size(x));
    % y(x < 0) = 0;
    % y(x == 0) = 1/2;
    y(x >= 0) = 1;
end


%% solution of differential equation 
function [N1,pulse] = solDiffEq(N1_0, t,t0,f,fwhm,P,lambda,Sr,tpa,gamma, excitationType, verbosity)
    SE = true;  % stimulated emission
    assert(t(1) == 0);
    assert(SE == true);

    dN_0 = (1-2*N1_0);
    
    %% physical constants
    h = 6.63e-34; % J.s
    c = 3e8; % m/s

    P_avg = P;

    switch excitationType 
        case 'CW'
            P_pek = P_avg;
            pulse = ones(size(t)); % integrates to [1/f] over 1 period
            pulseSqrdIntFunc = @(t) t; % analytical
            % pulseSqrdInt = cumtrapz(t,pulse.^2); % analytical
        case 'Gaussian'
            alfa = fwhm / (2*sqrt(log(2))); % alfa = fwhm/1.665
            P_pek = P_avg / (f*alfa) / sqrt(pi);
            pulse = exp(-((t-t0)/alfa).^2); % integrates to [sqrt(pi)*alfa] over 1 period
            pulseSqrdIntFunc = @(t) sqrt(pi/8) * alfa * erf(sqrt(2)*(t-t0)/alfa); % analytical
            % pulseSqrdInt = sqrt(pi/8) * alfa * ( erf(sqrt(2)*(t-t0)/alfa) + erf(sqrt(2)*t0/alfa) ) ; % analytical
        case 'Sech2'
            alfa = fwhm / (2*acosh(sqrt(2))); % alfa = fwhm/1.763
            P_pek = P_avg / (f*alfa) / 2;
            pulse = sech((t-t0)/alfa).^2; % integrates to [2*alfa] over 1 period
            pulseSqrdIntFunc = @(t) alfa * ( 2/3 * tanh((t-t0)/alfa) + 1/3 * sech((t-t0)/alfa).^2 .* tanh((t-t0)/alfa) ); % analytical
            % pulseSqrdInt = alfa/3 * ( (2 + sech((t-t0)/alfa).^2) .* tanh((t-t0)/alfa) + ...
            %                           (2 + sech((  t0)/alfa).^2) .* tanh((  t0)/alfa) ); % analytical
        case 'Rect'
            alfa = fwhm;
            P_pek = P_avg / (f*alfa);
            pulse = U(t-t0+alfa/2) - U(t-t0-alfa/2); % integrates to [alfa] over 1 period
            pulseSqrdIntFunc = @(t) (t-t0+alfa/2) .* U(t-t0+alfa/2) - (t-t0-alfa/2) .* U(t-t0-alfa/2); % analytical
    otherwise
        error('Invalid choice of excitationType');
    end
    
    %% other variables
    phi_avg = P_avg / (h*c/lambda) * Sr;
    phi_pek = P_pek / (h*c/lambda) * Sr;
    W_avg = tpa * phi_avg.^2;
    W_pek = tpa * phi_pek.^2;
    phi = phi_pek * pulse;
    W = W_pek * pulse.^2;
    pulseSqrdInt = pulseSqrdIntFunc(t) - pulseSqrdIntFunc(0);
    M = gamma*t + ((1+SE) * W_pek) * pulseSqrdInt; % analytical
    % M = cumtrapz(t, (1+SE)*W+gamma); % numerical
    
    
    %% Numerical integration
    dN = dN_0*exp(-M) +  exp(-M) .* cumtrapz( t, (gamma+(SE-1)*W).*exp(M) );
    
    %% Correct for infinity
    while true
        Q = find(isnan(dN) | isinf(dN));
        if isempty(Q)
            break
        else
            if verbosity >= 3
                fprintf('[corr] '); 
            end
            % keyboard
            for iq = 1:length(Q)
                k = Q(iq); 
                dN(k) = dN_0*exp(-M(k)) + trapz(t(1:k), (gamma+(SE-1)*W(1:k)) .* exp(M(1:k)-M(k)));
            end
        end
    end
    
    %% convert dN to N1
    N1 = 1/2*(1-dN);    

    %% Correct over and under-saturation
    if sum(N1<0)
        if verbosity >= 3
            fprintf('[N1<0] ');
        end
        % keyboard
        N1(N1<0) = 0;
    end
    
    if sum(N1>.5)
        if verbosity >= 3
            fprintf('[N1>0.5] ');
        end
        keyboard
        N1(N1>0.5) = 0.5;
    end
    
end
