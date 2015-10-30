function [t_ss,N1_ss,t,N1,pulse] = cianci_model(P, lambda, f, fwhm, Sr, tpa, gamma, N1_0, excitationType, verbosity)

    analytical = false;

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
    
    if analytical == true
        tCondensed = sort(unique([0, t0, 2*t0]));
        t = sort(unique([0,1/f,tCondensed]))' ;
    end

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
    assert(t(1) == 0);
    assert(strcmp(excitationType,'CW') || strcmp(excitationType,'Sech2'));

    dN_0 = (1-2*N1_0);
    
    %% physical constants
    h = 6.63e-34; % J.s
    c = 3e8; % m/s

    P_avg = P;

    switch excitationType 
        case 'CW'
            P_pek = P_avg;
            pulse = ones(size(t)); % integrates to [1/f] over 1 period
            pulseSqrdInt = @(t) t; % analytical
        case 'Gaussian'
            alfa = fwhm / (2*sqrt(log(2))); % alfa = fwhm/1.665
            P_pek = P_avg / (f*alfa) / sqrt(pi);
            pulse = exp(-((t-t0)/alfa).^2); % integrates to [sqrt(pi)*alfa] over 1 period
            pulseSqrdInt = @(t) sqrt(pi/8) * alfa * ( erf(sqrt(2)*(t-t0)/alfa) + ...
                                                          erf(sqrt(2)*(  t0)/alfa) ) ; % analytical
        case 'Sech2'
            alfa = fwhm / (2*acosh(sqrt(2))); % alfa = fwhm/1.763
            P_pek = P_avg / (f*alfa) / 2;
            pulse = sech((t-t0)/alfa).^2; % integrates to [2*alfa] over 1 period
            pulseSqrdInt = @(t) alfa * ( 2/3 * tanh((t-t0)/alfa) + 1/3 * sech((t-t0)/alfa).^2 .* tanh((t-t0)/alfa) + ...
                                             2/3 * tanh((  t0)/alfa) + 1/3 * sech((  t0)/alfa).^2 .* tanh((  t0)/alfa) ); % analytical
        case 'Rect'
            alfa = fwhm;
            P_pek = P_avg / (f*alfa);
            pulse = U(t-t0+alfa/2) - U(t-t0-alfa/2); % integrates to [alfa] over 1 period
            pulseSqrdInt = @(t) (t-t0+alfa/2) .* U(t-t0+alfa/2) - (t-t0-alfa/2) .* U(t-t0-alfa/2) - ...
                                    ( -t0+alfa/2) .* U( -t0+alfa/2) - ( -t0-alfa/2) .* U( -t0-alfa/2); % analytical
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
    M = @(t) gamma*t + ((1+1)*W_pek) * pulseSqrdInt(t); % analytical
    
    %% Numerical integration
    dN = dN_0*exp(-M(t)) +  gamma .* exp(-M(t)) .* cumtrapz(t, exp(M(t)));
 
%     dN = zeros(size(t));
%     dN(1) = dN_0 * exp(-M(t(1)));
%     for k = 2:length(dN)
%         dN(k) = dN_0*exp(-M(t(k))) +  gamma * trapz( t(1:k), exp(M(t(1:k))-M(t(k))) );
% %         mk = M(t(k));
% %         dN(k) = dN_0*exp(-M(t(k))) +  gamma * quadgk( @(t) exp(M(t)-mk), t(1), t(k) ) ; %* exp(-M(t(k)));
%     end
    
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
                dN(k) = dN_0*exp(-M(t(k))) + gamma * trapz(t(1:k), exp(M(t(1:k))-M(t(k))));
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
