function [t_ss,N1_ss,t,N1,pulse] = cianci_model(P, lambda, f, fwhm, Sr, tpa, gamma, N1_0, excitationType)

    t0 = 0 + fwhm*10;
    timeFinal = 1/f;
    dt1 = (1/f) / 10;
    dt2 = (1/gamma) / 10;
    dt = min([dt1, dt2*inf]);
    t = sort(unique([[0:dt:timeFinal], timeFinal]));
    
    tCondensed = sort(unique([0, t0, 2*t0, [t0-fwhm*2 : fwhm/2  : t0+fwhm*2], t0+fwhm*10]));

% tCondensed = sort(unique([0, t0, 2*t0, t0+fwhm*10]));
% t = [0,timeFinal];

    t = sort(unique([t, tCondensed]))';



    %% Solution
    [N1,pulse] = solDiffEq(N1_0, t,t0,f,fwhm,P,lambda,Sr,tpa,gamma, excitationType);
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
    y(x == 0) = 1/2;
    y(x > 0) = 1;
end


%% solution of differential equation 
function [N1,pulse] = solDiffEq(N1_0, t,t0,f,fwhm,P,lambda,Sr,tpa,gamma, excitationType)
    SE = true;  % stimulated emission
    assert(t(1) == 0);
    assert(SE == true);

    dN_0 = (1-2*N1_0);
    
    %% physical constants
    h = 6.63e-34; % J.s
    c = 3e8; % m/s

    switch excitationType 
        case 'CW'
            pulse = ones(size(t));
            phi = P/(h*c/lambda) * Sr * (pulse);
            W = tpa * phi.^2;
            % M = cumtrapz(t, (1+SE)*W+gamma); % numerical
            M = t .* ((1+SE)*W+gamma); % analytical
        case 'Gaussian'
            alfa = fwhm / (2*sqrt(log(2))); % alfa = fwhm/1.665
            pulse = exp(-((t-t0)/alfa).^2);
            phi = P/(h*c/lambda) * Sr * 1/f * (1/sqrt(pi) * 1/alfa * pulse);
            W = tpa * phi.^2;
            % M = cumtrapz(t, (1+SE)*W+gamma); % numerical
            W0 = tpa * (P/(h*c/lambda) * Sr * 1/f).^2;
            M = gamma*t + (1+SE)*W0 * (1/sqrt(8*pi) * 1/alfa) * ( erf(sqrt(2)*(t-t0)/alfa) + erf(sqrt(2)*t0/alfa) ) ; % analytical
        case 'Sech2'
            alfa = fwhm / (2*acosh(sqrt(2))); % alfa = fwhm/1.763
            pulse = sech((t-t0)/alfa).^2;
            phi = P/(h*c/lambda) * Sr * 1/f * (1/2 * 1/alfa * pulse);
            W = tpa * phi.^2;
            % M = cumtrapz(t, (1+SE)*W+gamma); % numerical
            W0 = tpa * (P/(h*c/lambda) * Sr * 1/f).^2;
            % YSA = 1/(6*alfa) + 1/(24*alfa) * (sech((t-t0)/alfa)).^3 .* (3*sinh((t-t0)/alfa) + sinh(3*(t-t0)/alfa));
            YSA = 1/(6*alfa)  *                              (tanh((t-t0)/alfa) -                       tanh(-t0/alfa)) ...
                + 1/(12*alfa) * ((sech((t-t0)/alfa)).^2 .* tanh((t-t0)/alfa) - (sech(t0/alfa)).^2 .* tanh(-t0/alfa));
            M = gamma*t + (1+SE)*W0 * YSA; % analytical
        case 'Rect'
            error('Analytical integration not implemented yet for RectPulse')
            alfa = fwhm;
            pulse = U(t-t0+alfa/2) - U(t-t0-alfa/2);
            phi = P/(h*c/lambda) * Sr * 1/f * (1/alfa * pulse);
            W = tpa * phi.^2;
            M = cumtrapz(t, (1+SE)*W+gamma); % numerical
    otherwise
        error('Invalid choice of excitationType');
    end

    %% Numerical integration
    dN = dN_0*exp(-M) +  exp(-M) .* cumtrapz( t, (gamma+(SE-1)*W).*exp(M) );

    %% Correct for infinity
    Q = find(isnan(dN) | isinf(dN));
    if ~isempty(Q)
        fprintf('[corr] ') 
        % keyboard
        for iq = 1:length(Q)
            k = Q(iq); 
            dN(k) = dN_0*exp(-M(k)) + trapz(t(1:k), (gamma+(SE-1)*W(1:k)) .* exp(M(1:k)-M(k)));
        end
    end
    
    %% convert dN to N1
    N1 = 1/2*(1-dN);    

    if sum(N1>.5)
        fprintf('\tN1 exceeds 0.5...')
        keyboard
    end
    
    %% Correct over and under-saturation
    N1(N1<0) = 0;
    N1(N1>0.5) = 0.5;


%     if sum(isnan(N1) | isinf(N1) | N1<0 | N1>0.5)
%         fprintf('\n') 
%         keyboard
%     end

    
end
