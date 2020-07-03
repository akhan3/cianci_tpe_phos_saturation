function [fluxDensity] = power2FluxDensity(P, lambda, beamWaist)
    %% physical constants
    h = 6.63e-34; % J.s
    c = 3e8; % m/s
    
    fluxDensity = P / (h*c/lambda) * gaussianBeamPSF(lambda, beamWaist);
end

function [Sr] = gaussianBeamPSF(lambda, beamWaist)
   %% Gaussian beam PSF
    w0 = beamWaist;
    z_rayleigh = pi * w0^2 / lambda;
    r = 0;    z = 0;
    [R,Z] = meshgrid(r,z);
    w = w0 * sqrt(1 + (Z./z_rayleigh).^2);
    S_gaussian = (2/pi)./(w.^2) .* exp(-2 * (R./w).^2);
    Sr = S_gaussian;
    assert(numel(Sr) == 1);
end