function pulseShapes()
clc
% clear

alfa = 2;
g_FWHM = 2 * sqrt(log(2)) * alfa;
s_FWHM = 2 * acosh(sqrt(2)) * alfa;
% [alfa g_FWHM s_FWHM] / 2

t0 = alfa*10;
dt = alfa/50;
t = [0 : dt : 100*alfa]';
t = unique([t; t0+[alfa/2; -alfa/2; 0;  g_FWHM/2; -g_FWHM/2; s_FWHM/2; -s_FWHM/2; ]]);

yg = 1/alfa * 1/sqrt(pi) * exp(-((t-t0)/alfa).^2);
ys = 1/alfa * 1/2 * (sech((t-t0)/alfa)).^2;
ys = ys.^2;
yp = 1/alfa * (U(t-t0+alfa/2)-U(t-t0-alfa/2));



YG = cumtrapz(t,yg);
YS = cumtrapz(t,ys);
YSA = 1/(6*alfa) + 1/(24*alfa) * (sech((t-t0)/alfa)).^3 .* (3*sinh((t-t0)/alfa) + sinh(3*(t-t0)/alfa));
% YSA = 1/(6*alfa) * (tanh((t-t0)/alfa) + tanh(t0/alfa)) + ...
%       1/(12*alfa) * ((sech((t-t0)/alfa)).^2 .* tanh((t-t0)/alfa) + (sech(t0/alfa)).^2 .* tanh(t0/alfa));
YP = cumtrapz(t,yp);


% figure('windowStyle','docked')
clf

% ph = plot(t,yg,'b-', t,ys,'r-', t,yp,'k-');
ph = plot(t,ys,'r-', t,YS,'bo', t,YSA,'g-');
legend({'sech^4(t)' , 'cumtrapz', 'Analytical'}, 'location','northoutside', 'orientation','horizontal')
% hold on;
% plot(t,YG,'b-.', t,YS,'k-', t,YP,'k-.', t,YSA,'m--')
% legend({'gaussian(t)' , 'sech^2(t)'},'location','northwest')
set(ph, 'markerfacecolor', 'w')
myplot

return

count = 0;
for n=1:100
    fprintf(1, repmat('\b',1,count)); %delete line before
    count = fprintf('| <-add line here: %3.2f*100',n/100);
    pause(.1)
end

end

function y = U(x)
    y = zeros(size(x));
%     y(x < 0) = 0;
    y(x == 0) = 1/2;
    y(x > 0) = 1;
end
