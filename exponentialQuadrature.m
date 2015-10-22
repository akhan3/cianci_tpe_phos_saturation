clc

tau = 1;
A = pi/pi;
fun = @(t) A * exp(-t/tau);
x = [0:tau/10:20*tau]';
y = A * exp(-x/tau);

AreaN_simps = [];
for k = 1:length(y)
    AreaN_simps(k) = simps(x(1:k),y(1:k));
end
AreaN = cumtrapz(x,y);
AreaA = A*tau * (1-exp(-x/tau));
AreaN = AreaN_simps;


clf;

subplot(211)
plot(x,y,'o-k', x,AreaN,'-r', x,AreaA,'--b');
lh = legend('e^{-x/\tau}', 'Numerical integration', 'Analytical integration');
set(lh,'location','southeast','fontsize',12);
myplot

subplot(212)
plot(x,y,'o-k', x,AreaN,'-r', x,AreaA,'--b');
lh = legend('e^{-x/\tau}', 'Numerical integration', 'Analytical integration');
set(lh,'location','southeast','fontsize',12);
ylim((A*tau) * [.999998 1.000002]);
myplot