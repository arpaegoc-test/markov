function FireSpread()
clc
clear
format long
cv = 0.4;
figure;

Lifetime = 30;
N=20000;
dt=Lifetime/N;

StateCount = 8;

p0 = zeros(1,StateCount)';
p0(1)=1;
p = ode5(@mssode_fire, 0:dt:Lifetime, p0);
plot(0:dt:Lifetime, p(:,2) ,'b' );

grid on;
set(gca,'YMinorGrid', 'off');
hold on;

