clc;
clear;
format long;

global SpareCount;
r = 0.07;      % 7% interest rate
m = 40;
i = 1:m;
pv = (1./(1+r).^i)';
c = 8;       %cost of 1 transformer
maxSpares = 8;
C_inv = (0:maxSpares) .* c;

f = @(n) 1 + n/0.05;


for SpareCount =0:maxSpares

    nStates = SpareCount + 2;

    v0 = zeros(1, nStates);
    dpdt_system = ode5(@rwdode, 0:0.05:m, v0');
    
    vAcc = dpdt_system(f(0:m),nStates)
    V = diff(vAcc)
    L(SpareCount+1) = sum(V .* pv);   % present value of future losses

    %Plots for the reward model
    %plot(0:0.05:m, dpdt_system(:, nStates)); hold on;
    %set(gca, 'YScale', 'log');
end;

C = L + C_inv;          % net present cost
M=[C_inv' L' C']

%plot for the net present cost
%plot( 0:maxSpares, M(:, 3))



%dvdt = ode45( @rwdode45, 0:1:m, v0);   %now obsolete
