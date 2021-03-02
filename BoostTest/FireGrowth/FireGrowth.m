function FireSpread()
clc
clear
format long


%   1 Non-Fire      
%   2 Sustained     b
%   3 Vigorous      g
%   4 Interactive   k
%   5 Remote        m
%   6 Full room     c


Lifetime = 45;
N=20000;
dt=Lifetime/N;

StateCount = 6;

fpath = 'FireGrowth_Markov.dat';
X = load( fpath  );
clr = {'r', 'b', 'g', 'k', 'm', 'c','b','k'};

T = X(:, 1);
for i = 3:7
    plot(T, X(:,i), clr{i-1});
    hold on;
end;

MaxProbabiltiies(X);

axis([0 T(end) 10^-8 10^0]);
%grid on;
set(gca,'YMinorGrid', 'off');
set(gca, 'YScale', 'log');

FireSeverity(X);

function MaxProbabiltiies(X)
for i = 1:6
    [y(i), v(i)] = max(X(:,i+1));
    m(i) = X(v(i),1);
end;

%plot([m(6) m(6)], [10^-8 y(6)],  'k' );
hold on;
[m; y]'


function FireSeverity(X)

figure;
plot(X(:,1),1 - X(:,2)- X(:,7));
set(gca, 'YScale', 'log');

%timespent = trapz( X(:,1), X(:,7)) 



function FlashoverPDF(X)
figure;

dt = diff(X(:,1));
dF = diff(X(:,2));

f = dF ./ dt;

plot( X(2:end,1), f);
hold on
plot( 0:45, exppdf(0:45,1/0.5), 'r')
axis([0 20 0 0.5]);