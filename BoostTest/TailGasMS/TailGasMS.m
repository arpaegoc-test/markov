%                       #    Sym   Col#     Life
%Booster fan            1    A      2       5
%Quench pumps           2    B,C    3       10
%Feedwater pump         1    D      4       10
%Circulation pumps      2    E,F    5       10
%Filter                 1    G      6       5

%this data set was generated with 15000 steps on the interval line
%with mission time of 20 years
%For 10 years, the probabilities appear at 7501th line of the file
%this can be used to test the results at just one point


function TailGas()

clc;
clear;


Z = load('TailGasMS.dat');
%       A B C D E F G  
cols = [2 3 3 4 5 5 6];
t = 1:1:length(Z);
T = Z(t,1);

A = 1-Z(t,2);
B = 1-Z(t,3);
C = 1-Z(t,3);
D = 1-Z(t,4);
E = 1-Z(t,5);
F = 1-Z(t,5);
G = 1-Z(t,6);

PlotComponents(T,A,B,C,D,E,F,G);

PlotSystemHumanMethod(T,A,B,C,D,E,F,G);
PlotSystemAutomatedMethod(T,A,B,C,D,E,F,G);




function PlotComponents(T,A,B,C,D,E,F,G)
figure;
hold on;
plot(T, 1-A, 'b');
plot(T, 1-D, 'r');
plot(T, 1-G, 'b-.'); 
set(gca, 'YScale', 'log')
PlotComponentsMCADG();

figure;
hold on;
plot(T, 1-B, 'b');
plot(T, 1-E,'r');
set(gca, 'YScale', 'log')
PlotComponentsMCBCEF();

function PlotSystemHumanMethod(T,A,B,C,D,E,F,G)
figure;
hold on;
p = A.* ( B + C - B.*C).*D.*(E + F - E.*F).*G;
plot(T,1-p);
grid on;
set(gca, 'YScale', 'log')

function PlotSystemAutomatedMethod(T,A,B,C,D,E,F,G)

% E1 = (A .* B .* D .* E .* G);
% E2 = (A .* B .* D .* F .* G);
% E3 = (A .* C .* D .* E .* G);
% E4 = (A .* C .* D .* F .* G);

cases=15;

UM = [
1	1	0	1	1	0	1
1	1	0	1	0	1	1
1	0	1	1	1	0	1
1	0	1	1	0	1	1
1	1	0	1	1	1	1
1	1	1	1	1	0	1
1	1	1	1	1	1	1
1	1	1	1	1	1	1
1	1	1	1	0	1	1
1	0	1	1	1	1	1
1	1	1	1	1	1	1
1	1	1	1	1	1	1
1	1	1	1	1	1	1
1	1	1	1	1	1	1
1	1	1	1	1	1	1   ];

weights = [1	1	1	1	-1	-1	-1	-1	-1	-1	1	1	1	1	-1];

for i = 1:length(T)
    umExpanded = repmat([A(i) B(i) C(i) D(i) E(i) F(i) G(i)],cases,1).*UM;
    x = find( umExpanded == 0);
    umExpanded(x) = 1;
    umProd = prod(umExpanded,2);
    weightedUM = umProd .* weights';
    unavail(i) = 1-sum(weightedUM);
end;
plot(T, unavail,'r-.');
idx = 7501;  % at 10 years
grid on;
set(gca, 'YScale', 'log')

function AtTime10Years(T,A1,B1,C1,D1,E1,F1,G1)
idx = 7501;  %at time 10 years

A = A1(idx); 
B = B1(idx);
C = C1(idx);
D = D1(idx);
E = E1(idx);
F = F1(idx);
G = G1(idx);

(1-[A B C D E F G])'   ;     %UNAVAILABILITY
[A B C D E F G]'        %AVAILABILITY

p = A.* ( B + C - B.*C).*D.*(E + F - E.*F).*G;
1-p



function PlotComponentsMCADG()

Z = load('TailGasMSMC.dat');
%       T A B C D E F G  
cols = [1 2 3 3 4 5 5 6];
idx =  [2 4 6]; %ADG

t = 1:1:length(Z);
T = Z(t,1);

for i = 1:length(idx)
    plot(T, Z(:,idx(i)),'o');
end;


function PlotComponentsMCBCEF()

Z = load('TailGasMSMC.dat');
%       T A B C D E F G  
cols = [1 2 3 3 4 5 5 6];
idx =  [3 5]; %ADG

t = 1:1:length(Z);
T = Z(t,1);

for i = 1:length(idx)
    plot(T, Z(:,idx(i)),'o');
end;