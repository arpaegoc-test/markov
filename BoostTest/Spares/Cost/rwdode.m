function dvdt = rwdode(t, v)

global SpareCount;
NTransformers = 12;
nSpares = SpareCount;
nStates = nSpares+2;
fr = 0.03;              % failure rate
M = 4;

%----Test Lisnianski
% NTransformers = 1;
% fr = 1;
% M = 200;
%----


%%%%%%%%%%% DEFINE THE TRANSITION RATE MATRIX a

L = NTransformers * fr;
spareNo = nSpares+1;
a = zeros(nStates, nStates);

for i = 1:nStates-1
    a(i+1, i) = L;      % FAILURE TRANSITION
    a(i, i+1) = spareNo * M;    %REPAIR TRANSITION
    spareNo = spareNo - 1;
end;

s = sum(a,2);
a = a-diag(s);

%%%%%%%%%%% DEFINE THE REWARD MATRIX r
Cp = 1;     % $1 per kwh
L = 10^5;   % Nominal capacity of 500/230 kV transformer:  KW
ens = L*8760;   %Energy Not Supplied (ENS) for 10^5 KW in a year: KW 
c = Cp*ens / 10^6;     %Loss in millions of $ for ENS
r01 = 0.05;    %Repair cost $50,000

%--------Test Lisnianski
% Cp = 3;     % $1 per kwh
% L = 10^5;   % Nominal capacity of 500/230 kV transformer:  KW
% ens = L*8760;   %Energy Not Supplied (ENS) for 10^5 KW in a year: KW 
% c = Cp*ens;     %Loss in $ for ENS
% r01 = 50000;    %
%----------

r = zeros(nStates, nStates);
r(1,1) = c;      % Loss in millions of $ for ENS
for i = 2:nStates
    r(i, i-1) = 0;      % 0 FOR A FAILURE TRANSITION
    r(i-1, i) = r01;   % Repair cost in millions of $
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


aSize = length(a);

for i = 1:aSize
    u(i) = r(i,i);
    for j = 1:aSize
        if(i ~= j)
            u(i) = u(i) + a(i,j)*r(i,j);
        end;
    end;
end;
u = u';

dvdt = u + (a*v) ;
