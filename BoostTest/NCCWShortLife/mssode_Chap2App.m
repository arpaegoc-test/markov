function dpdt = mssode(t, p)
global CommonCause;

if CommonCause == 1                 % Full CommonCause
    L = 7.57742E-06;
    Lc = 1.4688E-07;
elseif CommonCause == 0.5          % Sensitivity Case
    L = 0.002000284;
    Lc = 0.000002;
elseif CommonCause == 0        % No CommonCause for semi-Markov
    L = 7.7243e-6;
    Lc = 0;       
end;

M = 1/19.4;

Lh = 3.4243e-6;
C = Lc + Lh;

a = [0 0 0; L+C 0 M; C 2*L 0];
aSize = length(a);

for i = 1:aSize
    a(i,i) = -sum(a(i,:));
end;
% 

% for i = 1:aSize     % Horz
%     txt=' ';
%     for j = 1:aSize % Vert
%         if ( i~= j) 
%             txt = strcat( txt, num2str(a(j,i)), 'p_{', num2str(j), '}(t) + ');
%         end;
%     end;
%     txt = strcat( '\frac{dp_{', num2str(i), '}(t)}{dt} = ', txt, num2str(a(i,i)), 'p_{', num2str(i), '}(t)')
% end;

dpdt =  (p'*a)';
