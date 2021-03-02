function dpdt = mssode_fire(t, p)

nStates=8;
a = zeros(nStates, nStates);

a=[
0	0.2	0	0.2	0	0.1	0	0
0.2	0	0.2	0	0	0	0.1	0
0	0.2	0	0.5	0	0	0	0.1
0.2	0	0.5	0	0.667	0	0	0
0	0	0	0.333	0	0.2	0	0.2
0.067	0	0	0	0.2	0	0.08	0
0	0.067	0	0	0	0.08	0	0.2
0	0	0.067	0	0.2	0	0.2	0
];




aSize = length(a);

for i = 1:aSize
    a(i,i) = -sum(a(i,:));
end;


% for i = 1:aSize     % Horz
%     txt=' ';
%     for j = 1:aSize % Vert
%         if ( i~= j) 
%             txt = strcat( txt, num2str(a(j,i)), 'p', num2str(j), '(t) + ');
%         end;
%     end;
%     txt = strcat( 'dp', num2str(i), '(t)/dt = ', txt, num2str(a(i,i)), 'p', num2str(i), '(t)')
% end;

dpdt =  (p'*a)';
