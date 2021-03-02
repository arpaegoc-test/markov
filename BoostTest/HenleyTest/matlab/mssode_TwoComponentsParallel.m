function dpdt = mssode_TwoComponentsParallel(t, p)

c1f = 1/30;
c1r=  1/0.5;

c2f=    1/30;
c2r=    1/0.5;

a = [
        0   c2r c1r 0
        c2f 0   0   c1r
        c1f 0   0   c2r
        0   c1f c2f 0
    ];



% a = [
%         0   c1f c2f 0   0
%         c1r 0   0   c2f 0
%         c2r 0   0   0   c1f
%         0   0   c1r 0   0
%         0   c2r 0   0   0
%     ];

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
