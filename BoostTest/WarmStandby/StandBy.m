function StandBy()

clc;
clear;

%TestQuantities();


%figure;
%MarkovPlot();
SystemPlot('ExpExp');
SystemFromComponents('ExpExp');

function SystemPlot(dist)

    Z = load(cell2mat({dist 'StandbySystem.dat'}));
    t = 1:length(Z);
       
    T = Z(t,1);
    X = Z(t,2);
        
    clrs = {'r', 'r', 'k', 'g'};

    scaling = 1; max(X);
    plot(T, X/scaling, clrs{1});  
    hold on;
    grid on;
    title('system plot');
    smp=X(end)
    
function SystemFromComponents(dist)
    cv = 1.0:-0.1:0.5;
    
    for i = cv
figure;
        Z = load(cell2mat({dist 'StandbyComponent' num2str(i) '.dat'}));
        t = 1:350:length(Z);

        T = Z(t,1);
        X = Z(t,2);

        clrs = {'r', 'r', 'k-', 'g'};

        p = X.*X;
        scaling = 1; max(p);
        plot(T, p/scaling, clrs{3});  
        hold on;
        grid on;
        title('system plot');
        probrule=X(end).*X(end)
    end;
function AnalyticalTest()
mu = 1/0.3;
lamda = 1/30;

(  mu*(lamda+mu))/( lamda*(lamda+mu) + mu*mu);
        

function MarkovPlot()
    Lifespan = 0:0.01:90;
    aSize =  3;
    p0 = zeros(1,aSize)';
    p0(1)=1;
    dpdt = ode5(@mssode_Standby, Lifespan , p0);
    T = 1:250:length(Lifespan);
    plot( Lifespan(T), dpdt(T, 3),'b>-'  );
    xlabel('time');
    ylabel('System Unavailability');
    hold on;
    markov=dpdt(end,3)
    
    grid on;
    
function    TestQuantities()
    a1f=30;
    a1r=0.3;
    
    g = @(t) exppdf(t, a1r);
    G = @(t) expcdf(t, a1r);
    f = @(t) exppdf(t, a1f);
    F = @(t) expcdf(t, a1f);
    R = @(t) 1-F(t);
    
    t = 0.216;
    q01 = f(t);
    w0 = R(t);
    q10 = g(t)*R(t);
    q11 = g(t)*F(t);
    w1 =  (1-G(t))*R(t);
    
    [0 q01; q10 q11];
    [w0 w1];
    

function lognornmalParams()
    sd = sqrt(0.6);
    mean = 2.4;
    lambda = 1/mean;
    
    
        
    
    
    