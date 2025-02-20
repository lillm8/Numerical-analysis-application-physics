u_fram = true;
u_bak = false;
u_mitt = false;
u_symp = true;
u_ode = false;


h_order1 = 0.001;
h_order2 = sqrt(h_order1);
t_start = 0;
t_slut = 100;
a = 0.2;


% Euler-fram

if u_fram == true
    t_0 = t_start;
    t_1 = t_slut;
    % Tidssteg
    h = h_order1;
    % Antal steg
    N = round((t_1 - t_0) / h);
    
    % Framåt Euler (explicit Euler)
    q1 = zeros(1, N+1);
    q2 = zeros(1, N+1);
    p1 = zeros(1, N+1);
    p2 = zeros(1, N+1);
    
    % Initialvärden
    q1(1) = 1 - a;
    q2(1) = 0;
    p1(1) = 0;
    p2(1) = sqrt((1 + a) / (1 - a));
    
    for i = 1:N
        % Uppdatera q med värdet från föregående steg
        q1(i+1) = q1(i) + h * p1(i); 
        q2(i+1) = q2(i) + h * p2(i); 
        % Beräkna avståndet r
        r = sqrt(q1(i)^2 + q2(i)^2);
        % Uppdatera p med nuvarande q
        p1(i+1) = p1(i) - h * q1(i) / (r^3); 
        p2(i+1) = p2(i) - h * q2(i) / (r^3); 
    end
    
    % Plotta framåt Euler (röd)
    figure;
    plot(q1, q2, 'r');
    hold on;
end




% Euler-bak

if u_bak == true
    t_0 = t_start;
    t_1 = t_slut;
    % Tidssteg
    h = h_order1;
    % Antal steg
    N = round((t_1 - t_0) / h);


    % Återställ
    q1 = zeros(1, N);
    q2 = zeros(1, N);
    p1 = zeros(1, N);
    p2 = zeros(1, N);
    
    q1(1) = 1 - a;
    q2(1) = 0;
    p1(1) = 0;
    p2(1) = sqrt((1 + a) / (1 - a));
    
       
    for n = 1:N
        
        q1_ny = q1(n) + h*p1(n);
        q2_ny = q2(n) + h*p2(n);
        p1_ny = p1(n);
        p2_ny = p2(n);
        
    
        for iter = 1:1000
            r_ny = sqrt(q1_ny^2 + q2_ny^2);
            p1_ny = p1(n) - h * q1_ny / r_ny^(3);
            p2_ny = p2(n) - h * q2_ny / r_ny^(3);
    
            q1_ny = q1(n) + h*p1_ny;
            q2_ny = q2(n) + h*p2_ny;  
    
        end    
        
        q1(n+1) = q1_ny;
        q2(n+1) = q2_ny;
        p1(n+1) = p1_ny;
        p2(n+1) = p2_ny;
    
    end
    
    % Plotta bakåt Euler (blå)
    plot(q1, q2, 'b');
end




% ODE-45
if u_ode == true
    a = 0.5;
    % Initialvillkor
    q1_0 = 1 - a;
    q2_0 = 0;
    p1_0 = 0;
    p2_0 = sqrt((1 + a) / (1 - a));
    
    % Tillståndsvektor
    y0 = [q1_0; q2_0; p1_0; p2_0];
    
    % Definiera ODE-systemet för Keplers problem
    ode_syst = @(t, y) [ y(3); 
                         y(4); 
                        -y(1)/( (y(1)^2 + y(2)^2)^(3/2) ); 
                        -y(2)/( (y(1)^2 + y(2)^2)^(3/2) )];
    
    % Tidsintervall
    tspan = [0, 100];
    
    % Alternativ för ode45
    options = odeset('RelTol', 1e-3, 'AbsTol', 1e-4);
    
    % Lös systemet med ode45
    [t, y] = ode45(ode_syst, tspan, y0, options);
    
    plot(y(:,1), y(:,2), 'g');
end



% Implicita mittpunktsmetod
if u_mitt == true
    t_0 = t_start;
    t_1 = t_slut;
    % Tidssteg
    h = h_order2;
    % Antal steg
    N = round((t_1 - t_0) / h);    
    
    q1 = zeros(1, N);
    q2 = zeros(1, N);
    p1 = zeros(1, N);
    p2 = zeros(1, N);
    
    q1(1) = 1 - a;
    q2(1) = 0;
    p1(1) = 0;
    p2(1) = sqrt((1 + a) / (1 - a));
    
    for n = 1:N
        
        q1_ny = q1(n) + h*p1(n);
        q2_ny = q2(n) + h*p2(n);
        p1_ny = p1(n);
        p2_ny = p2(n);
        
    
        for iter = 1:1000
            r_ny = sqrt((q1_ny + q1(n))^2 + (q2_ny+q2(n))^2)/2;
            p1_ny = p1(n) - h * (q1_ny + q1(n))/ (2* r_ny^(3));
            p2_ny = p2(n) - h * (q2_ny + q2(n))/ (2* r_ny^(3));
    
            q1_ny = q1(n) + h*(p1_ny + p1(n))/2;
            q2_ny = q2(n) + h*(p2_ny + p2(n))/2;  
        end    
        
        q1(n+1) = q1_ny;
        q2(n+1) = q2_ny;
        p1(n+1) = p1_ny;
        p2(n+1) = p2_ny;
    
    end

    % Plotta mittpunkts-metod (svart)
    plot(q1, q2, 'k');
end






if u_symp == true
    t_0 = t_start;
    t_1 = t_slut;
    % Tidssteg
    h = h_order2;
    % Antal steg

    N = round((t_1 - t_0) / h);    
    
    q1 = zeros(1, N);
    q2 = zeros(1, N);
    p1 = zeros(1, N);
    p2 = zeros(1, N);
    
    q1(1) = 1 - a;
    q2(1) = 0;
    p1(1) = 0;
    p2(1) = sqrt((1 + a) / (1 - a));
    
    for n = 1:N
        
        q1_ny = q1(n) + h*p1(n);
        q2_ny = q2(n) + h*p2(n);
        p1_ny = p1(n);
        p2_ny = p2(n);
        
    
        for iter = 1:1000
            q1_ny = q1(n) + h*(p1_ny);
            q2_ny = q2(n) + h*(p2_ny);  
            
            r_ny = sqrt((q1(n))^2 + (q2(n))^2);
            p1_ny = p1(n) - h * (q1(n))/ (r_ny^(3/2));
            p2_ny = p2(n) - h * (q2(n))/ (r_ny^(3/2));
    
        end    
        
        q1(n+1) = q1_ny;
        q2(n+1) = q2_ny;
        p1(n+1) = p1_ny;
        p2(n+1) = p2_ny;
    end


    % Plotta symplektisk (lila)
    plot(q1, q2, 'm');


end

xlabel('q_1');
ylabel('q_2');
title('q1 vs q2');
grid on;
axis equal;

str = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s', ...
    'Graferna visar (när man ökar tiden exempelvis från t = 50 ---> t = 100 att', ...
    '', ...
    'Frammåt Euler: Energi ökar med tid (spiralerar ut).', ...
    'Bakåt Euler: Energin minskar med tid (spiralerar inåt), innan q1 och q2 hamnar i samma position och får en orimlig kontinuitet,', ...
    'Mittpunkts-metod: Energin är konstant (minska värdet av a för att göra tydligare).', ...
    'Symplektisk Euler: Energin är konstant över tid.');

disp(str);