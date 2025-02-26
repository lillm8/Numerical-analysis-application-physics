% Ordning 1:
u_fram = false;
u_bak = false;
u_symp = true;

% Ordning 2:
u_mitt = false;

% Övrigt:
u_ode = false;


h_order1 = 0.005;

h_order2 = 0.0005;
t_start = 0;
t_slut = 50;
a = 0.5;
figure;



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
    plot(q1, q2, 'r');
    hold on;
end




% Euler-bak

if u_bak == true
    t_0 = t_start;
    t_1 = t_slut;
    % Tidssteg
    h = h_order2;
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

        next = [q1(n); q2(n); p1(n); p2(n)];
    
        iter = 0;
        maxiter = 100;
        tol = 1e-8; 
    
    
        while (iter < maxiter)
            iter = iter + 1;

            r_next = next(1)^2 + next(2)^2;
              
            dF3_dq1 = h * (2 * next(1)^2 - next(2)^2)/(r_next^(5/2));
            dF3_dq2 = h * (3 * next(2) * next(1)) / (r_next^(5/2));
            dF4_dq1 = h * (3 * next(2) * next(1)) / (r_next^(5/2));
            dF4_dq2 = h * (2 * next(2)^2 - next(1)^2) / (r_next^(5/2));
            
            F = [-next(1) + q1(n) + next(3)*h; 
                -next(2) + q2(n) + next(4)*h;
                -next(3) + p1(n) - next(1)*h/(r_next^(3/2));
                -next(4) + p2(n) - next(2)*h/(r_next^(3/2))];
    
            J = [  -1,      0,      h, 0;
                    0,      -1,     0, h;
                 dF3_dq1, dF3_dq2, -1, 0;
                 dF4_dq1, dF4_dq2, 0, -1];
    
            % Newton step
            s = -J \ F;
            next = next + s;
    
            if norm(s, inf) < tol
                break;
            end   
        end
    
        % Om Newtons metod inte konvergerar, ge en varning
        if iter == maxiter
            warning('Newton-metoden nådde max antal iterationer vid steg %d', n);
        end
    
        q1(n+1) = next(1);
        q2(n+1) = next(2);
        p1(n+1) = next(3);
        p2(n+1) = next(4); 
    end

    
    % Plotta bakåt Euler (blå)
    plot(q1, q2, 'b');
    hold on;

end




% ODE-45
if u_ode == true
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
    options = odeset('RelTol', 1e-6,'AbsTol', 1e-6);
    
    % Lös systemet med ode45
    [t, y] = ode45(ode_syst, tspan, y0, options);
    
    plot(y(:,1), y(:,2), 'g');
    hold on;
    
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
        
        iter = 0;
        maxiter = 100;
        tol = 1e-8;
        
        % Dessa ska vara 4x1 vektorer på samma sätt som q1 är Nx1 vektorer.
        next = [q1(n); q2(n); p1(n); p2(n)];

        while (iter < maxiter)
            iter = iter + 1;
            
            next_halv = (next + [q1(n); q2(n); p1(n); p2(n)])/2;
            r_next_halv = next_halv(1)^2 + next_halv(2)^2;

            dF3_dq1 = (h/2) * (2 * next_halv(1)^2 - next_halv(2)^2)/(r_next_halv^(5/2));
            dF3_dq2 = (h/2) * (3 * next_halv(2) * next_halv(1))/(r_next_halv^(5/2));
            dF4_dq1 = (h/2) * (3 * next_halv(2) * next_halv(1))/(r_next_halv^(5/2));
            dF4_dq2 = (h/2) * (2 * next_halv(2)^2 - next_halv(1)^2)/(r_next_halv^(5/2));

            F =[-next(1) + q1(n) + h*next_halv(3); 
                -next(2) + q2(n) + h*next_halv(4);
                -next(3) + p1(n) - h * next_halv(1)/r_next_halv^(3/2);
                -next(4) + p2(n) - h * next_halv(2)/r_next_halv^(3/2)];
    
            J = [  -1,      0,      h/2, 0;
                    0,      -1,     0,  h/2;
                 dF3_dq1, dF3_dq2, -1,  0;
                 dF4_dq1, dF4_dq2,  0,  -1];

            s = -J \ F;
            next = next + s;
            
            if norm(s, inf) < tol     
                break;
            end

        end

        if iter == maxiter
            warning('Newton-metoden nådde max antal iterationer vid steg - Mittpunktsiteration. %d', n);
        end

        q1(n+1) = next(1);
        q2(n+1) = next(2);
        p1(n+1) = next(3);
        p2(n+1) = next(4);
    end

    % Plotta mittpunkts-metod (svart)
    plot(q1, q2, 'k');
    hold on;

end





% Symplektisk euler-metod

if u_symp == true
    t_0 = t_start;
    t_1 = t_slut;
    % Tidssteg
    h = h_order1;
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
    
    for n = 1:N-1
        % Uppdatera impulsen först (p)
        r = (q1(n)^2 + q2(n)^2)^(3/2);
        p1(n+1) = p1(n) - h * (q1(n) / r);
        p2(n+1) = p2(n) - h * (q2(n) / r);
    
        % Uppdatera positionen med nya impulsen
        q1(n+1) = q1(n) + h * p1(n+1);
        q2(n+1) = q2(n) + h * p2(n+1);
    end


    % Plotta symplektisk (lila)
    plot(q1, q2, 'm--');
    hold on;    
end



xlabel('q_1');
ylabel('q_2');
title('q1 vs q2');
grid on;
axis equal;


str = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s', ...
    'Graferna visar (när man ökar tiden exempelvis från t = 50 ---> t = 100 att', ...
    '', ...
    'Frammåt Euler: Energi ökar med tid (spiralerar ut).', ...
    'Bakåt Euler: Energin minskar med tid (spiralerar inåt), innan q1 och q2 hamnar i samma position och får en orimlig kontinuitet,', ...
    'Mittpunkts-metod: Energin är konstant (minska värdet av a för att göra tydligare).', ...
    'Symplektisk Euler: Energin är konstant över tid.');


disp(str);