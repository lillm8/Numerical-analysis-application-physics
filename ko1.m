u_fram = false;
u_bak = false;
u_ode = false;
u_mitt = true;
u_symp = false;

h_order1 = 0.005;
h_order2 = 0.0005;

h_f = h_order2;
h_b = h_order2;
h_o = h_order2;
h_m = h_order2;
h_s = h_order2;

t_0 = 0;
t_1 = 100;

N1_deg = round((t_1 - t_0)/h_order1);
N2_deg = round((t_1 - t_0)/h_order2);

N_f = N2_deg;
N_b = N2_deg;
N_o = N2_deg;
N_m = N2_deg;
N_s = N2_deg;

a = 0.5;

legend_lable = {};





figure;
% Euler-fram

if u_fram == true
    % Tidssteg
    h = h_f;
    % Antal steg
    N = N_f;

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
        r = sqrt(q1(i)^2 + q2(i)^2);
        % Uppdatera q med värdet från föregående steg
        q1(i+1) = q1(i) + h * p1(i); 
        q2(i+1) = q2(i) + h * p2(i); 
        % Beräkna avståndet r
        % Uppdatera p med nuvarande q
        p1(i+1) = p1(i) - h * q1(i) / (r^3); 
        p2(i+1) = p2(i) - h * q2(i) / (r^3);      
    end
    
    q1_f = q1;
    q2_f = q2;
    p1_f = p1;
    p2_f = p2;

    % Plotta framåt Euler (röd)
    plot(q1, q2, 'b');
    legend_lable{end + 1} = "Forward";
    hold on;
end




% Euler-bak

if u_bak == true
    % Tidssteg
    h = h_b;
    % Antal steg
    N = N_b;

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
        
        % Iteration och tolerans
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

    q1_b = q1;
    q2_b = q2;
    p1_b = p1;
    p2_b = p2;
    
    % Plotta bakåt Euler (blå)
    plot(q1, q2, 'r');
    legend_lable{end + 1} = "Backward";
    hold on;
end




% ODE-45
if u_ode == true
    % Antal steg
    N = N_o;

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
    legend_lable{end + 1} = "ODE45";
    hold on;
end



% Implicita mittpunktsmetod
if u_mitt == true
    % Tidssteg
    h = h_m;
    % Antal steg
    N = N_m;
    
    % Energy middle-point euler
    energy_mid = zeros(1,N);
    energy_mid_t = linspace(t_0,t_1, N);

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
 
    q1_m = q1;
    q2_m = q2;
    p1_m = p1;
    p2_m = p2;

    % Plotta mittpunkts-metod (svart)
    plot(q1, q2, 'k');
    legend_lable{end + 1} = "Middle-point";
    hold on;

end




% Symplektisk euler-metod

if u_symp == true
    % Tidssteg
    h = h_s;
    % Antal steg
    N = N_s;

    % Energi symp euler
    energy_symp = zeros(1,N);
    energy_symp_t = linspace(t_0,t_1, N);
    

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
   
    q1_s = q1;
    q2_s = q2;
    p1_s = p1;
    p2_s = p2;


    % Plotta symplektisk (lila)
    plot(q1, q2, 'm-');
    legend_lable{end + 1} = "Symplectic";
    hold on;    
end

xlabel('q_1');
ylabel('q_2');
title('q1 vs q2 for different methods');
legend(legend_lable)
grid on;
axis equal;



%____________________PLOTTA ENERGIERNA_____________________
figure;

% Energi frammåt euler
if u_fram == true
    energy_fram = zeros(1,N_f);
    energy_fram_t = linspace(t_0,t_1, N_f);

    for n = 1:N_f
        energy_fram(n) = (1/2)*(p1_f(n)^2 + p2_f(n)^2) - 1/sqrt(q1_f(n)^2 + q2_f(n)^2);
    end

    plot(energy_fram_t, energy_fram, "b");
    hold on;
end


% Energi bakåt euler
if u_bak == true
    energy_back = zeros(1,N_b);
    energy_back_t = linspace(t_0,t_1, N_b);

    for n = 1:N_b
        energy_back(n) = (1/2)*(p1_b(n)^2 + p2_b(n)^2) - 1/sqrt(q1_b(n)^2 + q2_b(n)^2);
    end

    plot(energy_back_t, energy_back, "r");
    hold on;
end


% Energi ODE-45
if u_ode == true    
    energy_ode = zeros(1,length(y(1)));
    energy_ode_t = t;

    for n = 1:length(t)
        energy_ode(n) = (1/2)*(y(n,3)^2 + y(n,4)^2) - 1/sqrt(y(n,1)^2 + y(n,2)^2);
    end

    plot(energy_ode_t, energy_ode, "g");
    hold on;
end

% Energi mittpunktsmetod
if u_mitt == true
    energy_mid = zeros(1,N_m);
    energy_mid_t = linspace(t_0,t_1,N_m);

    for n = 1:N_m
        energy_mid(n) = (1/2)*(p1_m(n)^2 + p2_m(n)^2) - 1/sqrt(q1_m(n)^2 + q2_m(n)^2);
    end

    plot(energy_mid_t, energy_mid, "k");
    hold on;
end

if u_symp == true
    energy_symp = zeros(1,N_s);
    energy_symp_t = linspace(t_0,t_1,N_s);

    for n = 1:N_s
        energy_symp(n) = (1/2)*(p1_s(n)^2 + p2_s(n)^2) - 1/sqrt(q1_s(n)^2 + q2_s(n)^2);
    end

    plot(energy_symp_t, energy_symp, "m-");
    hold on;
end

legend(legend_lable);
title("Energy vs time for different methods")
xlabel('Time');
ylabel('Energy');
grid on;
axis fill;

