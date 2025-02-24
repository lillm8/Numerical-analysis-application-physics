compute_energy = @(q1, q2, p1, p2) 0.5*(p1.^2 + p2.^2) - 1./sqrt(q1.^2 + q2.^2);

% Settings and booleans
% Ordning 1:
u_fram = false;   % Forward Euler
u_bak  = false;   % Backward Euler
u_symp = true;    % Symplectic Euler

% Ordning 2:
u_mitt = false;   % Implicit Midpoint

% Övrigt:
u_ode  = true;    % ODE45

h_order1 = 0.005;
h_order2 = 0.05;
t_start = 0;
t_slut  = 100;
a = 0.5;

% Plot the trajectories (q1 vs q2) as in the original code
figure;
hold on;

% Forward Euler (Explicit Euler)
if u_fram
    t0 = t_start;
    t1 = t_slut;
    h = h_order1;
    N = round((t1 - t0) / h);
    
    q1 = zeros(1, N+1); q2 = zeros(1, N+1);
    p1 = zeros(1, N+1); p2 = zeros(1, N+1);
    
    % Initial values
    q1(1) = 1 - a;
    q2(1) = 0;
    p1(1) = 0;
    p2(1) = sqrt((1 + a) / (1 - a));
    
    for i = 1:N
        q1(i+1) = q1(i) + h * p1(i);
        q2(i+1) = q2(i) + h * p2(i);
        r = sqrt(q1(i)^2 + q2(i)^2);
        p1(i+1) = p1(i) - h * q1(i) / r^3;
        p2(i+1) = p2(i) - h * q2(i) / r^3;
    end
    
    plot(q1, q2, 'r', 'DisplayName', 'Forward Euler');
    % Compute energy and time vector for energy plot:
    t_fram = linspace(t_start, t_slut, length(q1));
    energy_fram = compute_energy(q1, q2, p1, p2);
end

% Backward Euler (Implicit Euler)
if u_bak
    t0 = t_start;
    t1 = t_slut;
    h = h_order1;
    N = round((t1 - t0) / h);
    
    q1 = zeros(1, N+1); q2 = zeros(1, N+1);
    p1 = zeros(1, N+1); p2 = zeros(1, N+1);
    
    q1(1) = 1 - a;
    q2(1) = 0;
    p1(1) = 0;
    p2(1) = sqrt((1 + a) / (1 - a));
    
    for n = 1:N
        q1_ny = q1(n) + h * p1(n);
        q2_ny = q2(n) + h * p2(n);
        p1_ny = p1(n);
        p2_ny = p2(n);
        
        for iter = 1:100
            r_ny = sqrt(q1_ny^2 + q2_ny^2);
            p1_ny = p1(n) - h * q1_ny / r_ny^3;
            p2_ny = p2(n) - h * q2_ny / r_ny^3;
            q1_ny = q1(n) + h * p1_ny;
            q2_ny = q2(n) + h * p2_ny;
        end
        
        q1(n+1) = q1_ny;
        q2(n+1) = q2_ny;
        p1(n+1) = p1_ny;
        p2(n+1) = p2_ny;
    end
    
    plot(q1, q2, 'b', 'DisplayName', 'Backward Euler');
    t_bak = linspace(t_start, t_slut, length(q1));
    energy_bak = compute_energy(q1, q2, p1, p2);
end

% ODE45 (for comparison)
if u_ode
    % Initial conditions
    q1_0 = 1 - a;
    q2_0 = 0;
    p1_0 = 0;
    p2_0 = sqrt((1 + a) / (1 - a));
    y0 = [q1_0; q2_0; p1_0; p2_0];
    
    % Define the ODE system (Kepler's problem)
    ode_syst = @(t, y) [ y(3);
                         y(4);
                        -y(1)/( (y(1)^2 + y(2)^2)^(3/2) );
                        -y(2)/( (y(1)^2 + y(2)^2)^(3/2) )];
    
    % Time span and options
    tspan = [0, 100];
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    
    [t, y] = ode45(ode_syst, tspan, y0, options);
    plot(y(:,1), y(:,2), 'g', 'DisplayName', 'ODE45');
    
    t_ode = t;
    energy_ode = compute_energy(y(:,1), y(:,2), y(:,3), y(:,4));
    
    % Save final position (if needed in the output string)
    losningsvektor_ode45 = [y(end,1), y(end,2)];
end

% Implicit Midpoint Method
if u_mitt
    t0 = t_start;
    t1 = t_slut;
    h = h_order2;
    N = round((t1 - t0) / h);
    
    q1 = zeros(1, N+1); q2 = zeros(1, N+1);
    p1 = zeros(1, N+1); p2 = zeros(1, N+1);
    
    q1(1) = 1 - a;
    q2(1) = 0;
    p1(1) = 0;
    p2(1) = sqrt((1 + a) / (1 - a));
    
    for n = 1:N
        q1_ny = q1(n) + h * p1(n);
        q2_ny = q2(n) + h * p2(n);
        p1_ny = p1(n);
        p2_ny = p2(n);
        
        for iter = 1:1000
            r_ny = sqrt((q1_ny + q1(n))^2 + (q2_ny + q2(n))^2) / 2;
            p1_ny = p1(n) - h * (q1_ny + q1(n)) / (2 * r_ny^3);
            p2_ny = p2(n) - h * (q2_ny + q2(n)) / (2 * r_ny^3);
            q1_ny = q1(n) + h * (p1_ny + p1(n)) / 2;
            q2_ny = q2(n) + h * (p2_ny + p2(n)) / 2;
        end
        
        q1(n+1) = q1_ny;
        q2(n+1) = q2_ny;
        p1(n+1) = p1_ny;
        p2(n+1) = p2_ny;
    end
    
    plot(q1, q2, 'k', 'DisplayName', 'Implicit Midpoint');
    t_mitt = linspace(t_start, t_slut, length(q1));
    energy_mitt = compute_energy(q1, q2, p1, p2);
end

% Symplectic Euler Method
if u_symp
    t0 = t_start;
    t1 = t_slut;
    h = h_order1;
    N = round((t1 - t0) / h);
    
    q1 = zeros(1, N+1); q2 = zeros(1, N+1);
    p1 = zeros(1, N+1); p2 = zeros(1, N+1);
    
    q1(1) = 1 - a;
    q2(1) = 0;
    p1(1) = 0;
    p2(1) = sqrt((1 + a) / (1 - a));
    
    for n = 1:N
        q1_ny = q1(n) + h * p1(n);
        q2_ny = q2(n) + h * p2(n);
        p1_ny = p1(n);
        p2_ny = p2(n);
        
        for iter = 1:1000
            % In symplectic Euler, update q first with the new p
            q1_ny = q1(n) + h * p1_ny;
            q2_ny = q2(n) + h * p2_ny;
            r_ny = sqrt(q1(n)^2 + q2(n)^2);
            p1_ny = p1(n) - h * q1(n) / r_ny^3;
            p2_ny = p2(n) - h * q2(n) / r_ny^3;
        end
        
        q1(n+1) = q1_ny;
        q2(n+1) = q2_ny;
        p1(n+1) = p1_ny;
        p2(n+1) = p2_ny;
    end
    
    plot(q1, q2, 'm--', 'DisplayName', 'Symplectic Euler');
    t_symp = linspace(t_start, t_slut, length(q1));
    energy_symp = compute_energy(q1, q2, p1, p2);
    
    % Save final position (if needed in the output string)
    losningsvektor_symp = [q1(end), q2(end)];
end

xlabel('q_1');
ylabel('q_2');
title('Trajectory: q_1 vs q_2');
grid on;
axis equal;
hold off;

% Now, plot energy as a function of time for the selected methods

figure;

if u_fram
    plot(t_fram, energy_fram, 'r', 'DisplayName', 'Forward Euler');
end
if u_bak
    plot(t_bak, energy_bak, 'b', 'DisplayName', 'Backward Euler');
end
if u_mitt
    plot(t_mitt, energy_mitt, 'k', 'DisplayName', 'Implicit Midpoint');
end
if u_symp
    plot(t_symp, energy_symp, 'm--', 'DisplayName', 'Symplectic Euler');
end
if u_ode
    plot(t_ode, energy_ode, 'g', 'DisplayName', 'ODE45');
end

losningsvektor_ode45_str = '';
losningsvektor_symp_str = '';

if u_ode
    losningsvektor_ode45_str = mat2str(losningsvektor_ode45);
end
if u_symp
    losningsvektor_symp_str = mat2str(losningsvektor_symp);
end

str = sprintf(['Graferna visar (när man ökar tiden exempelvis från t = 50 ---> t = 100 att\n' ...
    '\nFrammåt Euler: Energi ökar med tid (spiralerar ut).\n' ...
    'Bakåt Euler: Energin minskar med tid (spiralerar inåt), innan q1 och q2 hamnar i samma position och får en orimlig kontinuitet,\n' ...
    'Mittpunkts-metod: Energin är konstant (minska värdet av a för att göra tydligare).\n' ...
    'Symplektisk Euler: Energin är konstant över tid.\n' ...
    '\nVärdet av (y1,y2) är =\nODE45: %s\nSymplektisk Euler: %s'], ...
    losningsvektor_ode45_str, losningsvektor_symp_str);
disp(str);