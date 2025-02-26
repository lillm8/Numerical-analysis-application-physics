t_0 = 0;
t_1 = 500;

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
tspan = [t_0, t_1];

% Alternativ för ode45
options = odeset('RelTol', 1e-15, 'AbsTol', 1e-18);

% Lös systemet med ode45
[t, y] = ode45(ode_syst, tspan, y0, options);


% Antal tidssteg
nu = 100;

tidstegen = linspace(0.005, 0.0005, nu);

% q1 = y(:,1), q2 = y(:,2)


errors_mit = zeros(1,nu);
errors_symp = zeros(1,nu);



for i = 1:length(tidstegen)

    % loop för mittpunkts_____________________________________

    t_0 = tspan(1);
    t_1 = tspan(2);
    % Tidssteg
    h = tidstegen(i);
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
    
        for iter = 1:20
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

    y_last_ode = [y(end,1), y(end,2)];
    y_last_mittpunkts = [q1(end), q2(end)];
    
    error = norm(y_last_mittpunkts - y_last_ode);
    errors_mit(i) = error;



    % Loop för symplektisk_________________________________

    
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
        
    
        for iter = 1:20
            q1_ny = q1(n) + h*(p1_ny);
            q2_ny = q2(n) + h*(p2_ny);  
            
            r_ny = sqrt((q1(n))^2 + (q2(n))^2);
            p1_ny = p1(n) - h * (q1(n))/ (r_ny^(3));
            p2_ny = p2(n) - h * (q2(n))/ (r_ny^(3));
        end    
        
        q1(n+1) = q1_ny;
        q2(n+1) = q2_ny;
        p1(n+1) = p1_ny;
        p2(n+1) = p2_ny;
    end
    
    y_last_symp = [q1(end), q2(end)];
    error = norm(y_last_symp - y_last_ode);
    errors_symp(i) = error;

end


mit_p = zeros(nu,1);
symp_p = zeros(nu,1);

for i = 1:(length(errors_mit)-1)
        mit_p(i) = log(errors_mit(i)/errors_mit(i+1))/log(tidstegen(i)/tidstegen(i+1));
        symp_p(i) = log(errors_symp(i)/errors_symp(i+1))/log(tidstegen(i)/tidstegen(i+1));
end

mit_p
symp_p