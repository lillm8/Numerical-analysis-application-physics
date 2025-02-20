a = 0.5;
% Starttid
t_0 = 0;
% Sluttid
t_1 = 100;
% Tidssteg
h = 0.05;
% Antal steg
N = (t_1 - t_0) / h;

% Tidsvektor
t = linspace(t_0, t_1, N+1);

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

% Återställ
q1 = zeros(1, N);
q2 = zeros(1, N);
p1 = zeros(1, N);
p2 = zeros(1, N);

q1(1) = 1 - a;
q2(1) = 0;
p1(1) = 0;
p2(1) = sqrt((1 + a) / (1 - a));


% Bakåt euler

for n = 1:N
    
    q1_ny = q1(n) + h*p1(n);
    q2_ny = q2(n) + h*p2(n);
    p1_ny = p1(n);
    p2_ny = p2(n);
    

    for iter = 1:1000
        r_ny = sqrt(q1_ny^2 + q2_ny^2);
        p1_ny = p1(n) - h * q1_ny / r_ny^(3/2);
        p2_ny = p2(n) - h * q2_ny / r_ny^(3/2);

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
tspan = [0, 500];

% Alternativ för ode45
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

% Lös systemet med ode45
[t, y] = ode45(ode_syst, tspan, y0, options);

plot(y(:,1), y(:,2), 'g');

legend('Framåt Euler', 'Bakåt Euler', 'Ode45 (exakt)');
xlabel('q_1');
ylabel('q_2');
title('Jämförelse mellan Framåt och Bakåt Euler');
grid on;
axis equal;
