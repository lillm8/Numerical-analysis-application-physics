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

figure;
plot(y(:,1), y(:,2), 'b');
xlabel('q_1');
ylabel('q_2');
grid on;
