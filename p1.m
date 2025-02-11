A = [0, 0, 1, 0;
    0, 0, 0, 1;
    -(q1^2 + q2^2)^(-3/2), 0, 0, 0;
    0, -(q1^2 + q2^2)^(-3/2), 0, 0;
];


a = 0.5;
% Starttid
t_0 = 0;
% Sluttid
t_1 = 100;
% Tidssteg
h = 0.05;
% Antal steg
N = (t_1 - t_0)/h;
% Frammåt euler
q1 = linspace(t_0,t_1,N);
q2 = linspace(t_0,t_1,N);
p1 = linspace(t_0,t_1,N);
p2 = linspace(t_0,t_1,N);
q1(0) = 1 - a;
q2(0) = 0;
p1(0) = 0;
p2(0) = sqrt((1 + a) / (1 - a));
y_{n+1} = y_n + h*f(y_n);


% Bakåt euler


% Implicita mittpunktsmetoden



% Symplektisk euler

