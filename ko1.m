a = 0.7;
% Starttid
t_0 = 0;
% Sluttid
t_1 = 1000;
% Tidssteg
h = 0.5;
% Antal steg
N = (t_1 - t_0)/h;
% Frammåt euler
t = linspace(t_0,t_1,N+1);
q1 = zeros(1, N);
q2 = zeros(1, N);
p1 = zeros(1, N);
p2 = zeros(1, N);

q1(1) = 1 - a;
q2(1) = 0;
p1(1) = 0;
p2(1) = sqrt((1 + a) / (1 - a));

for i = 1:N
    q1(i+1) = q1(i) + h*p1(i); 
    q2(i+1) = q2(i) + h*p2(i); 
    p1(i+1) = p1(i) - h*q1(i)*(q1(i)^2 + q2(i)^2)^(-3/2); 
    p2(i+1) = p2(i) - h*q2(i)*(q1(i)^2 + q2(i)^2)^(-3/2); 
end

figure;
plot(q1,q2, 'r-'); 
title('Euler');
xlabel('Tid');
ylabel('Position');
legend('q1', 'q2');



a1 = q1(2);
a2 = q2(2);

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
    

    for iter = 1:100     
        r_ny = sqrt(q1_ny^2 + q2_ny^2);
        p1_ny = p1(n) - h * q1_ny / r_ny^(3/2);
        p2_ny = p1(n) - h * q2_ny / r_ny^(3/2);

        q1_ny = q1(n) + h*p1_ny;
        q2_ny = q2(n) + h*p2_ny;  

    end    
    
    q1(n+1) = q1_ny;
    q2(n+1) = q2_ny;
    p1(n+1) = p1_ny;
    p2(n+1) = p2_ny;

end

plot(q1, q2, 'r-'); 
title('Euler fram [ÖVER TID]');
xlabel('q1');
ylabel('q2');
grid on;
axis equal;

% Implicita mittpunktsmetoden




% Symplektisk euler
