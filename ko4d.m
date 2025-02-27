N=200; % antal intervall
T=2; % sluttid
dx=1/N; % steglängd i rummet
dt=dx/2; % tidssteg, tänk på stabilitetsvillkoren
M=round(T/dt); % antal tidsteg
c = 1;
L = 1;

% allokering av minne
u=zeros(N+1,M+1); % u(n,m) lösningens värde vid tid (m-1)*dt i position n*dx
p=zeros(N+1,M+1); % p=u’

x = dx*(1:N+1)'; % x(n) är n*dx
E = zeros(1,M+1); % För att beräkna energin i varje tidssteg.

% Skapa matrisen A
v = ones(1,N+1)/dx^2;
vn = ones(1,N)/dx^2;
A = diag(v*-2) + diag(vn,1) + diag(vn,-1);
A(1,2) = 2/dx^2;
A(end,end-1) = 2/dx^2;

% Sätt begynnelsedata för u och p.
u(:, 1) = exp(-200*(x-0.5).^2); % B.V nr (6)
p(:, 1) = 0; % B.V nr (7)

% Räkna ut energin E vid tiden 0.
E(1) = 0.5*sum(p(:,1).^2) - 0.5*c^2*(u(:,1)'*(A*u(:,1)));

X = [0; x; L]; 

nframe = M+1; % kommando för film
mov(1:nframe) = struct('cdata',[],'colormap',[]);
figure;
% Lägg till randvärden så att längden blir N+1, vilket matchar X.
plot(X, [0; u(:,1); 0], 'b', 'Linewidth', 1); % Plot vid tiden t=0.
axis([0 1 -1 1])
set(gca, 'nextplot', 'replacechildren')
drawnow

mov(1) = getframe(gcf); % Första frame i filmen.
for m = 1:M % tidstegning med symplektisk Euler

    p(:, m+1) = p(:, m) + c^2 * dt * A * u(:, m);

    u(:, m+1) = u(:, m) + dt * p(:, m+1);

    % Plotta även lösningen från d’Alemberts formel
    figure(1);
    U = [0; u(:,m); 0];
    plot(X, U, 'b', 'Linewidth', 1)
    hold on;
    % Plotta även lösningen från d’Alemberts formel 
    g = @(x) exp(-200*(x-0.5).^2);
    u_dlambert = 0.5 * (g(x + m*dt) + g(x - m*dt));
    U_dlambert = [0; u_dlambert; 0];
    plot(X, U_dlambert, "r", 'LineWidth', 1);
    hold on;

    % Använd m*dt för att visa aktuell tid
    text(0.05, -0.8, sprintf('t=%.2f', m*dt))
    set(gca, 'nextplot', 'replacechildren')
    drawnow
    mov(m+1) = getframe(gcf);
    
    % Räkna ut energin av den numeriska lösningen vid detta tidssteg
    E(m+1) = 0.5*sum(p(:,m+1).^2) - 0.5*c^2*(u(:,m+1)'*(A*u(:,m+1)));
end

% Korrigerad plottning av energin: använd en tidsvektor med rätt längd.
figure(2); % Se till att vi plottar i rätt figur.
plot((0:m)*dt, E(1:m+1), "r", "LineWidth", 1); % Plotta hela energihistoriken.
xlabel('Tid');
ylabel('Energi');
title('Energi över Tid');