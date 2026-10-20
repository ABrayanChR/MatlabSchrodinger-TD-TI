clc, clear, close all

% anim = 1: se realiza la animación. anim = 0: no se realiza animación.
anim = 1;
simn = 1;
% simn = 1: se simula una partícula con una barrera de potencial.
% simn = 2: se simula una partícula con dos barreras de potencial.

% Propóngase una discretización según la siguiente malla
a = 0;
b = 500;
ti = 0;
tf = 2000;
dx = 1;
dt = 0.1;
N = ((b-a)/dx)+1; % Puntos del espacio
Nbar = ((tf-ti)/dt)+1; % Puntos del tiempo
x = a:dx:b;
t = ti:dt:tf;


% Se imponen las condiciones iniciales:

m = 1;
hbar = 1;
xo = 200;
q = 2; % momento lineal
sigma = 20; % varianza de la distribución espacial de la partícula
Vo = 0.7; % potencial de la barrera de potencial

% a y b corresponden a c y d de las transparencias
% mientras que c y d corresponden a e y f de las
% transparencias. Estas son las posiciones de las barreras de potencial.


c = 250;
d = 260;

% Definimos el potencial:
V = zeros(size(x));
V(x>=c & x<=d) = Vo;
if simn == 2
    e = 300;
    f = 310;
    V(x>=e & x<=f) = Vo;
end

% Definimos la condición inicial para el problema.

onda = zeros(N, Nbar);
onda(:, 1) = exp(1i * q .* x) .* exp(-((x - xo).^2) / (2 * sigma^2));


% Resolvemos la onda con la solución seleccionada

% Generamos la matriz A:
A = zeros(N, N);
for j = 1:N
    for k = 1:N
        if j == k
            A(j, k) = 2 * ((1i * 2 * m * dx^2) / (dt * hbar) - 1 - (m * dx^2) / (hbar^2) * V(j));
        elseif abs(j - k) == 1
            A(j, k) = 1;
        else
            A(j, k) = 0;
        end
    end
end

onda = solucion2(onda, V, N, Nbar, hbar, m, dx, dt, A);


if simn == 1

    figure;
    set(gca,'fontsize',8)
    ax1 = subplot(2,2,1);
    hold on
    plot(x, abs(onda(:,2)).^2)
    plot(x, V, 'r--')
    legend('|\psi(x,t)|^2', 'V(x)')
    xlabel('Posición');
    ylabel('Densidad de probabilidad |\psi(x,t)|^2');
    xlim([0, b]);
    ylim([0, 0.05]);
    grid on;
    hold off

    ax2 = subplot(2,2,2);
    set(gca,'fontsize',8)
    hold on
    plot(x, abs(onda(:,500)).^2)
    plot(x, V, 'r--')
    legend('|\psi(x,t)|^2', 'V(x)')
    xlabel('Posición');
    ylabel('Densidad de probabilidad |\psi(x,t)|^2');
    xlim([0, b]);
    ylim([0, 0.05]);
    grid on;
    hold off

    ax3 = subplot(2,2,3);
    set(gca,'fontsize',8)
    hold on
    plot(x, abs(onda(:,1000)).^2)
    plot(x, V, 'r--')
    legend('|\psi(x,t)|^2', 'V(x)')
    xlabel('Posición');
    ylabel('Densidad de probabilidad |\psi(x,t)|^2');
    xlim([0, b]);
    ylim([0, 0.05]);
    grid on;
    hold off

    ax4 = subplot(2,2,4);
    set(gca,'fontsize',8)
    hold on
    plot(x, abs(onda(:,1500)).^2)
    plot(x, V, 'r--')
    legend('|\psi(x,t)|^2', 'V(x)')
    xlabel('Posición');
    ylabel('Densidad de probabilidad |\psi(x,t)|^2');
    xlim([0, b]);
    ylim([0, 0.05]);
    grid on;
    hold off

end

if simn == 2

    figure;
    ax1 = subplot(3,2,1);
    set(gca,'fontsize',5)
    hold on
    plot(x, abs(onda(:,2)).^2)
    plot(x, V, 'r--')
    legend('|\psi(x,t)|^2', 'V(x)')
    xlabel('Posición');
    ylabel('Densidad de probabilidad |\psi(x,t)|^2');
    xlim([0, b]);
    ylim([0, 0.05]);
    grid on;
    hold off

    ax2 = subplot(3,2,2);
    set(gca,'fontsize',5)
    hold on
    plot(x, abs(onda(:,500)).^2)
    plot(x, V, 'r--')
    legend('|\psi(x,t)|^2', 'V(x)')
    xlabel('Posición');
    ylabel('Densidad de probabilidad |\psi(x,t)|^2');
    xlim([0, b]);
    ylim([0, 0.05]);
    grid on;
    hold off

    ax3 = subplot(3,2,3);
    set(gca,'fontsize',5)
    hold on
    plot(x, abs(onda(:,1000)).^2)
    plot(x, V, 'r--')
    legend('|\psi(x,t)|^2', 'V(x)')
    xlabel('Posición');
    ylabel('Densidad de probabilidad |\psi(x,t)|^2');
    xlim([0, b]);
    ylim([0, 0.05]);
    grid on;
    hold off

    ax4 = subplot(3,2,4);
    set(gca,'fontsize',5)
    hold on
    plot(x, abs(onda(:,1500)).^2)
    plot(x, V, 'r--')
    legend('|\psi(x,t)|^2', 'V(x)')
    xlabel('Posición');
    ylabel('Densidad de probabilidad |\psi(x,t)|^2');
    xlim([0, b]);
    ylim([0, 0.05]);
    grid on;
    hold off

    ax4 = subplot(3,2,5);
    set(gca,'fontsize',5)
    hold on
    plot(x, abs(onda(:,2000)).^2)
    plot(x, V, 'r--')
    legend('|\psi(x,t)|^2', 'V(x)')
    xlabel('Posición');
    ylabel('Densidad de probabilidad |\psi(x,t)|^2');
    xlim([0, b]);
    ylim([0, 0.05]);
    grid on;
    hold off

    ax5 = subplot(3,2,6);
    set(gca,'fontsize',5)
    hold on
    plot(x, abs(onda(:,2500)).^2)
    plot(x, V, 'r--')
    legend('|\psi(x,t)|^2', 'V(x)')
    xlabel('Posición');
    ylabel('Densidad de probabilidad |\psi(x,t)|^2');
    xlim([0, b]);
    ylim([0, 0.05]);
    grid on;
    hold off

end


if anim == 1
    figure;

    moduloonda = conj(onda) .* onda;
    for i = 1:Nbar
        plot(x,moduloonda(:,i))
        pause(0.01)
    end
end


% Para calcular la probabilidad de encontrar la partícula en un intervalo
% del espacio, tendremos que recorrer todos los puntos del espacio
% discretizado dentro de los pozos de potencial y sumar las probabilidades
% de encontrar la partícula en esos puntos ( cuya probabilidad es el
% cuadrado del valor absoluto de la función de onda ).


prob = probabilidad(x,t,onda,c,d,dx);






% Se proporcionan las dos funciones que resuelven el problema mediante la
% iteración y las funciones que determinan el potencial y la probabilidad
% de encontrar a la partícula en cierto intervalo



function [V] = Vcond(x,a,b,value)
V = 0;
if x <= b && x >= a
    V = value;
end
end


function [prob] = probabilidad(x,t,onda,a,b,dx)
    prob = zeros(length(t),1);
    for n = 1:length(t)
        for k = 1:length(x)
            if x(k) >= a && x(k) <= b
                prob(n) = prob(n) + abs(onda(k,n))^2*dx;
            end
        end
    end
    disp(prob')
end


function [onda] = solucion1(onda, V, N, Nbar, hbar, m, dx, dt)

for n = 1:Nbar-1
    for k = 2:N-1
        onda(n+1,k) = onda(n,k) - (1i/hbar)*dt * ( (-hbar^2/(2*m))*(onda(n,k-1)-2*onda(n,k)+onda(n,k+1))/dx^2 + V(k)*onda(n,k) );
    end
end

end



function [ondafinal] = solucion2(onda, V, N, Nbar, hbar, m, dx, dt, A)

% Resolvemos el sistema de ecuaciones.
for n = 2:Nbar
    omega = zeros(N, 1);
    omega(1) = onda(1);
    omega(N) = onda(N);

    for k = 2:N-1
        omega(k) = -onda(k-1) + 2 * ((1i * 2 * m * dx^2) / (dt * hbar) + 1 + (m * dx^2) / (hbar^2) * V(k)) * onda(k) - onda(k+1);
    end

    onda = A \ omega;

    % Normalizar onda
    normalizador = sqrt(sum(abs(onda).^2) * dx);
    onda = onda / normalizador;

    ondafinal(:, n) = onda;
end

end





