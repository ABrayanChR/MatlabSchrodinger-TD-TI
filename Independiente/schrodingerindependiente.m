clear, clc, close all

% ECUACIÓN DE SCHRÖDINGER INDEPENDIENTE DEL TIEMPO

% Solución de la ecuación de Schrödinger independiente del tiempo
% para un pozo de potencial con diversos potenciales. Resolución del
% problema de valores propios.


% Se definen las constantes y las variables
hbar = 1.054571817646156e-34;
M = 1;    %Kg

% Condiciones de contorno (e iniciales)

K = 1;
ya = 0;
yb = 0;
dya = 1;
L = 1;

% Se define el potencial (no normalizado) y las condiciones del problema con E y V normalizados.

V = @(x) (1/2)*K*x^2;

% Donde (hbar^2/(M*L^2)) es el factor usado para la normalización (está
% multiplicando para que poder usar valores normalizados de E en la
% búsqueda de autovalores) y donde 

k = @(x, E) (2*M/hbar^2)* ( E*(hbar^2/(M*L^2)) - (hbar^2/(M*L^2))* V(x) );


% Tomando a = 0 y b = L, podremos resolver el problema en el intervalo [0,L]

a = 0;
b = L;

% Ahora iteramos con el método para obtener todos los valores de la onda en
% cada posición con un intervalo h

N = 10000; % Número de puntos
a = a/L;
b = b/L;
h = (b - a) / (N-1);
x = a:h:b;

% Elegimos dos energías

% % % % % Útil para conocer los ceros de la función de onda
% % % % E = -1:1:1000;
% % % % vlue = [];
% % % %
% % % %
% % % % for i = 1:length(E)
% % % %     y = numerov(k, E(i), N, x, h, ya, dya);
% % % %     y = y(N);
% % % %     vlue = [vlue, y];
% % % % end
% % % % plot(E,vlue,'bo')



% Se elige un valor de la energía y se procede a calcular con los valores
% de la función de onda entre a y b

E1p = 1;
E2p = E1p + 1;
nvalores = 10;
incertidumbre = 1e-6;


[E, iteraciones, incertidumbres] = biseccion(k, N, x, h, ya, dya, E1p, E2p, incertidumbre, nvalores);

% Creamos una tabla para presentar resultados de las 10 energías calculadas


%Resultado analítico de E para la partícual encerrada en una caja
Ean = ( ( (1:nvalores).^2 )*pi^2 ) / 2;

T = table( [1:nvalores]',Ean', E', iteraciones', incertidumbres', 'VariableNames', {'n','E analítica (caja V=0)', 'E numérica', 'iteraciones', 'incertidumbre'}, 'RowNames', arrayfun(@num2str, 1:nvalores, 'UniformOutput', false));

disp(T)

writetable(T,'autovalores.csv');


% Para obtener las autofunciones, aplicamos numerov para obtenerlas
% numéricamente. Así para cada autovalor E(n) obtenido.

ondas = zeros(length(E),N);
for i = 1:length(E)
    ondas(i,:) = numerov(k, E(i), N, x, h, ya, dya); 
end

% Para normalizar estas ondas, simplemente realizamos la integral sobre
% cada valor (en este caso es una suma), debiendo de dar 1. Luego,
% dividimos esa onda por su integral.


for i = 1:length(E)
    integral = h * sum(ondas(i,1:N-1).^2);
    ondas(i,:) = ondas(i,:)/sqrt(integral);
end


%%%%%%%%%%%%%% Representación de la densidad de probabilidad

figure;
title("Representación de las funciones de onda y su probabilidad")
for i = 1:length(E)
    ax = subplot(5,2,i);
    set(ax,'fontsize',9)
    hold on
    plot(x, ondas(i,:), 'b')
    plot(x, abs(ondas(i,:)).^2, 'r--')
    if i == floor(length(E)/2)
        legend('\psi(x,t)','|\psi(x,t)|^2',Location='bestoutside')
    end
    xlabel('Posición');
    tituloax = sprintf("n = %d",i);
    title(ax,tituloax);
    grid on;
    hold off
end

%%%%%%%%%%%%%% Determinación del valor medio de la posición

xmedio = zeros(1,length(E));

for i = 1:length(E)
   
    xmedio(i) = h * sum(x(1:N-1) .* abs(ondas(i,1:N-1)).^2);

end

T = table([1:length(E)]', xmedio', 'VariableNames', {'n','xmedio'}, 'RowNames', arrayfun(@num2str, 1:length(E), 'UniformOutput', false));

writetable(T,'xmedio.csv');

disp(T)

%%%%%%%%%%%%%%%%%%%%%%%%%%% Análisis del error relativo y el error relativo medio


% Propongamos 3 valores distintos de h
hvalores = [0.1, 0.001, 0.0001];

for i = 1:length(hvalores)
    % No usamos clear aquí para evitar borrar funciones o handles
    hi = hvalores(i);
    x = a:hi:b;
    Ni = length(x);
    incertidumbre = 1e-6;
    E1p = 1;
    E2p = E1p + 1;
    
    % Obtenemos los 5 primeros autovalores
    [E,~, ~] = biseccion(k, Ni, x, hi, ya, dya, E1p, E2p, incertidumbre, 5);
    
    % Se calculan las autofunciones generadas por los autovalores
    ondas_num = zeros(2, Ni);
    ondas_num(1,:) = numerov(k, E(1), Ni, x, hi, ya, dya); 
    ondas_num(2,:) = numerov(k, E(5), Ni, x, hi, ya, dya); 
    
    % Normalizamos las autofunciones
    for idx = 1:2
        integral = hi * sum(ondas_num(idx,1:Ni-1).^2);
        ondas_num(idx,:) = ondas_num(idx,:) / sqrt(integral);
    end
    
    % Calculamos la solución teórica
    ondasteoricas = zeros(2, Ni);
    ondasteoricas(1,:) = sqrt(2/L) * sin(1*pi*x/L);
    ondasteoricas(2,:) = sqrt(2/L) * sin(5*pi*x/L);
    
    % Calculamos los errores relativos
    Err = zeros(2, Ni); % Error relativo
    for idx = [1,2]
        for j = 1:Ni
            if abs(ondasteoricas(idx,j)) > 0
                Err(idx,j) = abs(ondas_num(idx,j) - ondasteoricas(idx,j)) / abs(ondasteoricas(idx,j));
            else
                Err(idx,j) = 0;
            end
        end
    end
    
    % Calcular el error relativo medio para cada autofunción
    Erm = zeros(2,1);
    for idx = [1,2]
        Erm(idx) = mean(Err(idx,:));
    end
    
    fprintf("Erm(n=1, h=%f) = %f\n", hi, Erm(1))
    fprintf("Erm(n=5, h=%f) = %f\n\n", hi, Erm(1))
    
    % Graficar el error relativo para cada función
    figure;
    subplot(2,1,1)
    plot(x, Err(1,:), 'bo');
    titulo = sprintf("n = %d, h = %f",1, hi);
    title(titulo)
    xlabel("X")
    ylabel("Error relativo")
    yscale("log")
    legend('E(n=1)',Location='best')

    subplot(2,1,2)
    plot(x, Err(2,:), 'ro');
    titulo = sprintf("n = %d, h = %f",5, hi);
    title(titulo)
    xlabel("X")
    ylabel("Error relativo")
    yscale("log")
    legend('E(n=5)',Location='best')

    titulo = sprintf("simh%f.png",hi);
    saveas(gcf,titulo)

    

end



% Método de Numerov

function [y] = numerov(k, E, N, x, h, ya, dya)

y = zeros(1, N);
y(1) = ya;
y(2) = 2*h*dya;

for i = 2:N-1
    y(i+1) = (2*y(i)*(1 - (5*h^2)*k(x(i), E)/12) - y(i-1)*(1 + (h^2)*k(x(i-1), E)/12)) / (1 + (h^2)*k(x(i+1), E)/12);
end


end



% Método de bisección


function [E, iteracionesv, incertidumbres] = biseccion(k, N, x, h, ya, dya, E1, E2, incertidumbre, nvalores)
y1 = 10000;
y2 = 0;
E = zeros(1,nvalores);
iteracionesv = E;
incertidumbres = E;
for i = 1:nvalores
    iteraciones = 1;
    while abs(y1 - y2) > incertidumbre && abs(E1 - E2) > incertidumbre && iteraciones < 10000
        %Calculamos los valores finales de la función de onda para dos energías
        y1 = numerov(k, E1, N, x, h, ya, dya);
        y2 = numerov(k, E2, N, x, h, ya, dya);
        y1 = y1(N);
        y2 = y2(N);

        if y1*y2 > 0
            E2 = E2+10; %Se le suma 1 a E2 para buscar un valor válido para el que aplicar bisección
        end
        if y1*y2 < 0
            iteraciones = iteraciones + 1;
            E3 = (E1+E2)/2;
            y3 = numerov(k, E3, N, x, h, ya, dya);
            y3 = y3(N);
            if y1*y3 < 0
                E2 = E3;
                Elast = E2;
            end
            if y2*y3 < 0
                E1 = E3;
                Elast = E1;
            end
        end

    end
    E(i) = Elast;
    iteracionesv(i) = iteraciones;
    incertidumbres(i) = incertidumbre;
    y1 = 10000;
    y2 = 0;
end
end




