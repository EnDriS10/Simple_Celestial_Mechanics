clear; clc;

%% Datos

G_normal= 6.67E-11;

M_rr = 3.3011E23;
UA_m = 1.495978707E11;
M_s=30*8.64E4;

G = G_normal * (M_s^2) * M_rr^2 / (UA_m^3);

M_mercurio = 1;
D_mercurio = 0.387;
V_mercurio = 4.74E4;
%Incl_mercurio = (pi*7)/180;

M_marte = 1.944;
D_marte = 1.524;
V_marte = 2.413E4;
%Incl_marte = (pi*1.85)/180;

M_venus= 14.74;
D_venus= 0.723;
V_venus = 3.502E4;
%Incl_venus = (pi*3.394)/180;

M_tierra = 18.09;
D_tierra = 1;
V_tierra = 2.978E4;
%Incl_tierra = 0;

M_urano = 262.8;
D_urano = 19.229;
V_urano = 6.81E3;
%Incl_urano = (pi*0.774)/180;

M_neptuno = 3.1E2;
D_neptuno = 30.104;
V_neptuno = 5.43E3;
%Incl_neptuno = (pi*1.774)/180;

M_saturno = 1.722E3;
D_saturno = 9.537;
V_saturno = 9.67E3;
%Incl_saturno = (pi*2.488)/180;

M_jupiter= 5.752E3;
D_jupiter = 5.204;
V_jupiter = 1.307E4;
%Incl_jupiter = (pi*1.308)/180;

M_sol = 6.022E6;



%%

M= [M_mercurio, M_marte, M_venus, M_tierra, M_urano, M_neptuno, M_saturno, M_jupiter, M_sol];
D = [D_mercurio, D_marte, D_venus, D_tierra, D_urano, D_neptuno, D_saturno, D_jupiter, 0];
V = [V_mercurio, V_marte, V_venus, V_tierra, V_urano, V_neptuno, V_saturno, V_jupiter, 0]; 
%Incl = [Incl_mercurio, Incl_marte, Incl_venus, Incl_tierra, Incl_urano, Incl_neptuno, Incl_saturno, Incl_jupiter, 0];

%% Inicializar Posiciones

[r,v, p] = Inicializar(M, D, V);
    Tf=1E8;
    dt=1E-3;
[x,y] =Simul(M, r, v, p, Tf, dt);
%Graf(x,y,M)
%simularYAnimar(M, r, v, p, D, Tf, dt);

%%
plot(x(9,:),y(9,:))
 

%% Funciones
function [r_norm]= Norma(r)
   r_norm = sqrt( r(1,1)^2 + r(2,1)^2 );%+ r(3,1)^2 );
end

function [r,v,p] = Inicializar(M,D,V)

M_rr = 3.3011E23;
UA_m = 1.495978707E11;
M_s=30*8.64E4;
V_UA_d = V * M_s / UA_m;


r=zeros(2,length(D));
v=zeros(2,length(D));
p=zeros(2,length(D));


for i=1:length(D)

    Theta= 2*pi*rand(1);
    
    r(1,i) = D(i) * cos(Theta);
    r(2,i) = D(i) * sin(Theta);
    
    v(1,i) = -V_UA_d(i) * sin(Theta);
    v(2,i) = V_UA_d(i) * cos(Theta);

    p(:,i) = M(i).*v(:,i).*M_rr;
    
end
end


function [F] = ComputeF(M,r)
    G= 1.458602499152719e+16;    % G en UAs y multiplicada por el cuadrado de la masa de mercurio
    F=zeros(2,length(M));
        for i=1:length(M)
            for j=1:length(M)
                if i~=j
                    r_ij= r(:,i) - r(:,j);
                    d=Norma(r_ij);
                    F(:,i) = F(:,i) - G.*M(i).*M(j).*(r_ij)./(d)^3;
                end
            end
        end
end


%%

function [x,y] =Simul(M,r, v, p, Tf, dt)

    %Tf=1E5;
    %dt=1E-3;
    M_rr = 3.3011E23;
 
        
    x = zeros(9, Tf);
    y = zeros(9, Tf); 
 for step=1:Tf
    
    [F] = ComputeF(M, r);

    p = p + F.*dt;
    r = r + (p.*dt)./(M*M_rr);

    % Almacenar posiciones
    x(:, step) = r(1, :);
    y(:, step) = r(2, :);

 end
end

function fig = Graf(x,y,M)
figure;
hold on;
colors = lines(length(M));
for i = 1:length(M)
    
    if i == 1 
        plot(x(i,:), y(i,:), 'k', 'LineWidth', 2, 'DisplayName', 'Mercurio');    
    elseif i == 2 
        plot(x(i,:), y(i,:), 'r', 'LineWidth', 2, 'DisplayName', 'Marte');
    elseif i == 3 
        plot(x(i,:), y(i,:), 'y', 'LineWidth', 2, 'DisplayName', 'Venus');
    elseif i == 4  % Tierra
        plot(x(i,:), y(i,:), 'b', 'LineWidth', 2, 'DisplayName', 'Tierra');
    elseif i == 5 
        plot(x(i,:), y(i,:), 'c', 'LineWidth', 2, 'DisplayName', 'Urano');
    elseif i == 6  
        plot(x(i,:), y(i,:), 'Color', [0 0 0.5] , 'LineWidth', 2, 'DisplayName', 'Neptuno');
    elseif i == 7 
        plot(x(i,:), y(i,:), 'Color', [0 0 0.8] , 'LineWidth', 2, 'DisplayName', 'Saturno');
    elseif i == 8 
        plot(x(i,:), y(i,:), 'Color', [1 0.5 0] , 'LineWidth', 2, 'DisplayName', 'Jupiter');
    
    elseif i == 9  % Sol
        plot(x(i,:), y(i,:), 'yo', 'MarkerSize', 10, 'LineWidth', 3, 'DisplayName', 'Sol');
    end
end
xlabel('X (UA)');
ylabel('Y (UA)');
title('Órbitas Planetarias');
axis equal;
grid on;
legend('Location', 'best');
hold off;

end


%% Funcion principal de simulacion y animacion

function simularYAnimar(M, r, v, p, D, Tf, dt)
    M_rr = 3.3011E23;
    
    % Configurar animacion
    [x, y, planetas, trazas] = configurarAnimacion(r, M, D, Tf, dt);
    
    % Variables temporales para la simulacion
    r_current = r;
    p_current = p;
    
    % Bucle principal de simulacion
    for step = 1:Tf
        % Calcular fuerzas gravitacionales
        [F] = ComputeF(M, r_current);
        
        % Integrar ecuaciones de movimiento (metodo de Euler)
        p_current = p_current + F .* dt;
        r_current = r_current + (p_current .* dt) ./ (M * M_rr);
        
        % Almacenar posiciones
        x(:, step) = r_current(1, :);
        y(:, step) = r_current(2, :);
        
        % Actualizar animacion
        actualizarAnimacion(planetas, trazas, r_current, x, y, step, Tf);
    end
    
    fprintf('Simulacion completada: %d pasos\n', Tf);
end

function [x, y, planetas, trazas] = configurarAnimacion(r, M, D, Tf, dt)
    n = length(M);
    x = zeros(n, Tf);
    y = zeros(n, Tf);
    
    % Configurar figura
    figure('Name', 'Simulacion de Orbitas Planetarias', ...
           'NumberTitle', 'off', ...
           'Position', [100, 100, 800, 800]);
    hold on;
    grid on;
    axis equal;  % Mantener relacion de aspecto 1:1
    box on;
    xlabel('X (UA)', 'FontSize', 12);
    ylabel('Y (UA)', 'FontSize', 12);
    title('Animacion de Orbitas Planetarias', 'FontSize', 14);
    
    % Definir propiedades visuales
    nombres = {'Mercurio', 'Marte', 'Venus', 'Tierra', 'Urano', ...
               'Neptuno', 'Saturno', 'Jupiter', 'Sol'};
    
    % Para las lineas (trazas), usar solo colores, no marcadores
    colores_lineas = {'k', 'r', 'y', 'b', 'c', ...
                     [0 0 0.5], [0 0 0.8], [1 0.5 0], 'y'};
    
    % Para los marcadores de los planetas
    colores_marcadores = {'k', 'r', 'y', 'b', 'c', ...
                         [0 0 0.5], [0 0 0.8], [1 0.5 0], 'y'};
    
    Tam = [8, 8, 8, 8, 8, 8, 8, 8, 20];
    
    % Crear objetos graficos para planetas
    planetas = gobjects(n, 1);
    for i = 1:n
        if i == 9 % Sol
            planetas(i) = plot(r(1,i), r(2,i), 'yo', ...
                              'MarkerSize', Tam(i), ...
                              'MarkerFaceColor', 'y', ...
                              'LineWidth', 3, ...
                              'DisplayName', nombres{i});
        else
            planetas(i) = plot(r(1,i), r(2,i), 'o', ...
                              'Color', colores_marcadores{i}, ...
                              'MarkerSize', Tam(i), ...
                              'MarkerFaceColor', colores_marcadores{i}, ...
                              'LineWidth', 2, ...
                              'DisplayName', nombres{i});
        end
    end
    
    % Ajustar limites del grafico
    lim_max = max(D) * 1.2;
    xlim([-lim_max lim_max]);
    ylim([-lim_max lim_max]);
    
    % Anadir leyenda
    legend('Location', 'best', 'FontSize', 10);
    
    % Crear objetos para trazas (trayectorias)
    trazas = gobjects(n, 1);
    for i = 1:n
        trazas(i) = plot(NaN, NaN, '-', ...
                        'Color', colores_lineas{i}, ...
                        'LineWidth', 1, ...
                        'HandleVisibility', 'off');
    end
    
    % Anadir informacion adicional
    text(0.02, 0.98, sprintf('Pasos: %d\ndt: %.1e', Tf, dt), ...
         'Units', 'normalized', ...
         'VerticalAlignment', 'top', ...
         'BackgroundColor', 'w', ...
         'EdgeColor', 'k', ...
         'Margin', 2);
end

function actualizarAnimacion(planetas, trazas, r, x, y, step, Tf)
    n = length(planetas);
    
    % Actualizar posiciones de los planetas
    for i = 1:n
        set(planetas(i), 'XData', r(1,i), 'YData', r(2,i));
        
        % Actualizar trazas (ultimos 500 puntos)
        if step > 1
            inicio = max(1, step-500);
            set(trazas(i), 'XData', x(i, inicio:step), ...
                          'YData', y(i, inicio:step));
        end
    end
    

    % Actualizar titulo
    title(sprintf('Orbitas Planetarias - Paso: %d / %d (%.1f%%)', ...
                  step, Tf, 100*step/Tf), 'FontSize', 14);
    
    % Forzar actualizacion de la figura
    drawnow limitrate;
    
end