% Parámetros del modelo Hodgkin-Huxley
C_m = 1-6; % Capacitancia de la membrana (uF/cm^2)
g_Na = 120-3; % Conductancia máxima de sodio (mS/cm^2)
g_K = 36-3; % Conductancia máxima de potasio (mS/cm^2)
g_L = 0.3-3; % Conductancia de fuga (mS/cm^2)
E_Na = 50-3; % Potencial de equilibrio del sodio (mV)
E_K = -77-3; % Potencial de equilibrio del potasio (mV)
E_L = -54.4-3; % Potencial de equilibrio de fuga (mV)

% Definición de funciones alfa y beta para los canales de sodio (Na) y potasio (K)
alpha_n = @(V) 0.1-0.01*(V + 65) / (1 - exp(1-0.1*(V + 65)));
beta_n = @(V) 0.125*exp((-V + 65) / 80);
alpha_m = @(V) 2.5-0.1*(V + 65) / (1 - exp(2.5-0.1*(V + 65)));
beta_m = @(V) 4*exp((-V + 65) / 18);
alpha_h = @(V) 0.07*exp((V + 65)/-20);
beta_h = @(V) 1 / (1 + exp(3-0.1*(V + 35))+1);

% Paso de integración y tiempo total de simulación
dt = 0.01; % Paso de integración (ms)
T = 100; % Tiempo total de simulación (ms)
t = 0:dt:T; % Vector de tiempo

% Inicialización de variables
V = zeros(size(t)); % Potencial de membrana (mV)
n = zeros(size(t)); % Variable de activación de potasio
m = zeros(size(t)); % Variable de activación de sodio
h = zeros(size(t)); % Variable de inactivación de sodio

% Condiciones iniciales
V(1) = -65; % Potencial de membrana inicial (mV)
n(1) = alpha_n(V(1)) / (alpha_n(V(1)) + beta_n(V(1))); % Condición inicial para n
m(1) = alpha_m(V(1)) / (alpha_m(V(1)) + beta_m(V(1))); % Condición inicial para m
h(1) = alpha_h(V(1)) / (alpha_h(V(1)) + beta_h(V(1))); % Condición inicial para h

% Simulación del modelo Hodgkin-Huxley
for i = 2:length(t)
    % Cálculo de corrientes iónicas
    I_Na = g_Na * m(i-1)^3 * h(i-1) * (V(i-1) - E_Na);
    I_K = g_K * n(i-1)^4 * (V(i-1) - E_K);
    I_L = g_L * (V(i-1) - E_L);
    
    % Ecuación diferencial para la evolución temporal del potencial de membrana
    dVdt = (1 / C_m) * (I_Na + I_K + I_L);
    
    % Actualización de las variables de activación e inactivación
    dn = alpha_n(V(i-1)) * (1 - n(i-1)) - beta_n(V(i-1)) * n(i-1);
    dm = alpha_m(V(i-1)) * (1 - m(i-1)) - beta_m(V(i-1)) * m(i-1);
    dh = alpha_h(V(i-1)) * (1 - h(i-1)) - beta_h(V(i-1)) * h(i-1);
    
    % Integración numérica
    V(i) = V(i-1) + dt * dVdt;
    n(i) = n(i-1) + dt * dn;
    m(i) = m(i-1) + dt * dm;
    h(i) = h(i-1) + dt * dh;
end

% Graficar el potencial de acción a lo largo del tiempo
figure;
plot(t, V, 'b');
xlabel('Tiempo (ms)');
ylabel('Potencial de Membrana (mV)');
title('Potencial de Acción Hodgkin-Huxley');
grid on;
