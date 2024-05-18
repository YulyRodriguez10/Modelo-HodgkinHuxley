function hodgkin_huxley
    
    C_m = 1.0;    % Capacitancia de la membrana, uF/cm^2
    g_Na = 120.0; % Conductancia máxima de Na, mS/cm^2
    g_K = 36.0;   % Conductancia máxima de K, mS/cm^2
    g_L = 0.3;    % Conductancia de fuga, mS/cm^2
    E_Na = 50.0;  % Potencial de equilibrio de Na, mV
    E_K = -77.0;  % Potencial de equilibrio de K, mV
    E_L = -54.4;  % Potencial de equilibrio de fuga, mV

   
    t_end = 50; % Tiempo de simulación, ms
    dt = 0.01;  % Paso de tiempo, ms
    time = 0:dt:t_end; % Vector de tiempo

    
    V = -65; % Potencial de membrana, mV
    m = alpha_m(V) / (alpha_m(V) + beta_m(V));
    h = alpha_h(V) / (alpha_h(V) + beta_h(V));
    n = alpha_n(V) / (alpha_n(V) + beta_n(V));

    
    I_ext = zeros(size(time));
    I_ext(time >= 10 & time <= 40) = 10; % Corriente externa de 10 uA/cm^2 entre 10 ms y 40 ms

    
    V_trace = zeros(size(time));
    m_trace = zeros(size(time));
    h_trace = zeros(size(time));
    n_trace = zeros(size(time));

   
    for i = 1:length(time)
        % Corrientes iónicas
        I_Na = g_Na * m^3 * h * (V - E_Na);
        I_K = g_K * n^4 * (V - E_K);
        I_L = g_L * (V - E_L);

        
        dVdt = (I_ext(i) - I_Na - I_K - I_L) / C_m;

        
        V = V + dt * dVdt;

        
        dm = alpha_m(V) * (1 - m) - beta_m(V) * m;
        dh = alpha_h(V) * (1 - h) - beta_h(V) * h;
        dn = alpha_n(V) * (1 - n) - beta_n(V) * n;

        
        m = m + dt * dm;
        h = h + dt * dh;
        n = n + dt * dn;

        
        V_trace(i) = V;
        m_trace(i) = m;
        h_trace(i) = h;
        n_trace(i) = n;
    end

    
    figure;
    subplot(2,1,1);
    plot(time, V_trace);
    title('Potencial de membrana');
    xlabel('Tiempo (ms)');
    ylabel('V (mV)');

    subplot(2,1,2);
    plot(time, I_ext);
    title('Corriente de estímulo');
    xlabel('Tiempo (ms)');
    ylabel('I_{ext} (\muA/cm^2)');
end


function val = alpha_m(V)
    val = 0.1 * (V + 40) / (1 - exp(-(V + 40) / 10));
end

function val = beta_m(V)
    val = 4 * exp(-(V + 65) / 18);
end

function val = alpha_h(V)
    val = 0.07 * exp(-(V + 65) / 20);
end

function val = beta_h(V)
    val = 1 / (1 + exp(-(V + 35) / 10));
end

function val = alpha_n(V)
    val = 0.01 * (V + 55) / (1 - exp(-(V + 55) / 10));
end

function val = beta_n(V)
    val = 0.125 * exp(-(V + 65) / 80);
end

