clear; clc; close all;

%% ===============================
%  SIMULAZIONE COMPARATIVA CASI
%  Navier-Stokes vs Termo-Vorticale
%  ===============================

Tf = 200;
dt = 0.01;
tspan = 0:dt:Tf;

% Stato ridotto: X = [s0; s1c; theta1c; zeta1s]
X0 = [0.1; 0.05; 0.08; 0.02];

%% Matrice lineare del modello ridotto
A = [
    -0.08   0.15   0.00   0.05;
    -0.12   0.04   0.20   0.00;
     0.00  -0.25  -0.10   0.12;
     0.10   0.00  -0.18   0.03
];

%% Lista casi da simulare

cases(1).name = 'NS puro';
cases(1).eta = 0;
cases(1).mu = 0;
cases(1).c_omega = 1;
cases(1).A_geom = 1;
cases(1).delta = 0.02;

cases(2).name = 'Termo-vorticale puro';
cases(2).eta = 0.5;
cases(2).mu = 0;
cases(2).c_omega = 1;
cases(2).A_geom = 1;
cases(2).delta = 0.02;

cases(3).name = 'Memoria pura';
cases(3).eta = 0;
cases(3).mu = 0.5;
cases(3).c_omega = 1;
cases(3).A_geom = 1;
cases(3).delta = 0.02;

cases(4).name = 'Modello completo';
cases(4).eta = 0.5;
cases(4).mu = 0.5;
cases(4).c_omega = 1;
cases(4).A_geom = 1;
cases(4).delta = 0.02;

cases(5).name = 'TV senza saturazione';
cases(5).eta = 0.5;
cases(5).mu = 0;
cases(5).c_omega = 0;
cases(5).A_geom = 1;
cases(5).delta = 0.02;

cases(6).name = 'TV anisotropo';
cases(6).eta = 0.5;
cases(6).mu = 0.5;
cases(6).c_omega = 1;
cases(6).A_geom = 1.4;
cases(6).delta = 0.02;

%% Simulazione

results = struct();

for k = 1:length(cases)

    par = cases(k);

    fprintf('\n=====================================\n');
    fprintf('Simulazione caso: %s\n', par.name);
    fprintf('eta = %.3f, mu = %.3f, c_omega = %.3f, A_geom = %.3f\n', ...
        par.eta, par.mu, par.c_omega, par.A_geom);
    fprintf('=====================================\n');

    f = @(t,X) modello_ridotto_TV(t, X, A, par);

    opts = odeset('RelTol',1e-9,'AbsTol',1e-11);

    [t, X] = ode45(f, tspan, X0, opts);

    s0     = X(:,1);
    s1c    = X(:,2);
    theta  = X(:,3);
    zeta   = X(:,4);

    % Indicatori
    E = 0.5 * sum(X.^2, 2);
    Omega = sqrt(s1c.^2 + zeta.^2);
    R = sqrt(sum(X.^2, 2));

    results(k).name = par.name;
    results(k).t = t;
    results(k).X = X;
    results(k).E = E;
    results(k).Omega = Omega;
    results(k).R = R;

    results(k).E_mean = mean(E);
    results(k).E_max = max(E);
    results(k).Omega_mean = mean(Omega);
    results(k).Omega_max = max(Omega);
    results(k).Omega_std = std(Omega);
    results(k).R_mean = mean(R);
    results(k).R_max = max(R);

end

%% Tabella risultati

names = string({results.name})';

E_mean = [results.E_mean]';
E_max = [results.E_max]';
Omega_mean = [results.Omega_mean]';
Omega_max = [results.Omega_max]';
Omega_std = [results.Omega_std]';
R_mean = [results.R_mean]';
R_max = [results.R_max]';

T = table(names, E_mean, E_max, Omega_mean, Omega_max, Omega_std, R_mean, R_max);

disp(' ');
disp('========= TABELLA COMPARATIVA =========');
disp(T);

%% Grafico energia

figure;
hold on; grid on;
for k = 1:length(results)
    plot(results(k).t, results(k).E, 'LineWidth', 1.2);
end
xlabel('t');
ylabel('E(t)');
title('Confronto energia nei diversi casi');
legend({results.name}, 'Interpreter','none');

%% Grafico vorticità efficace

figure;
hold on; grid on;
for k = 1:length(results)
    plot(results(k).t, results(k).Omega, 'LineWidth', 1.2);
end
xlabel('t');
ylabel('\Omega(t)');
title('Confronto vorticità efficace');
legend({results.name}, 'Interpreter','none');

%% Spazio fase Omega - dOmega/dt

figure;
hold on; grid on;
for k = 1:length(results)
    Om = results(k).Omega;
    dOm = gradient(Om, results(k).t);
    plot(Om, dOm, 'LineWidth', 1.1);
end
xlabel('\Omega');
ylabel('d\Omega/dt');
title('Spazio fase della vorticità efficace');
legend({results.name}, 'Interpreter','none');

%% Traiettorie 3D ridotte

figure;
hold on; grid on;
for k = 1:length(results)
    X = results(k).X;
    plot3(X(:,2), X(:,3), X(:,4), 'LineWidth', 1.1);
end
xlabel('s_{1c}');
ylabel('\theta_{1c}');
zlabel('\zeta_{1s}');
title('Traiettorie nello spazio delle fasi ridotto');
legend({results.name}, 'Interpreter','none');
view(45,30);

%% Spettro della vorticità efficace

figure;
hold on; grid on;

for k = 1:length(results)

    Om = results(k).Omega;
    Om = Om - mean(Om);

    N = length(Om);
    Fs = 1/dt;

    Y = fft(Om);
    P2 = abs(Y/N);
    P1 = P2(1:floor(N/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);

    f = Fs*(0:floor(N/2))/N;

    plot(f, P1, 'LineWidth', 1.1);

end

xlabel('Frequenza');
ylabel('|FFT(\Omega)|');
title('Spettro della vorticità efficace');
xlim([0 5]);
legend({results.name}, 'Interpreter','none');

%% Salvataggio

save('risultati_confronto_TV_NS.mat', 'results', 'T', 'cases');

fprintf('\nSimulazione completata.\n');
fprintf('File salvato: risultati_confronto_TV_NS.mat\n');


%% ============================================================
%  FUNZIONE MODELLO RIDOTTO
%  ============================================================

function dXdt = modello_ridotto_TV(~, X, A, par)

    s0    = X(1);
    s1c   = X(2);
    theta = X(3);
    zeta  = X(4);

    eta = par.eta;
    mu = par.mu;
    c_omega = par.c_omega;
    A_geom = par.A_geom;
    delta = par.delta;

    % Parte lineare
    L = A * X;

    % Vorticità ridotta efficace
    omega_eff = sqrt(s1c^2 + zeta^2);

    % Saturazione termo-vorticale
    sat = 1 + c_omega * omega_eff^2;

    % Se c_omega = 0, sat = 1: caso senza saturazione
    if sat == 0
        sat = 1;
    end

    % Termine termo-vorticale ridotto
    TV = eta / sat * [
        0;
        -A_geom * theta * zeta;
         A_geom * s1c * zeta;
        -A_geom * theta * s1c
    ];

    % Termine memoria ridotto
    MEM = mu * [
        -0.10 * zeta;
         0.08 * s0;
        -0.12 * zeta;
         0.10 * theta
    ];

    % Saturazione cubica globale
    NONLIN = -delta * norm(X)^2 * X;

    dXdt = L + TV + MEM + NONLIN;

end