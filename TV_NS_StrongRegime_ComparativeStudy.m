clear; clc; close all;

%% ============================================================
%  STRONG-REGIME COMPARATIVE STUDY
%  Thermo-Vortical Reduced Extension of Navier-Stokes
%  File: TV_NS_StrongRegime_ComparativeStudy.m
%  ============================================================

Tf = 200;
dt = 0.01;
tspan = 0:dt:Tf;

%% Condizione iniziale forte
X0 = [1.0; 0.8; 1.2; 0.6];

%% Matrice lineare più instabile
A = [
    -0.04   0.25   0.00   0.08;
    -0.10   0.12   0.30   0.00;
     0.00  -0.35  -0.04   0.18;
     0.16   0.00  -0.22   0.10
];

%% Casi da confrontare

cases(1).name = 'NS puro';
cases(1).eta = 0;
cases(1).mu = 0;
cases(1).c_omega = 1;
cases(1).A_geom = 1;
cases(1).delta = 0.02;

cases(2).name = 'TV puro saturato';
cases(2).eta = 2.0;
cases(2).mu = 0;
cases(2).c_omega = 1;
cases(2).A_geom = 1;
cases(2).delta = 0.02;

cases(3).name = 'TV puro senza saturazione';
cases(3).eta = 2.0;
cases(3).mu = 0;
cases(3).c_omega = 0;
cases(3).A_geom = 1;
cases(3).delta = 0.02;

cases(4).name = 'Memoria pura';
cases(4).eta = 0;
cases(4).mu = 1.5;
cases(4).c_omega = 1;
cases(4).A_geom = 1;
cases(4).delta = 0.02;

cases(5).name = 'Modello completo';
cases(5).eta = 2.0;
cases(5).mu = 1.5;
cases(5).c_omega = 1;
cases(5).A_geom = 1;
cases(5).delta = 0.02;

cases(6).name = 'Completo anisotropo';
cases(6).eta = 2.0;
cases(6).mu = 1.5;
cases(6).c_omega = 1;
cases(6).A_geom = 1.4;
cases(6).delta = 0.02;

cases(7).name = 'Forte saturazione';
cases(7).eta = 2.0;
cases(7).mu = 1.5;
cases(7).c_omega = 10;
cases(7).A_geom = 1;
cases(7).delta = 0.02;

cases(8).name = 'Feedback destabilizzante';
cases(8).eta = -2.0;
cases(8).mu = 0;
cases(8).c_omega = 1;
cases(8).A_geom = 1;
cases(8).delta = 0.02;

%% Simulazioni

results = struct();

for k = 1:length(cases)

    par = cases(k);

    fprintf('\n=====================================\n');
    fprintf('Simulazione caso: %s\n', par.name);
    fprintf('eta = %.3f, mu = %.3f, c_omega = %.3f, A_geom = %.3f\n', ...
        par.eta, par.mu, par.c_omega, par.A_geom);
    fprintf('=====================================\n');

    f = @(t,X) modello_ridotto_TV_strong(t, X, A, par);

    opts = odeset( ...
        'RelTol', 1e-8, ...
        'AbsTol', 1e-10, ...
        'MaxStep', 0.05);

    [t, X] = ode45(f, tspan, X0, opts);

    s0    = X(:,1);
    s1c   = X(:,2);
    theta = X(:,3);
    zeta  = X(:,4);

    E = 0.5 * sum(X.^2, 2);
    Omega = sqrt(s1c.^2 + zeta.^2);
    R = sqrt(sum(X.^2, 2));

    dOmega = gradient(Omega, t);
    dE = gradient(E, t);

    results(k).name = par.name;
    results(k).t = t;
    results(k).X = X;
    results(k).E = E;
    results(k).Omega = Omega;
    results(k).R = R;
    results(k).dOmega = dOmega;
    results(k).dE = dE;

    results(k).E_mean = mean(E);
    results(k).E_max = max(E);
    results(k).E_final = E(end);

    results(k).Omega_mean = mean(Omega);
    results(k).Omega_max = max(Omega);
    results(k).Omega_final = Omega(end);
    results(k).Omega_std = std(Omega);

    results(k).R_mean = mean(R);
    results(k).R_max = max(R);
    results(k).R_final = R(end);

end

%% Tabella comparativa

names = string({results.name})';

E_mean = [results.E_mean]';
E_max = [results.E_max]';
E_final = [results.E_final]';

Omega_mean = [results.Omega_mean]';
Omega_max = [results.Omega_max]';
Omega_final = [results.Omega_final]';
Omega_std = [results.Omega_std]';

R_mean = [results.R_mean]';
R_max = [results.R_max]';
R_final = [results.R_final]';

T = table( ...
    names, ...
    E_mean, E_max, E_final, ...
    Omega_mean, Omega_max, Omega_final, Omega_std, ...
    R_mean, R_max, R_final);

disp(' ');
disp('========= TABELLA COMPARATIVA REGIME FORTE =========');
disp(T);

%% Figura 1: energia

figure;
hold on; grid on;
for k = 1:length(results)
    plot(results(k).t, results(k).E, 'LineWidth', 1.2);
end
xlabel('t');
ylabel('E(t)');
title('Regime forte: confronto energia');
legend({results.name}, 'Interpreter','none', 'Location','best');

%% Figura 2: vorticità efficace

figure;
hold on; grid on;
for k = 1:length(results)
    plot(results(k).t, results(k).Omega, 'LineWidth', 1.2);
end
xlabel('t');
ylabel('\Omega(t)');
title('Regime forte: confronto vorticità efficace');
legend({results.name}, 'Interpreter','none', 'Location','best');

%% Figura 3: spazio fase Omega-dOmega

figure;
hold on; grid on;
for k = 1:length(results)
    plot(results(k).Omega, results(k).dOmega, 'LineWidth', 1.1);
end
xlabel('\Omega');
ylabel('d\Omega/dt');
title('Regime forte: spazio fase della vorticità');
legend({results.name}, 'Interpreter','none', 'Location','best');

%% Figura 4: traiettorie 3D

figure;
hold on; grid on;
for k = 1:length(results)
    X = results(k).X;
    plot3(X(:,2), X(:,3), X(:,4), 'LineWidth', 1.1);
end
xlabel('s_{1c}');
ylabel('\theta_{1c}');
zlabel('\zeta_{1s}');
title('Regime forte: traiettorie nello spazio ridotto');
legend({results.name}, 'Interpreter','none', 'Location','best');
view(45,30);

%% Figura 5: spettro FFT vorticità

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
title('Regime forte: spettro della vorticità efficace');
xlim([0 5]);
legend({results.name}, 'Interpreter','none', 'Location','best');

%% Figura 6: energia vs vorticità

figure;
hold on; grid on;
for k = 1:length(results)
    plot(results(k).E, results(k).Omega, 'LineWidth', 1.1);
end
xlabel('E(t)');
ylabel('\Omega(t)');
title('Regime forte: diagramma energia-vorticità');
legend({results.name}, 'Interpreter','none', 'Location','best');

%% Figura 7: confronto normalizzato Omega

figure;
hold on; grid on;
for k = 1:length(results)
    Om = results(k).Omega;
    Omn = Om / max(Om);
    plot(results(k).t, Omn, 'LineWidth', 1.2);
end
xlabel('t');
ylabel('\Omega(t)/\Omega_{max}');
title('Regime forte: vorticità normalizzata');
legend({results.name}, 'Interpreter','none', 'Location','best');

%% Salvataggio

save('risultati_TV_NS_regime_forte.mat', 'results', 'T', 'cases', 'A', 'X0');

fprintf('\nSimulazione completata.\n');
fprintf('File salvato: risultati_TV_NS_regime_forte.mat\n');

%% ============================================================
%  FUNZIONE DEL MODELLO RIDOTTO IN REGIME FORTE
%  ============================================================

function dXdt = modello_ridotto_TV_strong(~, X, A, par)

    s0    = X(1);
    s1c   = X(2);
    theta = X(3);
    zeta  = X(4);

    eta = par.eta;
    mu = par.mu;
    c_omega = par.c_omega;
    A_geom = par.A_geom;
    delta = par.delta;

    L = A * X;

    omega_eff = sqrt(s1c^2 + zeta^2);

    sat = 1 + c_omega * omega_eff^2;

    TV = eta / sat * [
        0;
        -A_geom * theta * zeta;
         A_geom * s1c * zeta;
        -A_geom * theta * s1c
    ];

    MEM = mu * [
        -0.10 * zeta;
         0.08 * s0;
        -0.12 * zeta;
         0.10 * theta
    ];

    NONLIN = -delta * norm(X)^2 * X;

    dXdt = L + TV + MEM + NONLIN;

end