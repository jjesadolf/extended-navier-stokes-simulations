clear; clc; close all;

%% ============================================================
%  TV_NS_PDE_SigmaEff_Diagnostic_Stable_v2.m
%  Diagnostica stabile di sigma_eff per modello termo-vorticale
%  con proiezione incomprimibile corretta.
%  ============================================================

N = 24;
L = 2*pi;

x = linspace(0,L,N+1);
x(end) = [];

[X,Y,Z] = meshgrid(x,x,x);

nu = 0.08;
kappa = 0.08;

eta = 0.5;
c_omega = 10.0;

dt = 1e-4;
Tf = 2.0;
Nt = round(Tf/dt);

diagnostic_every = 200;

%% Griglia spettrale

k = [0:N/2 -N/2+1:-1]*(2*pi/L);
[KX,KY,KZ] = meshgrid(k,k,k);

K2 = KX.^2 + KY.^2 + KZ.^2;
K2(1,1,1) = 1;

%% Filtro spettrale anti-aliasing / stabilizzazione

Kmax = max(abs(k));
Kabs = sqrt(K2);

filter = exp(-36*(Kabs/(0.65*Kmax)).^12);
filter(1,1,1) = 1;

%% Condizioni iniziali moderate

u = 0.15*sin(Y) + 0.08*cos(Z);
v = 0.15*sin(Z) + 0.08*cos(X);
w = 0.15*sin(X) + 0.08*cos(Y);

[u,v,w] = project_incompressible(u,v,w,KX,KY,KZ,K2);

T = 0.2*sin(X).*cos(Y) + 0.1*cos(Z);

%% Storici

time_hist = [];
omega_max_hist = [];
omega_mean_hist = [];
sigma_mean_high_hist = [];
sigma_eff_mean_high_hist = [];
frac_control_hist = [];
S_omega_hist = [];
I_TV_hist = [];
R_control_hist = [];

fprintf('\n=============================================\n');
fprintf('  STABLE TV-NS PDE SIGMA_EFF DIAGNOSTIC v2\n');
fprintf('  eta = %.3f, c_omega = %.3f\n', eta, c_omega);
fprintf('=============================================\n');

%% Ciclo temporale

for n = 1:Nt

    t = n*dt;

    %% Derivate velocità

    [ux,uy,uz] = grad3(u,KX,KY,KZ);
    [vx,vy,vz] = grad3(v,KX,KY,KZ);
    [wx,wy,wz] = grad3(w,KX,KY,KZ);

    %% Vorticità omega = curl u

    omega_x = wy - vz;
    omega_y = uz - wx;
    omega_z = vx - uy;

    omega_abs = sqrt(omega_x.^2 + omega_y.^2 + omega_z.^2);

    %% Gradiente temperatura

    [Tx,Ty,Tz] = grad3(T,KX,KY,KZ);

    %% Termine termo-vorticale
    % F = -eta (grad T x omega)/(1+c_omega |omega|^2)

    [cx,cy,cz] = cross_field(Tx,Ty,Tz,omega_x,omega_y,omega_z);

    denom = 1 + c_omega*omega_abs.^2;

    Fx = -eta*cx./denom;
    Fy = -eta*cy./denom;
    Fz = -eta*cz./denom;

    %% Termine convettivo

    adv_u = u.*ux + v.*uy + w.*uz;
    adv_v = u.*vx + v.*vy + w.*vz;
    adv_w = u.*wx + v.*wy + w.*wz;

    rhs_u = -adv_u + Fx;
    rhs_v = -adv_v + Fy;
    rhs_w = -adv_w + Fz;

    %% Proiezione incomprimibile del RHS

    [rhs_u,rhs_v,rhs_w] = project_incompressible(rhs_u,rhs_v,rhs_w,KX,KY,KZ,K2);

    %% Passo semi-implicito viscoso

    uhat = fftn(u);
    vhat = fftn(v);
    what = fftn(w);

    rhs_uhat = fftn(rhs_u).*filter;
    rhs_vhat = fftn(rhs_v).*filter;
    rhs_what = fftn(rhs_w).*filter;

    uhat_new = (uhat + dt*rhs_uhat)./(1 + dt*nu*K2);
    vhat_new = (vhat + dt*rhs_vhat)./(1 + dt*nu*K2);
    what_new = (what + dt*rhs_what)./(1 + dt*nu*K2);

    u = real(ifftn(uhat_new));
    v = real(ifftn(vhat_new));
    w = real(ifftn(what_new));

    [u,v,w] = project_incompressible(u,v,w,KX,KY,KZ,K2);

    %% Evoluzione termica: advezione-diffusione

    [Tx,Ty,Tz] = grad3(T,KX,KY,KZ);

    adv_T = u.*Tx + v.*Ty + w.*Tz;

    That = fftn(T);
    rhs_That = fftn(-adv_T).*filter;

    That_new = (That + dt*rhs_That)./(1 + dt*kappa*K2);

    T = real(ifftn(That_new));

    %% Controllo NaN / Inf

    if any(~isfinite(u(:))) || any(~isfinite(v(:))) || ...
       any(~isfinite(w(:))) || any(~isfinite(T(:)))

        error('NaN o Inf generato al tempo t = %.6f. Riduci dt, eta o ampiezza iniziale.', t);

    end

    %% Diagnostica

    if mod(n,diagnostic_every)==0

        diag = compute_sigma_eff_diagnostic( ...
            u,v,w,T,eta,c_omega,KX,KY,KZ,K2);

        time_hist(end+1,1) = t;
        omega_max_hist(end+1,1) = diag.omega_max;
        omega_mean_hist(end+1,1) = diag.omega_mean;
        sigma_mean_high_hist(end+1,1) = diag.sigma_mean_high;
        sigma_eff_mean_high_hist(end+1,1) = diag.sigma_eff_mean_high;
        frac_control_hist(end+1,1) = diag.frac_control_high;
        S_omega_hist(end+1,1) = diag.S_omega;
        I_TV_hist(end+1,1) = diag.I_TV;
        R_control_hist(end+1,1) = diag.R_control;

        fprintf('t = %.4f | omega_max = %.5f | sigma_eff_high = %.5e | control = %.2f %%\n', ...
            t, diag.omega_max, diag.sigma_eff_mean_high, 100*diag.frac_control_high);

    end
end

%% Tabella finale

FinalTable = table( ...
    time_hist, ...
    omega_max_hist, ...
    omega_mean_hist, ...
    sigma_mean_high_hist, ...
    sigma_eff_mean_high_hist, ...
    frac_control_hist, ...
    S_omega_hist, ...
    I_TV_hist, ...
    R_control_hist);

disp(' ');
disp('========= DIAGNOSTICA SIGMA_EFF STABILE v2 =========');
disp(FinalTable);

%% Grafici

figure;
plot(time_hist, omega_max_hist, 'LineWidth', 1.4);
grid on;
xlabel('t');
ylabel('max |\omega|');
title('Massimo della vorticità');

figure;
plot(time_hist, sigma_mean_high_hist, 'LineWidth', 1.4);
hold on;
plot(time_hist, sigma_eff_mean_high_hist, 'LineWidth', 1.4);
grid on;
xlabel('t');
ylabel('media su zone ad alta vorticità');
title('\sigma vs \sigma_{eff}');
legend('\sigma','\sigma_{eff}','Location','best');

figure;
plot(time_hist, frac_control_hist, 'LineWidth', 1.4);
grid on;
xlabel('t');
ylabel('frazione');
title('Frazione zone critiche con \sigma_{eff}<0');
ylim([0 1]);

figure;
plot(time_hist, S_omega_hist, 'LineWidth', 1.4);
hold on;
plot(time_hist, eta*I_TV_hist, 'LineWidth', 1.4);
grid on;
xlabel('t');
ylabel('valore medio');
title('Stretching classico vs termine termo-vorticale');
legend('S_\omega','\eta I_{TV}','Location','best');

figure;
plot(time_hist, R_control_hist, 'LineWidth', 1.4);
grid on;
xlabel('t');
ylabel('R_{control}');
title('Rapporto di controllo termo-vorticale');

%% Salvataggio

save('risultati_sigma_eff_PDE_stable_v2.mat', ...
    'FinalTable', ...
    'time_hist', ...
    'omega_max_hist', ...
    'omega_mean_hist', ...
    'sigma_mean_high_hist', ...
    'sigma_eff_mean_high_hist', ...
    'frac_control_hist', ...
    'S_omega_hist', ...
    'I_TV_hist', ...
    'R_control_hist');

fprintf('\nSimulazione completata.\n');
fprintf('File salvato: risultati_sigma_eff_PDE_stable_v2.mat\n');

%% ============================================================
%  FUNZIONI LOCALI
%  ============================================================

function [fx,fy,fz] = grad3(f,KX,KY,KZ)

    fhat = fftn(f);

    fx = real(ifftn(1i*KX.*fhat));
    fy = real(ifftn(1i*KY.*fhat));
    fz = real(ifftn(1i*KZ.*fhat));

end

function [u,v,w] = project_incompressible(u,v,w,KX,KY,KZ,K2)

    uhat = fftn(u);
    vhat = fftn(v);
    what = fftn(w);

    % Proiezione di Leray:
    % uhat_perp = uhat - k (k·uhat)/|k|^2

    k_dot_uhat = KX.*uhat + KY.*vhat + KZ.*what;

    qhat = k_dot_uhat ./ K2;

    uhat = uhat - KX.*qhat;
    vhat = vhat - KY.*qhat;
    what = what - KZ.*qhat;

    uhat(1,1,1) = 0;
    vhat(1,1,1) = 0;
    what(1,1,1) = 0;

    u = real(ifftn(uhat));
    v = real(ifftn(vhat));
    w = real(ifftn(what));

end

function [cx,cy,cz] = cross_field(ax,ay,az,bx,by,bz)

    cx = ay.*bz - az.*by;
    cy = az.*bx - ax.*bz;
    cz = ax.*by - ay.*bx;

end

function diag = compute_sigma_eff_diagnostic(u,v,w,T,eta,c_omega,KX,KY,KZ,K2)

    %% Derivate velocità

    [ux,uy,uz] = grad3(u,KX,KY,KZ);
    [vx,vy,vz] = grad3(v,KX,KY,KZ);
    [wx,wy,wz] = grad3(w,KX,KY,KZ);

    %% Vorticità

    omega_x = wy - vz;
    omega_y = uz - wx;
    omega_z = vx - uy;

    omega_abs = sqrt(omega_x.^2 + omega_y.^2 + omega_z.^2);
    omega_safe = omega_abs + 1e-12;

    xi_x = omega_x ./ omega_safe;
    xi_y = omega_y ./ omega_safe;
    xi_z = omega_z ./ omega_safe;

    %% Tensore di deformazione S

    S11 = ux;
    S22 = vy;
    S33 = wz;

    S12 = 0.5*(uy + vx);
    S13 = 0.5*(uz + wx);
    S23 = 0.5*(vz + wy);

    %% sigma = xi · S xi

    sigma = ...
        S11.*xi_x.^2 + ...
        S22.*xi_y.^2 + ...
        S33.*xi_z.^2 + ...
        2*S12.*xi_x.*xi_y + ...
        2*S13.*xi_x.*xi_z + ...
        2*S23.*xi_y.*xi_z;

    %% Gradiente temperatura

    [Tx,Ty,Tz] = grad3(T,KX,KY,KZ);

    %% curl omega

    [ox_x,ox_y,ox_z] = grad3(omega_x,KX,KY,KZ);
    [oy_x,oy_y,oy_z] = grad3(omega_y,KX,KY,KZ);
    [oz_x,oz_y,oz_z] = grad3(omega_z,KX,KY,KZ);

    curl_omega_x = oz_y - oy_z;
    curl_omega_y = ox_z - oz_x;
    curl_omega_z = oy_x - ox_y;

    %% gradT x omega

    [gtxo_x,gtxo_y,gtxo_z] = cross_field( ...
        Tx,Ty,Tz,omega_x,omega_y,omega_z);

    %% Theta_T

    numerator = ...
        curl_omega_x.*gtxo_x + ...
        curl_omega_y.*gtxo_y + ...
        curl_omega_z.*gtxo_z;

    Theta_T = numerator ./ (omega_abs.^2 + 1e-12);

    %% sigma_eff

    sigma_eff = sigma - eta*Theta_T ./ (1 + c_omega*omega_abs.^2);

    %% Regione critica: alto 20% della vorticità

    Kcrit = prctile(omega_abs(:),80);
    high = omega_abs >= Kcrit;

    if nnz(high) == 0
        high = omega_abs >= max(omega_abs(:));
    end

    sigma_mean_high = mean(sigma(high));
    sigma_eff_mean_high = mean(sigma_eff(high));

    frac_control_high = sum(sigma_eff(high)<0)/nnz(high);

    %% Integrali medi globali

    S_omega_density = omega_abs.^2 .* sigma;
    I_TV_density = numerator ./ (1 + c_omega*omega_abs.^2);

    S_omega = mean(S_omega_density(:));
    I_TV = mean(I_TV_density(:));

    R_control = eta*I_TV/(abs(S_omega)+1e-12);

    %% Output

    diag.omega_max = max(omega_abs(:));
    diag.omega_mean = mean(omega_abs(:));
    diag.sigma_mean_high = sigma_mean_high;
    diag.sigma_eff_mean_high = sigma_eff_mean_high;
    diag.frac_control_high = frac_control_high;
    diag.S_omega = S_omega;
    diag.I_TV = I_TV;
    diag.R_control = R_control;

end