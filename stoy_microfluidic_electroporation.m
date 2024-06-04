%% Transfer Function Setup
should_rerun_tfest = 0;
should_rerun_tfest2 = 0;
should_rerun_tfest3 = 0;
should_rerun_tfest4 = 0;

%% Constants
U_0 = 100; % V
z_l = 200e-6; % m
z_d = 50e-9; % m
ep_l = 80;
ep_d = 80;
ep_0 = 8.8541878188E-12; % F m-1
e_l = U_0 / (z_l + 2 * z_d * ep_l / ep_d); % V m-1
sigma_l_mscm = 13.4; %mS cm-1 %500E-6; % S m-1

ep_m = 2;
ep_l = 80;
ep_c = 80;
R1 = 5E-6;
d = 5E-9;
R2 = R1 + d;
sigma_c = 0.5;
R1R2_RATIO = R1 / R2;

% R - cell radius
R = R1;
% C - membrane capacitance
C = 1.0e-2;
% lo - lambda o - extracellular medium conductivity
lo = sigma_l; % 2.0e-1; 
% lm - lambda m - membrane conductivity
lm = 5.0e-7;
% li - cytoplasmic conductivity
li = sigma_c; %2.0e-1;

%% Unit Conversions
J_TO_MJ = 1E3;
M_TO_CM = 1E2;
V_TO_KV = 1E-3;
V_M_TO_KV_CM = V_TO_KV / M_TO_CM;
mS_S = 1E-3;
CM_TO_M = 1/M_TO_CM;
mS_CM_TO_S_M = mS_S / CM_TO_M;

sigma_l = sigma_l_mscm * mS_CM_TO_S_M;

%% Timoshkin 
tau_f = (ep_l / ep_d + z_l / (2*z_d)) * (ep_0 * ep_d / sigma_l); 

dt = 1e-10;
t = 0:dt:1e-4;
e_0 = e_l * exp(-t ./ tau_f); % use for E in electroporation.m

w = 0.5 * sigma_l * e_l^2 * tau_f * (1-exp(-2 .* t ./ tau_f)); % Equation 9

%% Timoshkin Figure 2
figure;
subplot(1,2,1);
xlabel('t, s');
yyaxis left
loglog(t, e_0 .* (V_M_TO_KV_CM));
ylabel('E_0(t), kV/cm');

yyaxis right
loglog(t, w .* (J_TO_MJ / (M_TO_CM)^3));
ylabel('w(t), mJ/cm^3');

%% Timoshkin Full Differential Equation
% Components
taul = ep_0 * ep_l / sigma_l;
tau1 = ((2 * ep_m + ep_l) / sigma_c) * ep_0;
tau2 = ((2 * ep_l + ep_m) / sigma_l) * ep_0;
tau3 = ((ep_l - ep_m) / sigma_c) * ep_0;
tau4 = ((ep_l - ep_m) / sigma_l) * ep_0;

p = (tau1 + 0.5 * tau2 - R1R2_RATIO^3*(tau3+tau4)) / (0.5 * tau1 * tau2 - R1R2_RATIO^3 * tau3 * tau4);
q = (1 - R1R2_RATIO^3) / (0.5 * tau1 * tau2 - R1R2_RATIO^3 * tau3 * tau4);

tau5 = 2 / (p - sqrt(p^2 - 4*q));
tau6 = 2 / (p + sqrt(p^2 - 4*q));

kap = (1 / (tau5 * tau6)) - tau_f^-1 * (tau5^-1 + tau6^-1) + tau_f^-2;
lam = (1 / (tau5 * tau6)) - tau_f^-1 * (tau6^-1 - tau5^-1) - tau5^-2;
mu =  (1 / (tau5 * tau6)) - tau_f^-1 * (tau5^-1 - tau6^-1) - tau6^-2;

% Equations
A_back_1 = (exp(-t / tau_f) / kap) * (1 - ((taul + tau1)/tau_f) + ((taul * tau1) / tau_f^2));
A_back_2 = (exp(-t / tau5) / lam) *  (1 - ((taul + tau1)/tau5) +  ((taul * tau1) / tau5^2));
A_back_3 = (exp(-t / tau6) / mu) *   (1 - ((taul + tau1)/tau6) +  ((taul * tau1) / tau6^2));
A_back = A_back_1 - A_back_2 - A_back_3;
A = -3 * e_l * (tau1 * tau2 - 2 * (R1R2_RATIO)^3 * tau3 * tau4)^-1 * A_back;

B_back_1 = (exp(-t / tau_f) / kap) * (1 - ((taul + tau3)/tau_f) + ((taul * tau3) / tau_f^2));
B_back_2 = (exp(-t / tau5) / lam) *  (1 - ((taul + tau3)/tau5) +  ((taul * tau3) / tau5^2));
B_back_3 = (exp(-t / tau6) / mu) *   (1 - ((taul + tau3)/tau6) +  ((taul * tau3) / tau6^2));
B_back = B_back_1 - B_back_2 - B_back_3;
B = 3 * e_l * R1^3 * (tau1 * tau2 - 2 * (R1R2_RATIO)^3 * tau3 * tau4)^-1 * B_back;

E_p = -A + 2 * B * R1^-3; % Electric field
normalized = E_p / e_l; % Normalized

V_m = E_p * d; % Voltage

%% Timoshkin Figure 4
subplot(1,2,2);
xlabel('t, s');
plot([fliplr(-t), t], [V_m*0, 0, V_m(1:end-1)]);
hold on;
plot([fliplr(-t), t], [V_m*0, V_m*0+U_0]);
%set(gca,'Xscale','log');
ylabel('V_m, V');

if should_rerun_tfest
    uu = [V_m*0, V_m*0+U_0];
    yy = [V_m*0, 0, V_m(1:end-1)];
    sys = iddata(yy', uu', dt);
    tf = tfest(sys, 2, 2);
else
    load('timoshkin_fig4_tf.mat');
end

figure;
[mag,phase,wout] = bode(tf, {1e3, 1e9});
mag = squeeze(mag)*U_0;
hz = wout ./ 2 * pi;
plot(hz, mag);
set(gca, 'XScale', 'log')

%% Kotnik
% fs general function
fs = @(lo, li, lm, d, R) (3*lo * (3*d.*R.^2.*li+(3*d.^2*R-d.^3).*(lm-li))) ./ (2*R.^3*(lm+2*lo).*(lm+li/2)-2*(R-d).^3.*(lo-lm).*(li-lm));

% Generalizing Functions
tau = @(R, C, lo, lm, li) R*C / ((2*lo*li)/(2*lo+li) + R*lm/d);
tau_i = tau(R, C, lo, lm, li); % Instantiating tau with variables
V_mg = @(fs, E, R, theta, t, tau_i) fs.*E.*R.*cos(theta).*(1-exp(-t./tau_i)); % General function for Kotnik V_m
fs1 = fs(lo, li, lm, d, R); % Instantiating fs with variables
V_m1 = @(t, E) V_mg(fs1, E, R, 0, t, tau_i); % Use env variables, function of time and field
V_m2 = V_m1(t, e_0); % Kotnik Voltage plotted on arbitrary t, with decaying e_0 field

% Transfer Functions
if should_rerun_tfest2
    uu = [V_m*0, V_m*0+U_0];
    yy2 = [V_m2*0, 0, V_m2(1:end-1)];
    sys2 = iddata(yy2', uu', dt);
    tf2 = tfest(sys2, 2, 2);
else
    load('kotnik_tf.mat'); % referred to as tf2
end

[mag_k,phase_k,wout_k] = bode(tf2, {1e3, 1e9});
mag_k = squeeze(mag_k)*U_0;
hz_k = wout_k ./ 2 * pi;

%% Run Kotnik with fixed E
V_k_fixed = V_m1(t, fixed_e);

% Transfer Functions
if should_rerun_tfest4
    uu4 = [V_k_fixed*0, V_k_fixed*0+U_0];
    yy4 = [V_k_fixed*0, 0, V_k_fixed(1:end-1)];
    sys4 = iddata(yy4', uu4', dt);
    tf4 = tfest(sys4, 2, 2);
else
    load('kotnik_tf_fixed.mat'); % referred to as tf2
end

%% e_0 Transfer Function
e_0g = @(t) e_l * exp(-t ./ tau_f);
e_1 = e_0g(t);
% Transfer Functions
if should_rerun_tfest3
    uu3 = [e_1*0, e_1*0+U_0];
    yy3 = [e_1*0, 0, e_1(1:end-1)];
    sys3 = iddata(yy3', uu3', dt);
    tf3 = tfest(sys3, 2, 2);
else
    load('e0_tf.mat'); % referred to as tf3
end

[mag_e,phase_e,wout_e] = bode(tf3, {1e3, 1e9});
mag_e = squeeze(mag_e)*U_0;
hz_e = wout_e ./ 2 * pi;

% E_0 figure
figure; 
plot(hz_e, mag_e);
ylabel('E_0');
xlabel('f, Hz');
set(gca, 'XScale', 'log');

%% Schwann Equation 
% Fixed E for E field
fixed_e = U_0 / (z_l + 2*z_d);
schwann = @(E, R, wout, C_m, rho_i, rho_a) abs((3/2 .* E .* R)./(1+1i.*wout.*R.*C_m.*(rho_i + rho_a/2))); % Equation (1)
rho_a = 1 / lo;
rho_i = 1 / li;
schwann1 = schwann(mag_e, R, wout_e, C, rho_i, rho_a);
schwann2 = schwann(fixed_e, R, wout_e, C, rho_i, rho_a); % Schwann Fixed E

%% Full Frequency Response Plot
figure;
plot(hz_k, mag_k);
hold on;
plot(hz, mag);
set(gca, 'XScale', 'log')
% Schwann
plot(hz_e, schwann1);
plot(hz_e, schwann2);

% Bode Information, Convert 
[mag_k_f,phase_k_f,wout_k_f] = bode(tf4, {1e3, 1e9});
mag_k_f = squeeze(mag_k_f)*U_0;
hz_k_f = wout_k_f ./ 2 * pi;
plot(hz_k_f, mag_k_f);
legend('Kotnik (T E_0)', 'Timoshkin (T E_0)', 'Schwann (T E_0)', 'Schwann (Fixed E_0)', 'Kotnik (Fixed E_0)');

%% Normalized Plot
%figure;
%plot(hz_k, mag_k / max(mag_k));
%hold on;
%plot(hz, mag / max(mag));
%set(gca, 'XScale', 'log')
% Schwann
%plot(hz_e, schwann1 / max(schwann1));
%plot(hz_e, schwann2 / max(schwann2));



%% Plot Voltage vs. Time
figure;
plot([fliplr(-t), t], [V_m2*0, 0, V_m2(1:end-1)]);
xlabel('t, s')
ylabel('V_m')
title('Voltage vs Time');
hold on;
plot([fliplr(-t), t], [V_m*0, 0, V_m(1:end-1)]);
plot([fliplr(-t), t], [V_k_fixed*0, 0, V_k_fixed(1:end-1)]);
legend('Kotnik (T E_0)', 'Timoshkin (T E_0)', 'Kotnik (Fixed E_0)');