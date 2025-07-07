clear all
close all

%% Parameters : 
G = 6.674*10^(-20); % km^3/kg/s^2

uranus_mu_ring = 97.7 * 10^3; % km, circular orbit
uranus_mass = 86.811 * 10^24; % kg
uranus_GM = 5.7940 * 10^6; % km^3/s^2

% Ariel parameters
ariel_mass = 12.9 * 10^20; % kg
ariel_mean_radius = 578.9; % km
ariel_GM = G * ariel_mass; % km^3/s^2
ariel_sma = 190.9 * 10^3; % km
ariel_ecc = 0.0012;
uranus_2_ariel_peri_25km = ariel_sma * (1 - ariel_ecc) + ariel_mean_radius + 25; % km
uranus_2_ariel_apo_25km = ariel_sma * (1 + ariel_ecc) + ariel_mean_radius + 25; % km

% Miranda parameters
miranda_mass = 0.66 * 10^20; % kg
miranda_mean_radius = 235.7; % km
miranda_GM = G * miranda_mass; % km^3/s^2
miranda_sma = 129.9 * 10^3; % km
miranda_ecc = 0.0013;
uranus_2_miranda_peri_25km = miranda_sma * (1 - miranda_ecc) + miranda_mean_radius + 25; % km
uranus_2_miranda_apo_25km = miranda_sma * (1 + miranda_ecc) + miranda_mean_radius + 25; % km

%% Trans-Lunar injection via Hohmann: Uranus → Ariel
% Worst case (at apoapsis)
dv_uranus_2_ariel_worst = hohmann(uranus_GM, uranus_mu_ring, uranus_2_ariel_apo_25km, ...
                                  ariel_GM, uranus_2_ariel_apo_25km, (ariel_mean_radius + 25));
% Best case (at periapsis)
dv_uranus_2_ariel_best = hohmann(uranus_GM, uranus_mu_ring, uranus_2_ariel_peri_25km, ...
                                 ariel_GM, uranus_2_ariel_peri_25km, (ariel_mean_radius + 25));

%% Trans-Lunar injection via Hohmann: Uranus → Miranda
% Worst case (at apoapsis)
dv_uranus_2_miranda_worst = hohmann(uranus_GM, uranus_mu_ring, uranus_2_miranda_apo_25km, ...
                                    miranda_GM, uranus_2_miranda_apo_25km, (miranda_mean_radius + 25));
% Best case (at periapsis)
dv_uranus_2_miranda_best = hohmann(uranus_GM, uranus_mu_ring, uranus_2_miranda_peri_25km, ...
                                   miranda_GM, uranus_2_miranda_peri_25km, (miranda_mean_radius + 25));

%% Display results
fprintf('Δv (Uranus → Ariel) worst case (apoapsis):   %.4f km/s\n', dv_uranus_2_ariel_worst);
fprintf('Δv (Uranus → Ariel) best case  (periapsis):  %.4f km/s\n', dv_uranus_2_ariel_best);
fprintf('Δv (Uranus → Miranda) worst case (apoapsis): %.4f km/s\n', dv_uranus_2_miranda_worst);
fprintf('Δv (Uranus → Miranda) best case  (periapsis):%.4f km/s\n', dv_uranus_2_miranda_best);

%% Hohmann Transfer Function
function delta_v = hohmann(GM_1, r_1_1, r_1_2, GM_2, r_2_1, r_2_2)
    % Burn 1: depart origin orbit to transfer orbit
    delta_v_1 = sqrt(GM_1/r_1_1) * (sqrt(2*r_1_2/(r_1_1 + r_1_2)) - 1);

    % Burn 2: arrive at destination circular orbit
    delta_v_2 = sqrt(GM_2/r_2_2) * (1 - sqrt(2*r_2_1/(r_2_1 + r_2_2)));

    delta_v = abs(delta_v_1) + abs(delta_v_2);
end
