clear all
close all

%% Parameters : all km, kg, s
G = 6.674*10^(-20);

uranus_mu_ring = 97.7 * 10^3; % km, circular orbit
uranus_mass = 86.811 * 10^24; % kg
uranus_GM = 5.7940 * 10^6; % km^3/s^2

ariel_mass = 12.9 * 10^20; % kg
ariel_mean_radius = 578.9; % km
ariel_GM = G * ariel_mass; % km^3/s^2
ariel_sma = 190.9 * 10^3; % km
ariel_ecc = 0.0012;
ariel_v = velocity(uranus_GM,ariel_sma,ariel_sma);
uranus_2_ariel_peri_25km = ariel_sma * (1 - ariel_ecc) + ariel_mean_radius + 25; % km
uranus_2_ariel_apo_25km = ariel_sma * (1 + ariel_ecc) + ariel_mean_radius + 25; % km

miranda_mass = 0.66 * 10^20; % kg
miranda_mean_radius = 235.7; % km
miranda_GM = G * miranda_mass; % km^3/s^2
miranda_sma = 129.9 * 10^3; % km
miranda_ecc = 0.0013;
miranda_v = velocity(uranus_GM,miranda_sma,miranda_sma);
uranus_2_miranda_peri_25km = miranda_sma * (1 - miranda_ecc) + miranda_mean_radius + 25; % km
uranus_2_miranda_apo_25km = miranda_sma * (1 + miranda_ecc) + miranda_mean_radius + 25; % km

s2day = 1/(60*60*24);

%% Functions
function delta_v = hohmann(GM_1, r_1_1, r_1_2, GM_2, r_2_1, r_2_2)
    % Burn 1
    delta_v_1 = sqrt(GM_1/r_1_1) * (sqrt(2*r_1_2/(r_1_1+r_1_2))-1);

    % Burn 2
    delta_v_2 = sqrt(GM_2/r_2_2) * (1-sqrt(2*r_2_1/(r_2_1+r_2_2)));

    delta_v = abs(delta_v_1) + abs(delta_v_2);
end

function v = velocity(GM,r,a)
    v = sqrt(GM*(2/r-1/a)); 
end

function T = period(GM,a)
    T = 2*pi/sqrt(GM)*a^(3/2);
end

function [x_orb, y_orb] = CalculateOrbit(type, a, epsilon, r)
    t = 0:pi/50:2*pi; % plotting interval
    
    if type == 1
        x_orb = ((a*(1-epsilon^2))./(1+epsilon.*cos(t))).*cos(t);
        y_orb = ((a*(1-epsilon^2))./(1+epsilon.*cos(t))).*sin(t);
    elseif type == 0
        x_orb = r.*cos(t);
        y_orb = r.*sin(t);
    end
end

%% 
%%   Trans-Lunar injection via Hohmann: Uranus - Ariel
% worst case, @ apo
dv_uranus_2_ariel_worst = hohmann(uranus_GM,uranus_mu_ring,uranus_2_ariel_apo_25km,...
                                  ariel_GM, uranus_2_ariel_apo_25km,(ariel_mean_radius + 25));
% best case, @ per
dv_uranus_2_ariel_best = hohmann(uranus_GM,uranus_mu_ring,uranus_2_ariel_peri_25km,...
                                  ariel_GM, uranus_2_ariel_peri_25km,(ariel_mean_radius + 25));

%%   Trans-Lunar injection via Hohmann: Uranus - Miranda
% worst case, @ apo
dv_uranus_2_miranda_worst = hohmann(uranus_GM,uranus_mu_ring,uranus_2_miranda_apo_25km,...
                                  miranda_GM, uranus_2_miranda_apo_25km,(miranda_mean_radius + 25));
% best case, @ per
dv_uranus_2_miranda_best = hohmann(uranus_GM,uranus_mu_ring,uranus_2_miranda_peri_25km,...
                                  miranda_GM, uranus_2_miranda_peri_25km,(miranda_mean_radius + 25));

%%   UOP equatorialized traj to Ariel

uop_equa_r_a = 1.6875*10^6;
uop_equa_r_p = 1.25*10^5;
uop_equa_sma = (uop_equa_r_a + uop_equa_r_p)/2;
crispi_equa_sma = (uop_equa_r_a + ariel_mean_radius + ariel_sma + 25)/2;

dv_uop_equa_2_ariel = velocity(uranus_GM,uop_equa_r_a,crispi_equa_sma)-...
                      velocity(uranus_GM,uop_equa_r_a,uop_equa_sma);

v_uop_equa_peri = velocity(uranus_GM,ariel_mean_radius+ariel_sma+25,crispi_equa_sma);

T_equa_uop = period(uranus_GM,uop_equa_sma)*s2day;
T_equa_crispi = period(uranus_GM,crispi_equa_sma)*s2day;

%%   UOP side traj to Ariel

uop_side_r_a = sqrt((2.8750*10^6)^2+(1.75*10^6)^2);
uop_side_r_p = sqrt((0.9375*10^5)^2+(1.25*10^5)^2);
uop_side_sma = (uop_side_r_a + uop_side_r_p)/2;
crispi_side_sma = (uop_side_r_a + ariel_mean_radius + ariel_sma + 25)/2;

dv_uop_side_2_ariel = velocity(uranus_GM,uop_side_r_a,crispi_side_sma)-...
                      velocity(uranus_GM,uop_side_r_a,uop_side_sma);

v_uop_side_peri = velocity(uranus_GM,ariel_sma+ariel_mean_radius+25,crispi_side_sma);

T_side_uop = period(uranus_GM,uop_side_sma)*s2day;
T_side_crispi = period(uranus_GM,crispi_side_sma)*s2day;

%% 
