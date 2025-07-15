clear all
close all

%% Parameters : all km, kg, s
G = 6.674*10^(-20);

% Uranus param
uranus_mu_ring = 97.7 * 10^3; % km, circular orbit
uranus_mass = 86.811 * 10^24; % kg
uranus_GM = 5.7940 * 10^6; % km^3/s^2
uranus_radius = 25559; % km

% Ariel params
ariel_mass = 12.9 * 10^20; % kg
ariel_mean_radius = 578.9; % km
ariel_GM = G * ariel_mass; % km^3/s^2
ariel_sma = 190.9 * 10^3; % km
ariel_ecc = 0.0012;
ariel_v = velocity(uranus_GM,ariel_sma,ariel_sma);
uranus_2_ariel_peri_25km = ariel_sma * (1 - ariel_ecc) + ariel_mean_radius + 25; % km
uranus_2_ariel_apo_25km = ariel_sma * (1 + ariel_ecc) + ariel_mean_radius + 25; % km

% Miranda params
miranda_mass = 0.66 * 10^20; % kg
miranda_mean_radius = 235.7; % km
miranda_GM = G * miranda_mass; % km^3/s^2
miranda_sma = 129.9 * 10^3; % km
miranda_ecc = 0.0013;
miranda_v = velocity(uranus_GM,miranda_sma,miranda_sma);
uranus_2_miranda_peri_25km = miranda_sma * (1 - miranda_ecc) + miranda_mean_radius + 25; % km
uranus_2_miranda_apo_25km = miranda_sma * (1 + miranda_ecc) + miranda_mean_radius + 25; % km

s2day = 1/(60*60*24);

% Other moon param
umbriel_sma = 266 * 10^3; % km
titania_sma = 436.3 * 10^3; % km
oberon_sma = 583.5 * 10^3; % km
mab_sma = 97.74 * 10^3; % km

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

    t = linspace(0,2*pi,1000); % plotting interval
    if type == 1
        x_orb = ((a*(1-epsilon^2))./(1+epsilon.*cos(t))).*cos(t);
        y_orb = ((a*(1-epsilon^2))./(1+epsilon.*cos(t))).*sin(t);
    elseif type == 0
        x_orb = r.*cos(t);
        y_orb = r.*sin(t);
    end
end

function theta = true_anomoly(a,ecc,r)
    theta = acos( ((a*(1-ecc^2)/r)-1)/ecc );
end

%% 

crispi_a = 1.66*10^6;
crispi_sma = 8.40e+05;
crispi_p = 2*crispi_sma-1.66*10^6;
crispi_ecc = 1-(crispi_p/crispi_sma);

b = crispi_sma*sqrt(1-crispi_ecc^2);
c2uranus = crispi_sma - crispi_p;
crispi_b = [-c2uranus,-b,0];
crispi_b_mag = sqrt(b^2+c2uranus^2);
theta = pi-atan2(b,c2uranus);

v_crispi_b = velocity(uranus_GM,crispi_b_mag,crispi_sma);

new_dv_range = 0:0.5/1000:0.7;
new_flyby_perp_v_range = zeros(length(new_dv_range),1);
i = 1;

for dv_up = 0.02
    new_v_crispi = [v_crispi_b, dv_up,0];
    new_v_crispi_mag = sqrt(dv_up^2+v_crispi_b^2); %%%
    new_crispi_sma = 1 / (2/crispi_b_mag - (new_v_crispi_mag^2)/uranus_GM); %%%
    
    e_vec = cross(new_v_crispi,cross(crispi_b,new_v_crispi))/uranus_GM - ... 
            crispi_b/crispi_b_mag;
    
    new_crispi_ecc = norm(e_vec); %%%
    new_crispi_r_p = new_crispi_sma*(1-new_crispi_ecc^2); %%%
    
    new_true_a_crispi = true_anomoly(new_crispi_sma,new_crispi_ecc,ariel_sma);
    new_phi_fpa_crispi = acos((1+new_crispi_ecc*cos(new_true_a_crispi))/...
                     (sqrt(1+2*new_crispi_ecc*cos(new_true_a_crispi)+new_crispi_ecc^2)));
    
    new_v_crispi_intsec = velocity(uranus_GM,ariel_sma,new_crispi_sma);
    
    new_flyby_perp_v_crispi = new_v_crispi_intsec*sin(new_phi_fpa_crispi);
    new_flyby_perp_v_range(i) = new_flyby_perp_v_crispi;
    i=i+1;
end