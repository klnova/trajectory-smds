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
ariel_r_hill = 3220; % km
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
function v = velocity(GM,r,a)
    v = sqrt(GM*(2/r-1/a)); 
end

function T = period(GM,a)
    T = 2*pi * sqrt((a^3)/GM);
end

function theta = true_anomoly(a,ecc,r)
    theta = acos( ((a*(1-ecc^2)/r)-1)/ecc );
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




%%   UOP
uop_equa_r_a = [-1.53125*10^6, -1.625*10^6];
uop_equa_r_p = [0.3125*10^5, 0.625*10^5];
uop_equa_sma = norm(uop_equa_r_a - uop_equa_r_p)/2;
uop_equa_ecc = 1-(norm(uop_equa_r_p)/uop_equa_sma);
uop_equa_T = period(uranus_GM,uop_equa_sma)*s2day;
%% Flyby

% velocity at 25km from ariel's surface
v_range = 1:0.01:10; alt_ariel = 25;

flyby_sma_range = zeros(length(v_range),1);
v_inf_range = zeros(length(v_range),1);
d_range = zeros(length(v_range),1);
phi_range = zeros(length(v_range),1);
r_p = alt_ariel + ariel_mean_radius;
i = 1;

% plot flyby parameters based on flyby velocity 
for v = v_range
    flyby_sma = 1 / (2/r_p - (v^2) /ariel_GM);
    v_inf = sqrt(-ariel_GM/flyby_sma);
    phi = 2*asin(1 / (1 + (r_p*v_inf^2)/ariel_GM) );
    d = r_p * sqrt(1 + (2*ariel_GM) / (r_p*v_inf^2) );

    flyby_sma_range(i) = flyby_sma;
    v_inf_range(i) = v_inf;
    phi_range(i) = phi;
    d_range(i) = d;
    
    i = i + 1;
end

plot(v_range,v_inf_range)
title("v inf in [km/s] over flyby velocity [km/s]")
grid on

figure
plot(v_range,d_range)
title("D [km] over flyby velocity [km/s]")
grid on

figure
plot(v_range,phi_range)
title("phi [rad] over flyby velocity [km/s]")
grid on




