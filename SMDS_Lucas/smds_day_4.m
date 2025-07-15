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
%%   Time for arc in orbit to semi-minor axis
uop_equa_r_a = 1.6875*10^6;
uop_equa_r_p = 1.25*10^5;
uop_equa_sma = (uop_equa_r_a + uop_equa_r_p)/2;
uop_equa_ecc = 1-(uop_equa_r_p/uop_equa_sma);

b = uop_equa_sma*sqrt(1-uop_equa_ecc^2);
c2uranus = uop_equa_sma - uop_equa_r_p;
uop_equa_b = [-c2uranus,-b,0];
uop_equa_b_mag = sqrt(b^2+c2uranus^2);
theta = 2*pi-(pi-atan2(b,c2uranus));

% Time in orbit calc 
T_uop = 2*pi*sqrt((uop_equa_sma^3)/uranus_GM);
E = atan(sqrt((1-uop_equa_ecc)/(1+uop_equa_ecc))*tan(theta/2))*2 + 2*pi;
M = E - (uop_equa_ecc*sin(E));

t = ((M/(2*pi)) * T_uop) - T_uop/2;

%%   Orbital diagram layout
uop_equa_r_a = 1.6875*10^6;
uop_equa_r_p = 1.25*10^5;
uop_equa_sma = (uop_equa_r_a + uop_equa_r_p)/2;
uop_equa_ecc = 1-(uop_equa_r_p/uop_equa_sma);

% Plot Uranus
pos = [-uranus_radius -uranus_radius uranus_radius*2 uranus_radius*2]; 
rectangle('Position',pos,'Curvature',[1 1],'FaceColor',[0.1840 0.7450 0.9370],'EdgeColor',[0.1840 0.7450 0.9370]);
hold on;
[x_uranus_plot, y_uranus_plot] = CalculateOrbit(0, uranus_radius, 0, uranus_radius);
plot(x_uranus_plot, y_uranus_plot,'Color',[0.1840 0.7450 0.9370]);
hold on;
% Plot Ariel
[x_ariel_plot, y_ariel_plot] = CalculateOrbit(0, ariel_sma, 0, ariel_sma);
plot(x_ariel_plot, y_ariel_plot,'Color',[0.2310 0.6660 0.1960]);
hold on;
% Plot Miranda
[x_miranda_plot, y_miranda_plot] = CalculateOrbit(0, miranda_sma, 0, miranda_sma);
plot(x_miranda_plot, y_miranda_plot,'Color',[0.0660 0.4430 0.7450]);
hold on;
% Plot Umbriel
[x_umbriel_plot, y_umbriel_plot] = CalculateOrbit(0, umbriel_sma, 0, umbriel_sma);
plot(x_umbriel_plot, y_umbriel_plot,'Color',[0.0660 0.4430 0.7450]);
hold on;
% Plot Titania
[x_titania_plot, y_titania_plot] = CalculateOrbit(0, titania_sma, 0, titania_sma);
plot(x_titania_plot, y_titania_plot,'Color',[0.0660 0.4430 0.7450]);
hold on;
% Plot Oberon
[x_oberon_plot, y_oberon_plot] = CalculateOrbit(0, oberon_sma, 0, oberon_sma);
plot(x_oberon_plot, y_oberon_plot,'Color',[0.0660 0.4430 0.7450]);
hold on;

% Plot UOP
[x_uop_plot, y_uop_plot] = CalculateOrbit(1, uop_equa_sma, uop_equa_ecc, uop_equa_sma);
plot(x_uop_plot, y_uop_plot,'Color',[0.5210 0.0860 0.8190]);
hold on;

%%   CRISPI orbit 1
% Define orbit
crispi_equa_sma = (uop_equa_r_a + uranus_radius)/2;
crispi_equa_ecc = 1-(uranus_radius/crispi_equa_sma);
r_hill = 3220;

dv_crispi_equa_2_ariel = abs(velocity(uranus_GM,uop_equa_r_a,crispi_equa_sma)-...
                          velocity(uranus_GM,uop_equa_r_a,uop_equa_sma));

v_crispi_equa_intsec = velocity(uranus_GM,ariel_sma,crispi_equa_sma);
true_a_crispi = true_anomoly(crispi_equa_sma,crispi_equa_ecc,ariel_sma);
phi_fpa_crispi = acos((1+crispi_equa_ecc*cos(true_a_crispi))/(sqrt(1+2*crispi_equa_ecc*cos(true_a_crispi)+crispi_equa_ecc^2)));
flyby_perp_v_crispi = v_crispi_equa_intsec*sin(phi_fpa_crispi);

% Times
T_equa_uop = period(uranus_GM,uop_equa_sma)*s2day;
T_equa_crispi = period(uranus_GM,crispi_equa_sma)*s2day;
 
E = atan(sqrt((1-crispi_equa_ecc)/(1+crispi_equa_ecc))*tan((-true_a_crispi+2*pi)/2))*2 + 2*pi;
M = E - (crispi_equa_ecc*sin(E));
t_cruise = ((M/(2*pi)) * T_equa_crispi) - T_equa_crispi/2;

t_science = 2*sqrt(((r_hill)^2-(25+ariel_mean_radius)^2))/flyby_perp_v_crispi;


% Plot CRISPI
[x_crispi_1_plot, y_crispi_1_plot] = CalculateOrbit(1, crispi_equa_sma, crispi_equa_ecc, crispi_equa_sma);
plot(x_crispi_1_plot, y_crispi_1_plot,'Color',[0.8190 0.0150 0.5450]);
hold on;

% Plot burn
quiver(-uop_equa_r_a,0,0,2*10^5,'Color','red','LineWidth',1.8);

title("Orbit 1");
legend('Uranus',"Ariel's orbit", '','','','','UOP trajectory','CRISPI trajectory','Burn: 0.3681 km/s delta-V');
axis equal;
xlim([-18*10^5 7*10^5])

%%   CRISPI orbit 2



