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

%%
%%   Orbital diagram layout

hold on
% Plot Uranus
pos = [-uranus_radius -uranus_radius uranus_radius*2 uranus_radius*2]; 
rectangle('Position',pos,'Curvature',[1 1],'FaceColor',[0.1840 0.7450 0.9370],'EdgeColor',[0.1840 0.7450 0.9370]);
scatter(0,0,'filled','MarkerFaceColor',[0.1840 0.7450 0.9370],'SizeData',nan);

% Plot Ariel
[x_ariel_plot, y_ariel_plot] = CalculateOrbit(0, ariel_sma, 0, ariel_sma);
plot(x_ariel_plot, y_ariel_plot,'Color',[0.2310 0.6660 0.1960]);

% Plot Miranda
[x_miranda_plot, y_miranda_plot] = CalculateOrbit(0, miranda_sma, 0, miranda_sma);
plot(x_miranda_plot, y_miranda_plot,'Color',[0.0660 0.4430 0.7450]);

% Plot Umbriel
[x_umbriel_plot, y_umbriel_plot] = CalculateOrbit(0, umbriel_sma, 0, umbriel_sma);
plot(x_umbriel_plot, y_umbriel_plot,'Color',[0.0660 0.4430 0.7450]);

% Plot Titania
[x_titania_plot, y_titania_plot] = CalculateOrbit(0, titania_sma, 0, titania_sma);
plot(x_titania_plot, y_titania_plot,'Color',[0.0660 0.4430 0.7450]);

% Plot Oberon
[x_oberon_plot, y_oberon_plot] = CalculateOrbit(0, oberon_sma, 0, oberon_sma);
plot(x_oberon_plot, y_oberon_plot,'Color',[0.0660 0.4430 0.7450]);


%%   UOP
uop_equa_r_a = [-1.53125*10^6, -1.625*10^6];
uop_equa_r_p = [0.3125*10^5, 0.625*10^5];
uop_equa_sma = norm(uop_equa_r_a - uop_equa_r_p)/2;
uop_equa_ecc = 1-(norm(uop_equa_r_p)/uop_equa_sma);
uop_equa_T = period(uranus_GM,uop_equa_sma)*s2day;

% Plot 
[x_uop_plot, y_uop_plot] = CalculateOrbit(1, uop_equa_sma, uop_equa_ecc, uop_equa_sma);
plot(x_uop_plot, y_uop_plot,'Color',[0.5210 0.0860 0.8190]);

%% 
traj_num = 'traj_1';
switch (traj_num)
    case 'traj_1'
        %% Trajectory 1 - No orbit transfer, no burn , deploy direct from UOP
        
        % Deployement
        true_a_uop_flyby_1 = 2*pi + -true_anomoly(uop_equa_sma,uop_equa_ecc,ariel_sma); % from the bottom
        phi_fpa_uop_flyby_1 = acos((1+uop_equa_ecc*cos(true_a_uop_flyby_1))/...
                                   (sqrt(1+2*uop_equa_ecc*cos(true_a_uop_flyby_1)+uop_equa_ecc^2)));
        scatter([ariel_sma*cos(true_a_uop_flyby_1), ariel_sma*cos(2*pi - true_a_uop_flyby_1)], ...
                [ariel_sma*sin(true_a_uop_flyby_1), ariel_sma*sin(2*pi - true_a_uop_flyby_1)], ...
                50,'x','r')

        % UOP-Ariel flyby       
        v_uop_flyby_1 = norm([velocity(uranus_GM,ariel_sma,uop_equa_sma)*cos(phi_fpa_uop_flyby_1)-ariel_v, ...    
                             velocity(uranus_GM,ariel_sma,uop_equa_sma)*sin(phi_fpa_uop_flyby_1)]);
        
        % Time 
        E_flyby_1 = atan(sqrt((1-uop_equa_ecc)/(1+uop_equa_ecc))*tan(true_a_uop_flyby_1/2))*2 + 2*pi;
        M_flyby_1 = E_flyby_1 - (uop_equa_ecc*sin(E_flyby_1));

        % From one ariel intersection to the other
        t_max_cruise_flyby_1 = 2* (uop_equa_T - ((M_flyby_1/(2*pi)) * uop_equa_T)); 

        t_science_flyby_1 = 2*norm([ariel_r_hill, 25 + ariel_mean_radius])/ ... 
                    (velocity(uranus_GM,ariel_sma,uop_equa_sma)*sin(phi_fpa_uop_flyby_1))*s2day;

        title("Orbit 1 - No orbit transfer, deploy direct from UOP");
        legend('Uranus',"Ariel's orbit", 'Other moons','','','',"UOP's trajectory","Deployment options");
        xlim([-24*10^5 7*10^5]);
        axis equal;

    case 'traj_2'
        %% Trajectory 2 - Burn at apoapsis, crash into Uranus 

        % UOP-Ariel flyby
        true_a_uop_flyby_2 = 2*pi + -true_anomoly(uop_equa_sma,uop_equa_ecc,ariel_sma); % from the bottom
        phi_fpa_uop_flyby_2 = acos((1+uop_equa_ecc*cos(true_a_uop_flyby_1))/...
                                   (sqrt(1+2*uop_equa_ecc*cos(true_a_uop_flyby_1)+uop_equa_ecc^2)));
        
        v_uop_flyby_2 = sqrt( (velocity(uranus_GM,ariel_sma,uop_equa_sma)*cos(phi_fpa_uop_flyby_1)-ariel_v)^2 + ...
                              (velocity(uranus_GM,ariel_sma,uop_equa_sma)*sin(phi_fpa_uop_flyby_1))^2);
        
        % Time 
        
        title("Orbit 1 - No orbit transfer, deploy direct from UOP");
        legend('Uranus',"Ariel's orbit", 'Other moons','','','',"UOP's trajectory");
        axis equal;
    case 'traj_3'
        %% Trajectory 3 - Burn at periapsis
end