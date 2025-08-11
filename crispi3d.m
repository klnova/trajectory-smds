%% # CRISPI-MS Delta V Calculations - 3D Version
clc
clear
%%
%% ## Uranus, Major Moons & Ring Parameters [verified]
%% Uranus [source: https://nssdc.gsfc.nasa.gov/planetary/factsheet/uranusfact.html]
mu_Uranus = 5.7940e+06; % Gravitational parameter [km^3/s^2]
r_Uranus = 25559; % radius [km]

% Moons [source: https://nssdc.gsfc.nasa.gov/planetary/factsheet/uraniansatfact.html]
% Oberon
r_ob = 761.4; % radius [km]
a_ob = 583.5e+03; % semi-major axis [km]
P_ob = 13.463234; % orbital period [days]
epsilon_ob = 0.0014; % eccentricity
r_ob_orb = a_ob*(1 + epsilon_ob); % radius of orbit [km] (come back and address elliptical)
i_ob = 0.058; % inclination [radians] - approximate

% Titania
r_ti = 788.9; % radius [km]
a_ti = 436.30e+03; % semi-major axis [km]
P_ti = 8.705867; % orbital period [days]
epsilon_ti = 0.0011; % eccentricity
r_ti_orb = a_ti*(1 + epsilon_ti); % radius of orbit [km] (come back and address elliptical)
i_ti = 0.34; % inclination [radians] - approximate

% Umbriel
r_umb = 584.7; % radius [km]
a_umb = 266.00e+03; % semi-major axis [km]
P_umb = 4.144176; % orbital period [days]
epsilon_umb = 0.0039; % eccentricity
r_umb_orb = a_umb*(1 + epsilon_umb); % radius of orbit [km] (come back and address elliptical)
i_umb = 0.128; % inclination [radians] - approximate

% Ariel
r_ar = 581.1; % radius [km]
a_ar = 190.90e+03; % semi-major axis [km]
P_ar = 2.520379; % orbital period [days]
epsilon_ar = 0.0012; % eccentricity
r_ar_orb = a_ar*(1 + epsilon_ar); % radius of orbit [km] (come back and address elliptical)
p_ar = a_ar+(1*epsilon_ar^2);
V_ar_orb = sqrt((mu_Uranus)/r_ar_orb); % orbital velocity [km/s]
i_ar = 0.26; % inclination [radians] - approximate

% Miranda
r_mir = 240; % radius [km]
a_mir = 129.90e+03; % semi-major axis [km]
P_mir = 1.413479; % orbital period [days]
epsilon_mir = 0.0013; % eccentricity
r_mir_orb = a_mir*(1 + epsilon_mir); % radius of orbit [km] (come back and address elliptical)
V_mir_orb = sqrt((mu_Uranus)/r_mir_orb); % orbital velocity [km/s]
i_mir = 4.34; % inclination [radians] - approximate

% Mu Ring [source: https://nssdc.gsfc.nasa.gov/planetary/factsheet/uranringfact.html]
r_mu = 114.7e+03; % radius [km] (the outer radius of the ring is taken here, ring is 17,000 km in width)
i_mu = 0; % inclination [radians] - rings are in equatorial plane

%%
%% ## UOP Post-Equatorialized Orbit
% The following parameters were obtained from visual inspection of the final equatorialized trajectory shown in the UOP mission concept appendix c
% source: https://science.nasa.gov/wp-content/uploads/2023/10/uranus-orbiter-and-probe-appendix-c.pdf

r_a = 1.6625e+06; % apoapsis [km]
r_p = 1.25e+05; % periapsis [km]
major = r_p + r_a; % major axis [km]
a = major/2; % semi-major axis [km]
epsilon = 1 - (r_p/a); % eccentricity
F = 0-2*(a-r_p); % second ellipse focus (for plotting)
p = r_p*(1 + epsilon); % semilatus rectum
H = sqrt(p*mu_Uranus); % angular momentum of orbit [kg*m^2/s]
V_p = H/r_p;
V_a = H/r_a;
P_UOP = 2*pi*sqrt(a^3/mu_Uranus)/60/60/24; % Period of orbit [days]
i_UOP = 0; % inclination [radians] - equatorial orbit

%%
%% ## Sphere of influence
% Sphere of Influece calculation - how close do we need to be to Ariel to
% be within its SOI?

GM_Ariel = 83.43; % [km^3/s^2]
V_esc_Ariel = sqrt(2*GM_Ariel/r_ar + 35);

a_SOI = (r_ar_orb);
r_SOI_Ar = a_SOI*(((12.9*10^20)/(86.811*10^24))^(2/5));

%%
%% ## Impulsive Tangential Burn Calculations
% Initial Elliptical Orbit to Moons/Rings
% Apoapsis burn to reduce periapsis

% Ariel
[Ariel_delta_v, CRISPI_V_p, period_initial_UOP, period_final_Ariel] = CalculateDeltaVPeriapsis(r_p, (r_ar_orb - r_ar - 35), r_a);
V_flyby_Ariel = abs(V_ar_orb - CRISPI_V_p );
Ariel_delta_v_total = Ariel_delta_v + abs(7 - V_flyby_Ariel);

% Miranda
[Miranda_delta_v, CRISPI_V_p, period_initial_UOP, period_final_Miranda] = CalculateDeltaVPeriapsis(r_p, (r_mir_orb + r_mir + 35), r_a);
V_flyby_Miranda = abs(V_mir_orb - CRISPI_V_p );
Miranda_delta_v_total = Miranda_delta_v + abs(7 - V_flyby_Miranda);

% Mab
[Mab_delta_v, CRISPI_V_p, period_initial_UOP, period_final_Mab] = CalculateDeltaVPeriapsis(r_p, r_mu, r_a);

%%
%% ## 3D P
figure()
% UOP Post-Equatorialization Orbit
[x_UOP, y_UOP, z_UOP] = CalculateOrbit3D(1, a, epsilon, r_Uranus, i_UOP);
plot3(x_UOP, y_UOP, z_UOP, 'LineWidth', 2)
axis equal
grid on
hold on

% Uranus
[x_Ur, y_Ur, z_Ur] = CalculateOrbit3D(0, a, 0, r_Uranus, 0);
plot3(x_Ur, y_Ur, z_Ur, 'LineWidth', 2)

% Moon Orbits:
% Oberon
[x_ob, y_ob, z_ob] = CalculateOrbit3D(1, a_ob, epsilon_ob, r_ob, i_ob);
plot3(x_ob, y_ob, z_ob, 'LineWidth', 1.5)

% Titania
[x_ti, y_ti, z_ti] = CalculateOrbit3D(1, a_ti, epsilon_ti, r_ti, i_ti);
plot3(x_ti, y_ti, z_ti, 'LineWidth', 1.5)

% Umbriel
[x_umb, y_umb, z_umb] = CalculateOrbit3D(1, a_umb, epsilon_umb, r_umb, i_umb);
plot3(x_umb, y_umb, z_umb, 'LineWidth', 1.5)

% Ariel
[x_ar, y_ar, z_ar] = CalculateOrbit3D(1, a_ar, epsilon_ar, r_ar, i_ar);
plot3(x_ar, y_ar, z_ar, 'LineWidth', 1.5)

% Miranda
[x_mir, y_mir, z_mir] = CalculateOrbit3D(1, a_mir, epsilon_mir, r_mir, i_mir);
plot3(x_mir, y_mir, z_mir, 'LineWidth', 1.5)

% Mu Ring
[x_mu, y_mu, z_mu] = CalculateOrbit3D(0, 0, 0, r_mu, i_mu);
plot3(x_mu, y_mu, z_mu, 'LineWidth', 2)

%%
%% ## Early crossover flyby
disp('early crossover flyby')
delta_v = -0.5;
E_orbit = (((V_a + delta_v)^2)/2) - ((mu_Uranus)/r_a);
V_crossover_Ariel = sqrt(2*(E_orbit+(mu_Uranus/r_ar_orb)));
H = r_a*(V_a + delta_v);
p = H^2/mu_Uranus;
epsilon_new = 1 - p/r_a;
a = r_a/(epsilon_new+1);
nu_Ariel = acos((1/epsilon_new)*(((a*(1-epsilon_new^2))/(r_ar_orb))-1));
nu_deg = rad2deg(nu_Ariel);
r_p_at_Ariel = a/(1-epsilon_new);
phi = acos((1+epsilon_new*cos(nu_Ariel))/(sqrt(1+2*epsilon_new*cos(nu_Ariel)+epsilon_new^2)));
V_rel_perp = sin(phi)*V_crossover_Ariel;
V_rel_tan = cos(phi)*V_crossover_Ariel;
P_flyby = 2*pi*sqrt(a^3/mu_Uranus)/60/60/24; % Period of orbit [days]

% Time of flight to flyby Ariel
u_cross = acos((epsilon_new+cos((nu_Ariel)))/(1+epsilon_new*(cos(nu_Ariel))));
M_cross = u_cross - epsilon_new*sin(u_cross);
tof_CRISPI = (P_flyby/(2*pi))*(M_cross);
tof_CRISPI_from_apo = P_flyby/2 - (P_flyby/(2*pi))*(M_cross);

M_UOP = linspace(0,2*pi,1000);
[E_UOP, nu_UOP] = invKepler(M_UOP,epsilon);

% CRISPI orbit with slight inclination for visibility
i_CRISPI = 0.1; % small inclination for 3D visualization
[x_new, y_new, z_new] = CalculateOrbit3D(1, a, epsilon_new, 0, i_CRISPI);
plot3(x_new, y_new, z_new, 'LineWidth', 2)

legend('UOP','Uranus','Oberon Orbit','Titania Orbit','Umbriel Orbit','Ariel Orbit','Miranda Orbit','Mu Ring Outer Range','CRISPI Orbit','Location','NorthWest')
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
title('CRISPI Orbital Path - 3D View')
view(45, 30) % Set 3D viewing angle

%%
%% 
figure()
% plot
[x_UOP, y_UOP, z_UOP] = CalculateOrbit3D(1, a, epsilon, r_Uranus, i_UOP);
plot3(x_UOP, y_UOP, z_UOP, 'LineWidth', 2)
set(gca, 'color', [0.05, 0.05, 0.1]) % Dark space background
axis equal
grid on
hold on

% Uranus as a sphere
plot3(0, 0, 0, '.', 'MarkerSize', 32, 'Color', 'cyan')

% Ariel orbit
[x_ar, y_ar, z_ar] = CalculateOrbit3D(1, a_ar, epsilon_ar, r_ar, i_ar);
plot3(x_ar, y_ar, z_ar, 'LineWidth', 2)

% Mu Ring
[x_mu, y_mu, z_mu] = CalculateOrbit3D(0, 0, 0, r_mu, i_mu);
plot3(x_mu, y_mu, z_mu, 'LineWidth', 2)

% CRISPI orbit
[x_new, y_new, z_new] = CalculateOrbit3D(1, a, epsilon_new, 0, i_CRISPI);
plot3(x_new, y_new, z_new, 'LineWidth', 2)

legend('UOP','Uranus','Ariel','Mu Ring (outer radius)','CRISPI','Location','NorthWest')
title('CRISPI Orbital Path - 3D Visualization','FontSize', 20)
xlabel('X [km]','FontSize', 16)
ylabel('Y [km]','FontSize', 16)
zlabel('Z [km]','FontSize', 16)

%view
view(45, 30)
camlight('headlight')
lighting gouraud
colormap(cool)


%% # 3D Orbit 
function [x_orb, y_orb, z_orb] = CalculateOrbit3D(type, a, epsilon, r, inclination)
    t = linspace(0, 2*pi, 1000);
    
    if type == 1
        % Elliptical orbit
        x_orb = ((a*(1-epsilon^2))./(1+epsilon.*cos(t))).*cos(t);
        y_orb = ((a*(1-epsilon^2))./(1+epsilon.*cos(t))).*sin(t);
        z_orb = zeros(size(t)); 
        
        %inclination rotation
        if inclination ~= 0
            for i = 1:length(t)
                y_temp = y_orb(i) * cos(inclination) - z_orb(i) * sin(inclination);
                z_temp = y_orb(i) * sin(inclination) + z_orb(i) * cos(inclination);
                y_orb(i) = y_temp;
                z_orb(i) = z_temp;
            end
        end
        
    elseif type == 0
        % orb
        x_orb = r.*cos(t);
        y_orb = r.*sin(t);
        z_orb = zeros(size(t));
        
        % rotation/given
        if inclination ~= 0
            for i = 1:length(t)
                y_temp = y_orb(i) * cos(inclination) - z_orb(i) * sin(inclination);
                z_temp = y_orb(i) * sin(inclination) + z_orb(i) * cos(inclination);
                y_orb(i) = y_temp;
                z_orb(i) = z_temp;
            end
        end
    end
end

%% 
function [delta_v, V_p_new, period_initial, period_final] = CalculateDeltaVPeriapsis(r_p_initial, r_p_final, r_a)
    mu_Uranus = 5.7940e+06; % Gravitational constant [km^3/s^2]
    
    % Old orbit
    a_old = (r_p_initial + r_a)/2;
    eps_old = 1 - (r_p_initial/a_old);
    p_old = a_old*(1-eps_old^2);
    H_old = sqrt(p_old*mu_Uranus);
    V_a_old = H_old/r_a;
    V_p_old = H_old/r_p_initial;
    
    % New orbit
    a_new = (r_p_final + r_a)/2;
    eps_new = 1 - (r_p_final/a_new);
    p_new = a_new*(1-eps_new^2);
    H_new = sqrt(p_new*mu_Uranus);
    V_a_new = H_new/r_a;
    V_p_new = H_new/r_p_final;
    
    % Delta V
    delta_v = abs(V_a_old - V_a_new);
    
    % Period of Orbit
    period_initial = 2*pi*sqrt(a_old^3/mu_Uranus)/60/60/24;
    period_final = 2*pi*sqrt(a_new^3/mu_Uranus)/60/60/24;
end

function [delta_v, period_initial, period_final] = CalculateDeltaVApoapsis(r_a_initial, r_a_final, r_p)
    mu_Uranus = 5.7940e+06; % Gravitational constant [km^3/s^2]
    
    % Old orbit
    a_old = (r_a_initial + r_p)/2;
    eps_old = (r_a_initial/a_old) - 1;
    p_old = a_old*(1-eps_old^2);
    H_old = sqrt(p_old*mu_Uranus);
    V_p_old = H_old/r_p;
    V_a_old = H_old/r_a_initial;
    
    % New orbit
    a_new = (r_a_final + r_p)/2;
    eps_new = (r_a_final/a_new) - 1;
    p_new = a_new*(1-eps_new^2);
    H_new = sqrt(p_new*mu_Uranus);
    V_p_new = H_new/r_p;
    v_a_new = H_old/r_a_final;
    
    % Delta V
    delta_v = abs(V_p_old - V_p_new);
    
    % Period of Orbit
    period_initial = 2*pi*sqrt(a_old^3/mu_Uranus)/60/60/24;
    period_final = 2*pi*sqrt(a_new^3/mu_Uranus)/60/60/24;
end

function [E, nu] = invKepler(M, e)
    % Solve Kepler's equation for eccentric anomaly E
    E = M; % Initial guess
    for i = 1:10 % Newton-Raphson iterations
        E = E - (E - e*sin(E) - M)./(1 - e*cos(E));
    end
    
    % Calculate true anomaly
    nu = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
end