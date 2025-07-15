clear all
close all

%% Parameters : all km, kg, s
G = 6.674*10^(-20);

uranus_mu_ring = 97.7 * 10^3; % km, circular orbit
uranus_mass = 86.811 * 10^24; % kg
uranus_GM = 5.7940 * 10^6; % km^3/s^2
uranus_radius = 25559; % km

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

function theta = true_anomoly(a,ecc,r)
    theta = acos( ((a*(1-ecc^2)/r)-1)/ecc );
end



%% 
%%   UOP equatorialized traj to Miranda

uop_equa_mir_r_a = 1.6875*10^6;
uop_equa_mir_r_p = 1.25*10^5;
uop_equa_mir_sma = (uop_equa_mir_r_a + uop_equa_mir_r_p)/2;
crispi_equa_mir_sma = (uop_equa_mir_r_a + miranda_sma + miranda_mean_radius + 25)/2;

dv_uop_equa_2_miranda = velocity(uranus_GM,uop_equa_mir_r_a,crispi_equa_mir_sma)-...
                        velocity(uranus_GM,uop_equa_mir_r_a,uop_equa_mir_sma);

v_uop_mir_equa_peri = velocity(uranus_GM,miranda_sma+miranda_mean_radius+25,crispi_equa_mir_sma);

T_mir_equa_uop = period(uranus_GM,uop_equa_mir_sma)*s2day;
T_mir_equa_crispi = period(uranus_GM,crispi_equa_mir_sma)*s2day;

%%   Ideal

ideal_v = 8.9;
ideal_sma = 1/ ((2/(ariel_sma+ariel_mean_radius+25))-(ideal_v^2)/uranus_GM);
velocity(uranus_GM,ariel_mean_radius+ariel_sma+25,abs(ideal_sma));

%%   UOP Ariel nodes

uop_equa_r_a = 1.6875*10^6;
uop_equa_r_p = 1.25*10^5;
uop_equa_sma = (uop_equa_r_a + uop_equa_r_p)/2;
uop_equa_ecc = 1-(uop_equa_r_p/uop_equa_sma);
crispi_equa_sma = (uop_equa_r_a + ariel_sma + ariel_mean_radius + 25)/2;

dv_uop_equa_2_ariel = abs(velocity(uranus_GM,uop_equa_r_a,crispi_equa_sma)-...
                          velocity(uranus_GM,uop_equa_r_a,uop_equa_sma));

v_uop_equa_intsec = velocity(uranus_GM,ariel_sma,uop_equa_sma);

true_a_noburn = true_anomoly(uop_equa_sma,uop_equa_ecc,ariel_sma);
phi_fpa_noburn = acos((1+uop_equa_ecc*cos(true_a_noburn))/(sqrt(1+2*uop_equa_ecc*cos(true_a_noburn)+uop_equa_ecc^2)));

flyby_perp_v_noburn = v_uop_equa_intsec*sin(phi_fpa_noburn);

T_equa_uop = period(uranus_GM,uop_equa_sma)*s2day;
T_equa_crispi = period(uranus_GM,crispi_equa_sma)*s2day;

%%   Increase eccentricity - Ariel nodes

uop_equa_r_a = 1.6875*10^6;
uop_equa_r_p = 1.25*10^5;
uop_equa_sma = (uop_equa_r_a + uop_equa_r_p)/2;
uop_equa_ecc = 1-(uop_equa_r_p/uop_equa_sma);
crispi_equa_sma = (uop_equa_r_a + uranus_radius)/2;
crispi_equa_ecc = 1-(uranus_radius/crispi_equa_sma);

dv_crispi_equa_2_ariel = abs(velocity(uranus_GM,uop_equa_r_a,crispi_equa_sma)-...
                          velocity(uranus_GM,uop_equa_r_a,uop_equa_sma));

v_crispi_equa_intsec = velocity(uranus_GM,ariel_sma,crispi_equa_sma);
true_a_crispi = true_anomoly(crispi_equa_sma,crispi_equa_ecc,ariel_sma);
phi_fpa_crispi = acos((1+crispi_equa_ecc*cos(true_a_crispi))/(sqrt(1+2*crispi_equa_ecc*cos(true_a_crispi)+crispi_equa_ecc^2)));

flyby_perp_v_crispi = v_crispi_equa_intsec*sin(phi_fpa_crispi);

T_equa_uop = period(uranus_GM,uop_equa_sma)*s2day;
T_equa_crispi = period(uranus_GM,crispi_equa_sma)*s2day;

%%% delta-v vs. speed plot

dv_range = -1*(0:0.5/1000:0.5);
flyby_perp_v_range = zeros(length(dv_range),1);
i = 1;

for dv = dv_range
    crispi_sma = 1/(2/uop_equa_r_a-(((dv + velocity(uranus_GM,uop_equa_r_a,uop_equa_sma))^2)/uranus_GM));
    crispi_ecc = (uop_equa_r_a/crispi_sma) - 1;
    v_crispi_intsec = velocity(uranus_GM,ariel_sma,crispi_sma);

    true_a = true_anomoly(crispi_sma,crispi_ecc,ariel_sma);
    phi_fpa = acos((1+crispi_ecc*cos(true_a))/(sqrt(1+2*crispi_ecc*cos(true_a)+crispi_ecc^2)));
    flyby_perp_v = v_crispi_intsec*sin(phi_fpa);
    
    flyby_perp_v_range(i) = flyby_perp_v;
    i=i+1;
end 

plot(-1*dv_range,flyby_perp_v_range);

title("Relative flyby velocity with respect to \Delta v")
xlabel("$\Delta v \: [km/s]$",'Interpreter','latex')
ylabel("$Flyby \: velocity \: [km/s]$",'Interpreter','latex')

grid on;
grid minor;
ax=gca;
ax.TickLabelInterpreter = 'latex';