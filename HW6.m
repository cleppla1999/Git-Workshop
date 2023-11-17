% HW6 Script 
% 11/2/2023
% Carson Leppla
% Housekeeping 
clc; clear; close all;

%% Homework Givens
mu_sun = 1.32712428 * 10^11;
mu_earth = 3.986004418 * 10^5;
Au = 149597870.7;

%% Problem 1 
Date1 = 2460012.5;
R1 = [-1.451526*10^8, 2.872782*10^7, 1.245389*10^7];
V1 = [-6.766088,-26.810318,-11.620955];

Date2 = 2460146.5;
R2 = [2.454556 * 10^7, -9.611903 * 10^7, - 4.480278 * 10^7];
V2 = [33.882657, 7.872920, 1.398796];

TOF = (Date2 - Date1)*(60^2*24);
[a,e,p,dtheta] = LambertSolver(mu_sun,R1,R2,"less",TOF);

% solve for the dv's using f and g functions 
% use f and functions to determine dv1 and dv2
theta1 = acosd(1/e*(p/norm(R1)- 1));
theta2 = acosd(1/e*(p/norm(R2)-1));

% Find theta2 at the current state 
[f,g,f_dot,g_dot] = fg(norm(R2),norm(R1),p,mu_sun,dtheta);

V1f = (R2-f*R1)/g;
V2i = f_dot*R1+g_dot*V1f;

dV1 = norm(V1f - V1);
dV2 = norm(V2 - V2i);
%% Test example 
% R1 = [-654 13605 1997];
% V1 = [-5.53 0.849 0.6830];
% 
% R2 = [7284 -19341 -3264];
% V2 = [3.07 2.63 0.444];
% TOF = 5 * 60^2;
% 
%[a,e,p,delta_theta,dV1,dV2] = LambertAndDVSolver(mu_earth,[R1 V1],[R2 V2],"greater",TOF);
%% Problem 2
% Part A done use LAmberAndDVSolver
% Part B
Earth_Data = importdata("HW6_Ephem_Earth.txt");
for i = 4:length(Earth_Data.textdata)
    Earth_Epoch(i-3) = str2double(Earth_Data.textdata{i});
end

Mars_Data = importdata("HW6_Ephem_Mars.txt");
for i = 4:length(Mars_Data.textdata)
    Mars_Epoch(i-3) = str2double(Mars_Data.textdata{i});
end

% Grab the TOF, R1, V1, R2, V2 Vectors 
X1 = Earth_Data.data(:,:);
X2 = Mars_Data.data(:,:);

for i = 1:length(X2)
    for j = 1:length(X1)
    TOF = (Mars_Epoch(i)- Earth_Epoch(j))*60^2*24;
    [~,~,~,~,dV1(i,j),dV2(i,j)] = LambertAndDVSolver(mu_sun,X1(j,:),X2(i,:),"less",TOF);
    end
end

%% Plot the Vinifity Mag for Earth Departure
x = Earth_Epoch;
y = Mars_Epoch;
[X,Y] = meshgrid(x,y);

figure;
contour(X,Y,dV1,50);
colorbar;
title("Vinifinity for Earth Departure");
ylabel("Mars Epoch");
xlabel("Earth Epoch");

figure;
contour(X,Y,dV2,50);
colorbar;
title("Vinifinity for Mars Arival");
ylabel("Mars Epoch");
xlabel("Earth Epoch");

%% Functions
function [a,e,p,delta_theta,dV1,dV2] = LambertAndDVSolver(mu,X1,X2,transfer_angle,TOF)
%LAMBERTSOLVER function to solve lamberts problem 
% inputs: r1, r2, transfer_angle
% outputs: at,et
r1 = X1(1:3);
V1 = X1(4:6);

r2 = X2(1:3);
V2 = X2(4:6);
% Everything grabbed from the state vectors 
r1_norm = norm(r1);
r2_norm = norm(r2);
% Step 1: Calculate the transfer angle
delta_theta = abs(acosd(dot(r1,r2)/(r1_norm*r2_norm)));
if transfer_angle == "greater"
    if delta_theta > 180
        delta_theta = delta_theta;
    else
        delta_theta = 360 - delta_theta;
    end
elseif transfer_angle == "less"
    if delta_theta < 180
        delta_theta = delta_theta;
    else
        delta_theta = 360 - delta_theta;
    end
end

%Step 2: Calcualte geometric quantities 
c = sqrt(r1_norm^2+r2_norm^2-2*r1_norm*r2_norm*cosd(delta_theta));
s = 0.5*(r1_norm+r2_norm+c);

%Step 3: DEtermine if transfer uses elliptical or hyperblic orbit
% Check TOF with parablolic orbit
% Find TOFp 
if delta_theta < 180
    TOFp = 1/3*sqrt(2/mu)*(s^(3/2)-(s-c)^(3/2));
else 
    TOFp = 1/3*sqrt(2/mu)*(s^(3/2)+(s-c)^(3/2));
end

if TOF < TOFp
    warning("TOF does not condone using an elliptical orbit")
end

% Step 4: Determine if transfer requres shorter or longer TOF via
% comparsion to TOF along the minimum energy transfer
a_m = s/2;
n_m = sqrt(mu/(a_m^3));
alpha_m = pi;
Betha_m_0 = 2*asin(sqrt((s-c)/s));

if delta_theta < 180
    Betha_m = Betha_m_0;
else
    Betha_m = -Betha_m_0;
end

% Calcualte TOF_min
TOF_min = 1/n_m *((alpha_m-Betha_m)-(sin(alpha_m)-sin(Betha_m)));

% use deltha_theta and TOF to find the right qudrant
% Step 5: Iteratively solve the TOF equation to calculate a 
% intilize TOF diff
% a starting as a min
% da will be set to 50 km
diff = 10000;
a = a_m;
da = 1;
% Count variable
k = 1;

while diff > 1*10^-5 
    % starting a = a_m
    if k > 1
        a = a + da;
    end
    % Find alpha0
     alpha0 = 2*asin(sqrt(s/(2*a)));

    % Find betha0
     betha0 = 2*asin(sqrt((s-c)/(2*a)));

     %Check the quadrants 
     % Change beta based on delta theat
     if delta_theta < 180
         betha = betha0;
     else
         betha = -betha0;
     end

     %Check the quadrants based on TOF 
     if TOF < TOF_min
         alpha = alpha0;
     else
         alpha = 2*pi - alpha0;
     end

     % recompute n
     n = sqrt(mu/(a^3));
     % Compute new TOF
     TOFi = 1/n*(alpha-betha-(sin(alpha)-sin(betha)));

     % Compute the TOF difference and update the count
     diff = TOF - TOFi;
     k = k + 1;
end

% Now compute the eccentricity
p = (4*a*(s-r1_norm)*(s-r2_norm))/(c^2)*sin((alpha+betha)/2)^2;
e = sqrt(1-p/a);

% solve for the dv's using f and g functions 
% use f and functions to determine dv1 and dv2
% theta1 = acosd(1/e*(p/norm(R1)- 1));
% theta2 = acosd(1/e*(p/norm(R2)-1));

% Find theta2 at the current state 
[f,g,f_dot,g_dot] = fg(r2_norm,r1_norm,p,mu,delta_theta);

V1f = (r2-f*r1)/g;
V2i = f_dot*r1+g_dot*V1f;

dV1 = norm(V1f - V1);
dV2 = norm(V2 - V2i);
end

function [f,g,f_dot,g_dot] = fg(r,r0,p,mu,dtheta)
%FG function to determine f and g and f_dot and g_dot 
% inputs r,r0,mu,dtheta
f = 1 - r/p *(1-cosd(dtheta));
g = (r*r0)/sqrt(mu*p)*sind(dtheta);
f_dot = sqrt(mu/p)*tand(dtheta/2)*((1-cosd(dtheta))/p - 1/r - 1/r0);
g_dot = 1 - (r0/p)*(1-cosd(dtheta));
end

function [a,e,p,delta_theta] = LambertSolver(mu,r1,r2,transfer_angle,TOF)
%LAMBERTSOLVER function to solve lamberts problem 
% inputs: r1, r2, transfer_angle
% outputs: at,et
r1_norm = norm(r1);
r2_norm = norm(r2);
% Step 1: Calculate the transfer angle
delta_theta = abs(acosd(dot(r1,r2)/(r1_norm*r2_norm)));
if transfer_angle == "greater"
    if delta_theta > 180
        delta_theta = delta_theta;
    else
        delta_theta = 360 - delta_theta;
    end
elseif transfer_angle == "less"
    if delta_theta < 180
        delta_theta = delta_theta;
    else
        delta_theta = 360 - delta_theta;
    end
end

%Step 2: Calcualte geometric quantities 
c = sqrt(r1_norm^2+r2_norm^2-2*r1_norm*r2_norm*cosd(delta_theta));
s = 0.5*(r1_norm+r2_norm+c);

%Step 3: DEtermine if transfer uses elliptical or hyperblic orbit
% Check TOF with parablolic orbit
% Find TOFp 
if delta_theta < 180
    TOFp = 1/3*sqrt(2/mu)*(s^(3/2)-(s-c)^(3/2));
else 
    TOFp = 1/3*sqrt(2/mu)*(s^(3/2)+(s-c)^(3/2));
end

if TOF < TOFp
    warning("TOF does not condone using an elliptical orbit")
end

% Step 4: Determine if transfer requres shorter or longer TOF via
% comparsion to TOF along the minimum energy transfer
a_m = s/2;
n_m = sqrt(mu/(a_m^3));
alpha_m = pi;
Betha_m_0 = 2*asin(sqrt((s-c)/s));

if delta_theta < 180
    Betha_m = Betha_m_0;
else
    Betha_m = -Betha_m_0;
end

% Calcualte TOF_min
TOF_min = 1/n_m *((alpha_m-Betha_m)-(sin(alpha_m)-sin(Betha_m)));

% use deltha_theta and TOF to find the right qudrant
% Step 5: Iteratively solve the TOF equation to calculate a 
% intilize TOF diff
% a starting as a min
% da will be set to 50 km
diff = 10000;
a = a_m;
da = 1;
% Count variable
k = 1;

while diff > 1*10^-5
    % starting a = a_m
    if k > 1
        a = a + da;
    end
    % Find alpha0
     alpha0 = 2*asin(sqrt(s/(2*a)));

    % Find betha0
     betha0 = 2*asin(sqrt((s-c)/(2*a)));

     %Check the quadrants 
     % Change beta based on delta theat
     if delta_theta < 180
         betha = betha0;
     else
         betha = -betha0;
     end

     %Check the quadrants based on TOF 
     if TOF < TOF_min
         alpha = alpha0;
     else
         alpha = 2*pi - alpha0;
     end

     % recompute n
     n = sqrt(mu/(a^3));
     % Compute new TOF
     TOFi = 1/n*(alpha-betha-(sin(alpha)-sin(betha)));

     % Compute the TOF difference and update the count
     diff = TOF - TOFi;
     k = k + 1;
end

% Now compute the eccentricity
p = (4*a*(s-r1_norm)*(s-r2_norm))/(c^2)*sin((alpha+betha)/2)^2;
e = sqrt(1-p/a);

end


