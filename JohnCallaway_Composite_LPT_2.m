%% John Callaway
% LPT Code
% ME EN 6520 Extra Assignment
%Text inputs to inputs.txt
%Line 1 is material properties E1, E2, G12, v12
%Line 2 is strength for failure critereon
%Line 3 is thermal coefficients
%Line 4 is ply thickness
%Line 5 is ply orientation layup
%Line 6 is loading - N_x, N_y, N_xy, M_x, M_y, M_xy
%Line 7 is temperature change

clear, clc
format short g

F = csvread ('inputs.txt');
fileID = fopen('inputs.txt');
C = textscan(fileID,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','delimiter',',');
fclose(fileID);

for k = 1:length(C)
ply(k) = C{k}(5);
end
ply(find(isnan(ply))) = [];

% Name constants from text file:
% Material properties:
E1 = F(1,1);
E2 = F(1,2);
G = F(1,3);
v = F(1,4);

% Strength Properties:
S_L_plus = F(2,1);
S_L_minus = F(2,2);
S_T_plus = F(2,3);
S_T_minus = F(2,4);
S_LT = F(2,5);

% Coefficients of thermal expansion
alpha_1 = F(3,1);
alpha_2 = F(3,2);

% Ply thickness\
t = F(4,1);

% Laminate stacking
L = ply;

% Forces and moments
N_x = F(6,1);
N_y = F(6,2);
N_xy = F(6,3);
M_x = F(6,4);
M_y = F(6,5);
M_xy = F(6,6);

% Temperature change
Temp = F(7,1);

%1. Calculate Q and S matrices:

R = [1 0 0
    0 1 0
    0 0 2];

S = [1/E1 -v/E1 0
    -v/E1 1/E2 0
    0 0 1/G];
    
Q = inv(S);

% 1. Calculate Q_bar and S_bar matrices
for i = 1:length(L)
    theta = L(i);
    
    T = [(cosd(theta)).^2 (sind(theta)).^2 2*sind(theta).*cosd(theta)
(sind(theta)).^2 (cosd(theta)).^2 -2*sind(theta).*cosd(theta)
-sind(theta).*cosd(theta) sind(theta).*cosd(theta) (cosd(theta)).^2-(sind(theta)).^2];

    Q_bar(:,(i*3-2):i*3) = T\Q*R*T/(R);
    S_bar(:,(i*3-2):i*3) = inv(Q_bar(:,(i*3-2):i*3));
    
end
    
% Calculate Z values for plies:

for i= 1:length(L)+1
        z(i) = -length(L)/2*t+((i-1)*t);
end

% 2. Calculate ABD matrix and inverse

A = zeros(3);
B = zeros(3);
D = zeros(3);
ABD = zeros(6);

for i= 1:3
     for j= 1:length(L)
            A(i,1) = A(i,1)+Q_bar(i,(j*3)-2)*(z(j+1)-z(j));
            A(i,2) = A(i,2)+Q_bar(i,(j*3)-1)*(z(j+1)-z(j));
            A(i,3) = A(i,3)+Q_bar(i,(j*3)-0)*(z(j+1)-z(j));
            B(i,1) = B(i,1)+0.5*(Q_bar(i,(j*3)-2)*(z(j+1)^2-z(j)^2));
            B(i,2) = B(i,2)+0.5*(Q_bar(i,(j*3)-1)*(z(j+1)^2-z(j)^2));
            B(i,3) = B(i,3)+0.5*(Q_bar(i,(j*3)-0)*(z(j+1)^2-z(j)^2));
            D(i,1) = D(i,1)+(1/3)*(Q_bar(i,(j*3)-2)*(z(j+1)^3-z(j)^3));
            D(i,2) = D(i,2)+(1/3)*(Q_bar(i,(j*3)-1)*(z(j+1)^3-z(j)^3));
            D(i,3) = D(i,3)+(1/3)*(Q_bar(i,(j*3)-0)*(z(j+1)^3-z(j)^3));
     end
end

% make values 0 that should be 0:
for i= 1:3
    for j= 1:3
        if A(i,j) < 1*10^-12 && A(i,j) > -1*10^-12
            A(i,j) = 0;
        end
        if B(i,j) < 1*10^-12 && B(i,j) > -1*10^-12
            B(i,j) = 0;
        end
        if D(i,j) < 1*10^-12 && D(i,j) > -1*10^-12
            D(i,j) = 0;
        end
    end
end

ABD(1:3,1:3) = A;
ABD(4:6,1:3) = B;
ABD(1:3,4:6) = B;
ABD(4:6,4:6) = D;
ABDinv = ABD^-1;

% 3. Apparent laminate stiffness properties:
E_x = 1/(length(L)*t*ABDinv(1,1));
E_y = 1/(length(L)*t*ABDinv(2,2));
v_xy = -ABDinv(1,2)/ABDinv(1,1);
G_xy = 1/(length(L)*t*ABDinv(3,3));
E_fx = 12/((length(L)*t)^3*ABDinv(4,4));
E_fy = 12/((length(L)*t)^3*ABDinv(5,5));

% 4. Midplane strains and curvatures
% Calculate thermal expansion coefficients in x-y:
% Calculate N and M thermal
alpha = [alpha_1,alpha_2,0]';
N_T = [0,0,0]';
M_T = [0,0,0]';
for i = 1:length(L)
    theta = L(i);
    
    T = [(cosd(theta)).^2 (sind(theta)).^2 2*sind(theta).*cosd(theta)
(sind(theta)).^2 (cosd(theta)).^2 -2*sind(theta).*cosd(theta)
-sind(theta).*cosd(theta) sind(theta).*cosd(theta) (cosd(theta)).^2-(sind(theta)).^2];
    alpha_xy = inv(T)*alpha;
    Q_bar_T = T\Q*R*T/(R);
    N_T = N_T+Temp*Q_bar_T*alpha_xy*(z(i+1)-z(i));
    M_T = M_T+0.5*Temp*Q_bar_T*alpha_xy*(z(i+1)^2-z(i)^2);
end
N = [N_x,N_y,N_xy]'+N_T;
M = [M_x,M_y,M_xy]'+M_T;
NM = [N
    M];

strain_mid = ABDinv*NM;

eps_x_mid = strain_mid(1);
eps_y_mid = strain_mid(2);
eps_xy_mid = strain_mid(3);
K_x = strain_mid(4);
K_y = strain_mid(5);
K_xy = strain_mid(6);

% 5. Calculation of strains at all points in laminate and plot

% Calculation of mechanical strain
eps_mech_x = eps_x_mid + z*K_x;
eps_mech_y = eps_y_mid + z*K_y;
eps_mech_xy = eps_xy_mid + z*K_xy;

eps_mech = [eps_mech_x;eps_mech_y;eps_mech_xy];

figure(1)
plot(eps_mech_x,z,'r-',eps_mech_y,z,'b-',eps_mech_xy,z,'g-')
set(gca,'ydir','reverse')
legend('\epsilon_x','\epsilon_y','\epsilon_x_y')
grid on
xlabel('Strain, []')
ylabel('z,[m]')
title('Part 5 - Strain')

% 6. Calculation of stresses at top and bottom of ply:

% Get repeating z values at ply interfaces
for i=1:length(L)+1
    zplot((i*2)-1)=z(i);
    zplot(i*2)=z(i);
end
    zplot(2*(length(L)+1)) = [];
    zplot(1) = [];
    

for i = 1:length(L)
    theta = L(i);
    T = [(cosd(theta)).^2 (sind(theta)).^2 2*sind(theta).*cosd(theta)
(sind(theta)).^2 (cosd(theta)).^2 -2*sind(theta).*cosd(theta)
-sind(theta).*cosd(theta) sind(theta).*cosd(theta) (cosd(theta)).^2-(sind(theta)).^2];
    alpha_xy = inv(T)*alpha;
    eps_therm = alpha_xy*Temp;
    eps1 = eps_mech(:,i)-eps_therm;
    eps2 = eps_mech(:,i+1)-eps_therm;
    stress_top(:,i) = Q_bar(:,3*i-2:3*i)*(eps1);
    stress_bottom(:,i) = Q_bar(:,3*i-2:3*i)*(eps2);
end

for i = 1:length(L)
    stress_xy(:,i*2-1) = stress_top(:,i);
    stress_xy(:,i*2) = stress_bottom(:,i);
end

% 7. Plots of global stresses for each ply:

figure(2)
plot(stress_xy(1,:),zplot,'-*')
legend('\sigma_x')
set(gca,'ydir','reverse')
grid on
xlabel('\sigma_x,[Pa]')
ylabel('z,[m]')
title('\sigma_x')
    
figure(3)
plot(stress_xy(2,:),zplot,'-*')
legend('\sigma_y')
set(gca,'ydir','reverse')
grid on
xlabel('\sigma_y, [Pa]')
ylabel('z,[m]')
title('\sigma_y')
    
figure(4)
plot(stress_xy(3,:),zplot,'-*')
legend('\sigma_x_y')
set(gca,'ydir','reverse')
grid on
xlabel('\sigma_x_y,[Pa]')
ylabel('z, [m]')
title('\sigma_x_y')

% 8. Material coordinates on top and bottom of each ply
    
for i = 1:length(L)
   theta = L(i);
   
   T = [(cosd(theta)).^2 (sind(theta)).^2 2*sind(theta).*cosd(theta)
    (sind(theta)).^2 (cosd(theta)).^2 -2*sind(theta).*cosd(theta)
    -sind(theta).*cosd(theta) sind(theta).*cosd(theta) (cosd(theta)).^2-(sind(theta)).^2];
   
    alpha_xy = inv(T)*alpha;
    eps_therm = alpha_xy*Temp;
    eps1 = eps_mech(:,i)-eps_therm;
    eps2 = eps_mech(:,i+1)-eps_therm;
    stress_top(:,i) = Q_bar(:,3*i-2:3*i)*(eps1);
    stress_bottom(:,i) = Q_bar(:,3*i-2:3*i)*(eps2);
    stress_top_12(:,i) = T*stress_top(:,i);
    stress_bottom_12(:,i) = T*stress_bottom(:,i);
end

for i = 1:length(L)
    stress_12(:,i*2-1) = stress_top_12(:,i);
    stress_12(:,i*2) = stress_bottom_12(:,i);
end

% 9. Hashin Failure Criteria:

Fiber_Failure = ((stress_12(1,:).^2)./(S_L_plus*S_L_minus))+...
        ((1/S_L_plus-1/S_L_minus).*stress_12(1,:))+...
        ((stress_12(3,:).^2)./(S_LT^2));
Resin_Failure = ((stress_12(2,:).^2)./(S_T_plus*S_T_minus))+...
        ((1/S_T_plus-1/S_T_minus).*stress_12(2,:))+...
        ((stress_12(3,:).^2)./(S_LT^2));
    
mode = 'No Failure';
M=0;

i = 1;
while i <= length(L)
    theta2((i*2)-1)=L(i);
    theta2(i*2)=L(i);
    i = i+1;
end

if max(Fiber_Failure) >= 1
    mode = 'Fiber Failure';
    [M, I] = max(Fiber_Failure);
    z_fail_f = zplot(I);
    ply_fail_f = theta2(I);
end
if max(Resin_Failure) >= 1
    [M2, I2] = max(Resin_Failure);
    if M2 > M
        mode = 'Resin Failure';
        z_fail_r = zplot(I2);
        ply_fail_r = theta2(I2);
    end
end

% Outputs to .txt file:
i = 1;
fid = fopen( 'outputs_JohnCallaway.txt', 'wt' );
fprintf(fid,'John Callaway \n');
fprintf(fid,'LPT Code - Spring 2017 \n');
fprintf(fid,'\n Part 1 \n');
for i = 1:length(L)
    fprintf(fid,'\n Ply %i: \n',i);
    fprintf(fid,'\n The Q matrix is:\n');
    fprintf(fid,'%8.3g %8.3g %8.3g\n',Q);
    fprintf(fid,'\n The Q bar Matrix is:\n');
    fprintf(fid,'%8.3g %8.3g %8.3g\n',Q_bar(:,i*3-2:i*3));
    fprintf(fid,'\n The S matrix is:\n');
    fprintf(fid,'%8.3g %8.3g %8.3g\n',S);
    fprintf(fid,'\n The S bar Matrix is:\n');
    fprintf(fid,'%8.3g %8.3g %8.3g\n',S_bar(:,i*3-2:i*3));
end
fprintf(fid,'\n Part 2: \n');
fprintf(fid,'\n The ABD matrix is: \n');
fprintf(fid,'%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g \n',ABD);
fprintf(fid,'\n The ABD inverse matrix is: \n');
fprintf(fid,'%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g \n',ABDinv);

fprintf(fid,'\n Part 3: \n');
fprintf(fid,'\n E_x is: \n');
fprintf(fid,'%8.2g \n', E_x);
fprintf(fid,'\n E_y is: \n');
fprintf(fid,'%8.2g \n', E_y);
fprintf(fid,'\n v_xy is: \n');
fprintf(fid,'%8.2g \n', v_xy);
fprintf(fid,'\n G_xy is: \n');
fprintf(fid,'%8.2g \n', G_xy);
fprintf(fid,'\n E_fx is: \n');
fprintf(fid,'%8.2g \n', E_fx);
fprintf(fid,'\n E_fy is: \n');
fprintf(fid,'%8.2g \n', E_fy);

fprintf(fid,'\n Part 4: \n');
fprintf(fid,'\n Midplane Strain in x is: \n');
fprintf(fid,'%8.4g \n',eps_x_mid);
fprintf(fid,'\n Midplane Strain in y is: \n');
fprintf(fid,'%8.4g \n',eps_y_mid);
fprintf(fid,'\n Midplane Strain in xy is: \n');
fprintf(fid,'%8.4g \n',eps_xy_mid);
fprintf(fid,'\n Curvature in x is: \n');
fprintf(fid,'%8.4g \n',K_x);
fprintf(fid,'\n Curvature in y is: \n');
fprintf(fid,'%8.4g \n',K_y);
fprintf(fid,'\n Curvature in xy is: \n');
fprintf(fid,'%8.4g \n',K_xy);

fprintf(fid,'\n Part 6: \n');
fprintf(fid,'\n Global Coordinate Stresses: \n');
for i = 1:length(L)
    fprintf(fid,'\n Ply %i Top:',i); 
    fprintf(fid,' %8.4g %8.4g %8.4g \n',stress_top(:,i));
    fprintf(fid,'\n Ply %i Bottom:',i); 
    fprintf(fid,' %8.4g %8.4g %8.4g \n',stress_bottom(:,i));
end

fprintf(fid,'\n Part 8: \n');
fprintf(fid,'\n Material Coordinate Stresses: \n');
for i = 1:length(L)
    fprintf(fid,'\n Ply %i Top:',i); 
    fprintf(fid,' %8.4g %8.4g %8.4g \n',stress_top_12(:,i));
    fprintf(fid,'\n Ply %i Bottom:',i); 
    fprintf(fid,' %8.4g %8.4g %8.4g \n',stress_bottom_12(:,i));
end

fprintf(fid,'\n Part 9: \n');

if strcmp(mode,'No Failure')
    fprintf(fid,'\n No failure occured \n');
end
if strcmp(mode,'Fiber Failure')
    fprintf(fid,'Ply %i failed first, its orientation is %d degrees \n',ceil(I/2),ply_fail_f);
    fprintf(fid,'The location of the failure is %8.3g m \n',z_fail_f');
    fprintf(fid,'The failure occured in the fibers \n');
end
if strcmp(mode,'Resin Failure')
    fprintf(fid, 'Ply %i failed first, its orientation is %d degrees \n',ceil(I2/2),ply_fail_r);
    fprintf(fid, 'The location of the failure is %8.3g m \n',z_fail_r');
    fprintf(fid, 'The failure occured in the resin \n');
end