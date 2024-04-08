%In this code I will produce the second stage of the multipole graph, i.e.
%   the plots of the dipole + octupole magnetic field components

%Firstly, I will define some constants in cgs units


clear all
clf

R_sol = 6.957 * 10^10; 
M_sol = 1.989 * 10^33;
G = 6.674 * 10^-8;
year = 3.1536 * 10^7;
day = 86400;

beta = 0.6;  %Beta Parameter

%Now, I will define values for a typical test stars

M_star = M_sol;
R_star = 4 * R_sol;       
Period = 7 * day;
M_acc_rate = ((10^-8) * M_sol) / year;

B_dipole = 1000;   %Dipole magnetic field strength
B_octupole = 4000;  %Octupole magnetic field strength

u = (B_dipole) * ((R_star^3)/2);   %NOT SURE ABOUT THIS

%Now, I will define the truncation radius for this test star, which will
%   also be the largest possible value of rm

R_t = (beta * (u^(4/7)) * ((2*G*M_star)^(-1/7)) * (M_acc_rate^(-2/7))) / R_sol;

   




B_dipole = 1000;
B_octupole = 5000;

r_null = (3/4 * B_octupole/B_dipole)^0.5  %In solar radii


rm = linspace (1.1, 2.5, 3);
rm2 = linspace (1.1, r_null ,3);

theta1 = 0:pi/10:pi;
phi1 = 0:2*pi/10:2*pi;
gamma1 = (1 - phi1) * 360;
theta_dip = 0;
theta_oct = [0, 180];    %For Parallel magnetic field
dS = 10e-2;

val = 100;






for p = 1:numel(phi1)

    for t = 1:(numel(theta1) - 1)

        for i = 2 %1:numel(rm)
    

            [x_coords, y_coords, z_coords, BB1] = multipole(rm(i), theta1(t), phi1(p), gamma1(p), B_dipole, B_octupole, theta_dip, theta_oct(1), dS, R_star);
            [x_coords2, y_coords2, z_coords2, BB2] = multipole2(rm2(i), phi1(p), gamma1(p), B_dipole, B_octupole, theta_dip, theta_oct(1), dS, R_star);
            [x_coords3, y_coords3, z_coords3, BB3] = multipole3(rm2(i), phi1(p), gamma1(p), B_dipole, B_octupole, theta_dip, theta_oct(1), dS, R_star);
        
        
        
            X3D = [x_coords, -x_coords, x_coords, -x_coords];
            Y3D = [-y_coords, y_coords, y_coords, -y_coords];
            Z3D = [z_coords, z_coords, -z_coords, -z_coords];

        
            X3D_2 = [x_coords2, -x_coords2, x_coords2, -x_coords2];
            Y3D_2 = [-y_coords2, y_coords2, y_coords2, -y_coords2];
            Z3D_2 = [z_coords2, z_coords2, -z_coords2, -z_coords2];

        
            X3D_3 = [x_coords3, -x_coords3, x_coords3, -x_coords3];
            Y3D_3 = [-y_coords3, y_coords3, y_coords3, -y_coords3];
            Z3D_3 = [z_coords3, z_coords3, -z_coords3, -z_coords3];
 
        
        
        
        
            figure(1)
%             subplot(2,1,1)

        
            plot3(X3D, Y3D, Z3D, 'k-', LineWidth=1);
            plot3(X3D_2, Y3D_2, Z3D_2, 'b-', LineWidth=1);
            plot3(X3D_3, Y3D_3, Z3D_3, 'r-', LineWidth=1);
        
            xlabel('X axis', 'FontSize', 15)
            ylabel('Y axis', 'FontSize', 15)
            zlabel('Z axis', 'FontSize', 15)
            title('Parallel Multipole Closed Magnetic Field Lines 3D', 'FontSize', 18)
            ax.FontSize = 15;
            grid on
            axis equal
            hold on


        end
    end

end




theta2 = 0:pi/100:pi;
phi2 = pi:2*pi/100:3*pi;
gamma2 = (1 - phi2) * 360;
BM = zeros(numel(theta2), numel(phi2));
R = 1;

for i = 1:numel(phi2)
    
    [B_r, B_theta, B_phi] = field1(R, theta2, phi2, gamma2, B_dipole, B_octupole, theta_dip, theta_oct(1), R_star);

    
    BM(:,i) = B_r;

    BBBB(:,i) = (B_r.^2 + B_theta.^2 + B_phi.^2).^0.5;

end

[x,y,z] = sphere(val);
hs1 = surf(x,y,z, flipud(BM)/1000, 'EdgeColor','none');
colormap("turbo");
colorbar;
q1 = get(hs1);
a = colorbar;
a.FontSize = 15; 
a.TickLabelInterpreter = 'Latex';
a.Label.String = 'kG';











%Now for the Dipole


for p = 1:numel(phi1)

    for t = 1:(numel(theta1) - 1)

        for i = 2 %1:numel(rm)
    
            [x_coords, y_coords, z_coords, BB] = dipole(rm(i), theta1(t), phi1(p), gamma1(p), B_dipole, theta_dip, dS, R_star);
        
        
            X3D = [x_coords, -x_coords, x_coords, -x_coords];
            Y3D = [-y_coords, y_coords, y_coords, -y_coords];
            Z3D = [z_coords, z_coords, -z_coords, -z_coords];

       
       
%             subplot(2,1,2)
            figure(5)

        
            plot3(X3D, Y3D, Z3D, 'k-', LineWidth=1);
        
            xlabel('X axis', 'FontSize', 15)
            ylabel('Y axis', 'FontSize', 15)
            zlabel('Z axis', 'FontSize', 15)
            title('Dipole Closed Magnetic Field Lines 3D', 'FontSize', 18)
            ax.FontSize = 15;
            grid on
            axis equal
            hold on


        end
    end

end






BD = zeros(numel(theta2), numel(phi2));
R = 1;

for i = 1:numel(phi2)
    
    [B_r, B_theta, B_phi] = field1(R, theta2, phi2, gamma2, B_dipole, 0, theta_dip, theta_oct(1), R_star);
    
    BD(:,i) = B_r;

end

[x,y,z] = sphere(val);
hs1 = surf(x,y,z, flipud(BD)/1000, 'EdgeColor','none');
colormap("turbo");
colorbar;
q1 = get(hs1);
a = colorbar;
a.FontSize = 15; 
a.TickLabelInterpreter = 'Latex';
a.Label.String = 'kG';





figure(2)

subplot(1,2,1)
set(gca, 'TicklabelInterpreter', 'latex')
set(gca, 'FontSize', 28)
tt = [-90 90];
pp = [0 360];
image(pp, -tt, flipud(BM/1000), 'CDataMapping', 'scaled')
title('Parallel Multipole Magnetic Map', 'FontSize', 18)
xlabel("Longitude (degrees)",'Interpreter','latex')
ylabel("Latitude (degrees)",'Interpreter','latex')
a = colorbar;
a.Label.String = 'kG';
colormap('turbo')



subplot(1,2,2)
set(gca, 'TicklabelInterpreter', 'latex')
set(gca, 'FontSize', 28)
tt = [-90 90];
pp = [0 360];
image(pp, -tt, flipud(BD/1000), 'CDataMapping', 'scaled')
title('Dipole Magnetic Map', 'FontSize', 18)
xlabel("Longitude (degrees)",'Interpreter','latex')
ylabel("Latitude (degrees)",'Interpreter','latex')
a = colorbar;
a.Label.String = 'kG';
colormap('turbo')











































theta3 = (0:120/500:120)*pi/180;
phi3 = 0:2*pi/500:2*pi;
gamma3 = (1 - phi3) * 360;

B1 = zeros(numel(theta3), numel(phi3));
B2 = zeros(numel(theta3), numel(phi3));
B3 = zeros(numel(theta3), numel(phi3));


for i = 1:numel(phi3)


    [B_r, B_theta, B_phi] = field1(R, theta3, phi3, gamma3, B_dipole, B_octupole, theta_dip, theta_oct(1), R_star);

    B1(:,i) = B_r;
    B2(:,i) = B_theta;
    B3(:,i) = B_phi;

end






pos = [0, 30, 60, 90, 120];
RR = linspace(0, 120, 501);
Az = linspace(0, 360, 501);

P = B1/1000;
P1 = rescale(B2/1000, min(P, [], 'all'), max(P, [], 'all'));

figure(3)
set(gca, 'FontSize', 16)

subplot(1,3,1)
[~,c] = polarPcolor(RR, Az, P, 'Nspokes', 5, 'circlesPos', pos, 'ncolor', 10);
caxis([-2.2,2.2])
title('B_R');

subplot(1,3,2)
[~,b] = polarPcolor(RR, Az, B2/1000, 'Nspokes', 5, 'circlesPos', pos, 'ncolor', 10);
caxis([-2.2,2.2])
title('B_{theta}');

subplot(1,3,3)
title('B_{phi}');
[~,a] = polarPcolor(RR, Az, B3/1000, 'Nspokes', 5, 'circlesPos', pos, 'ncolor', 10);
colormap('turbo')
% caxis([-2.2,2.2])










































































%Function to plot the circle

function [x_units, y_units] = circle(x, y, r, N)
    
    angle = 0:(pi/N):(2*pi);
    x_units = x + (r * cos(angle));
    y_units = y + (r * sin(angle));

end









%This is the function to plot the multipole coordinates, I will try to
%integrate the closed lines given starting from the theta=0 line.

function [x_coords, y_coords, z_coords, BB] = multipole(r, theta, phi, gamma, B_dipole, B_octupole, theta_dip, theta_oct, dS, R_star)



    x_coords = [r.*sin(theta).*cos(phi)];
    y_coords = [r.*sin(theta).*sin(phi)];
    z_coords = [r.*cos(theta)];


    beta = theta_oct - theta_dip;


    mu_d = (R_star^3 * B_dipole)/2;
    mu_o = (R_star^5 * B_octupole)/4;

    while r > 1
    
    
    
        B_r = (B_dipole * ((1/r)^3) * cos(theta) * cos(theta_dip)) + 0.5 * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(theta_oct)^2) - 3) * cos(theta) * cos(theta_oct));
        B_theta = (0.5 * B_dipole * ((1/r)^3) * sin(theta) * cos(theta_dip)) + (3/8) * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(theta_oct)^2) - 1) * sin(theta) * cos(theta_oct)); 
        
        u_r_d = (mu_d*sin(theta)*cos(phi)*cos(gamma)*sin(beta)) + (mu_d*sin(theta)*sin(phi)*sin(gamma)*sin(beta)) + (mu_d*cos(theta)*cos(beta));
        u_theta_d = (mu_d*cos(theta)*cos(phi)*cos(gamma)*sin(beta)) + (mu_d*cos(theta)*sin(phi)*sin(gamma)*sin(beta)) - (mu_d*sin(theta)*cos(beta));
        u_phi_d = (-mu_d*sin(phi)*cos(gamma)*sin(beta)) + (mu_d*cos(phi)*sin(gamma)*sin(beta));
    
        u_r_o = (mu_o*sin(theta)*cos(phi)*cos(gamma)*sin(beta)) + (mu_o*sin(theta)*sin(phi)*sin(gamma)*sin(beta)) + (mu_o*cos(theta)*cos(beta));
        u_theta_o = (mu_o*cos(theta)*cos(phi)*cos(gamma)*sin(beta)) + (mu_o*cos(theta)*sin(phi)*sin(gamma)*sin(beta)) - (mu_o*sin(theta)*cos(beta));
        u_phi_o = (-mu_o*sin(phi)*cos(gamma)*sin(beta)) + (mu_o*cos(phi)*sin(gamma)*sin(beta));
    
    
        B_phi = (((B_dipole/2) * (1/r)^3) * (-u_phi_d/mu_d)) + (((B_octupole/4) * (1/r)^5) * (-(15/2)*((u_r_o/mu_o)^2)*(u_phi_o/mu_o) + (1/16)*((u_phi_o/mu_o))));

        B = (B_r.^2 + B_theta.^2 + B_phi.^2).^0.5;
        dr = (B_r ./ B) .* dS;
        dtheta = (B_theta ./ (r.*B)) .* dS;
        dphi = (B_phi ./ (r.*B.*sin(theta))) .* dS;
    
        
        r = r - dr;
        theta = theta - dtheta;
        phi = phi + dtheta;



        x_coords = [x_coords; r.*sin(theta).*cos(phi)];
        y_coords = [y_coords; r.*sin(theta).*sin(phi)];       
        z_coords = [z_coords; r.*cos(theta)];

        


    
    end

    BB = B;

end







%Below, the two functions to plot the multipole component inside the null
%radius 




function [x_coords2, y_coords2, z_coords2, BB2] = multipole2(r, phi, gamma, B_dipole, B_octupole, theta_dip, theta_oct, dS, R_star)

    theta = acos((3/5 + (2/5) * (B_dipole/B_octupole) * (1/r)^2)^0.5);

    x_coords2 = [r.*sin(theta).*cos(phi)];
    y_coords2 = [r.*sin(theta).*sin(phi)];
    z_coords2 = [r.*cos(theta)];
    

    beta = theta_oct - theta_dip;

    mu_d = (R_star^3 * B_dipole)/2;
    mu_o = (R_star^5 * B_octupole)/4;

    while r > 1
    
    
    
        B_r = (B_dipole * ((1/r)^3) * cos(theta) * cos(theta_dip)) + 0.5 * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(0)^2) - 3) * cos(theta) * cos(theta_oct));
        B_theta = (0.5 * B_dipole * ((1/r)^3) * sin(theta) * cos(theta_dip)) + (3/8) * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(0)^2) - 1) * sin(theta) * cos(theta_oct)); 


        u_r_d = (mu_d*sin(theta)*cos(phi)*cos(gamma)*sin(beta)) + (mu_d*sin(theta)*sin(phi)*sin(gamma)*sin(beta)) + (mu_d*cos(theta)*cos(beta));
        u_theta_d = (mu_d*cos(theta)*cos(phi)*cos(gamma)*sin(beta)) + (mu_d*cos(theta)*sin(phi)*sin(gamma)*sin(beta)) - (mu_d*sin(theta)*cos(beta));
        u_phi_d = (-mu_d*sin(phi)*cos(gamma)*sin(beta)) + (mu_d*cos(phi)*sin(gamma)*sin(beta));
    
        u_r_o = (mu_o*sin(theta)*cos(phi)*cos(gamma)*sin(beta)) + (mu_o*sin(theta)*sin(phi)*sin(gamma)*sin(beta)) + (mu_o*cos(theta)*cos(beta));
        u_theta_o = (mu_o*cos(theta)*cos(phi)*cos(gamma)*sin(beta)) + (mu_o*cos(theta)*sin(phi)*sin(gamma)*sin(beta)) - (mu_o*sin(theta)*cos(beta));
        u_phi_o = (-mu_o*sin(phi)*cos(gamma)*sin(beta)) + (mu_o*cos(phi)*sin(gamma)*sin(beta));
    
    
        B_phi = (((B_dipole/2) * (1/r)^3) * (-u_phi_d/mu_d)) + (((B_octupole/4) * (1/r)^5) * (-(15/2)*((u_r_o/mu_o)^2)*(u_phi_o/mu_o) + (1/16)*((u_phi_o/mu_o))));

    
        B = (B_r.^2 + B_theta.^2 + B_phi.^2).^0.5;
        dr = (B_r ./ B) .* dS;
        dtheta = (B_theta ./ (r.*B)) .* dS;
        dphi = (B_phi ./ (r.*B.*sin(theta))) .* dS;
    
        
        r = r + dr;
        theta = theta + dtheta;
        phi = phi + dtheta;



        x_coords2 = [x_coords2; r.*sin(theta).*cos(phi)];
        y_coords2 = [y_coords2; r.*sin(theta).*sin(phi)];       
        z_coords2 = [z_coords2; r.*cos(theta)];


    
    end

    BB2 = B;


end







function [x_coords3, y_coords3, z_coords3, BB3] = multipole3(r, phi, gamma, B_dipole, B_octupole, theta_dip, theta_oct, dS, R_star)

    theta = acos((3/5 + (2/5) * (B_dipole/B_octupole) * (1/r)^2)^0.5);

    x_coords3 = [r.*sin(theta).*cos(phi)];
    y_coords3 = [r.*sin(theta).*sin(phi)];
    z_coords3 = [r.*cos(theta)];

    beta = theta_oct - theta_dip;

    mu_d = (R_star^3 * B_dipole)/2;
    mu_o = (R_star^5 * B_octupole)/4;

    while r > 1
    
    
    
        B_r = (B_dipole * ((1/r)^3) * cos(theta) * cos(theta_dip)) + 0.5 * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(0)^2) - 3) * cos(theta) * cos(theta_oct));
        B_theta = (0.5 * B_dipole * ((1/r)^3) * sin(theta) * cos(theta_dip)) + (3/8) * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(0)^2) - 1) * sin(theta) * cos(theta_oct));


        u_r_d = (mu_d*sin(theta)*cos(phi)*cos(gamma)*sin(beta)) + (mu_d*sin(theta)*sin(phi)*sin(gamma)*sin(beta)) + (mu_d*cos(theta)*cos(beta));
        u_theta_d = (mu_d*cos(theta)*cos(phi)*cos(gamma)*sin(beta)) + (mu_d*cos(theta)*sin(phi)*sin(gamma)*sin(beta)) - (mu_d*sin(theta)*cos(beta));
        u_phi_d = (-mu_d*sin(phi)*cos(gamma)*sin(beta)) + (mu_d*cos(phi)*sin(gamma)*sin(beta));
    
        u_r_o = (mu_o*sin(theta)*cos(phi)*cos(gamma)*sin(beta)) + (mu_o*sin(theta)*sin(phi)*sin(gamma)*sin(beta)) + (mu_o*cos(theta)*cos(beta));
        u_theta_o = (mu_o*cos(theta)*cos(phi)*cos(gamma)*sin(beta)) + (mu_o*cos(theta)*sin(phi)*sin(gamma)*sin(beta)) - (mu_o*sin(theta)*cos(beta));
        u_phi_o = (-mu_o*sin(phi)*cos(gamma)*sin(beta)) + (mu_o*cos(phi)*sin(gamma)*sin(beta));
    
    
        B_phi = (((B_dipole/2) * (1/r)^3) * (-u_phi_d/mu_d)) + (((B_octupole/4) * (1/r)^5) * (-(15/2)*((u_r_o/mu_o)^2)*(u_phi_o/mu_o) + (1/16)*((u_phi_o/mu_o))));
    
        B = (B_r.^2 + B_theta.^2 + B_phi.^2).^0.5;
        dr = (B_r ./ B) .* dS;
        dtheta = (B_theta ./ (r.*B)) .* dS;
        dphi = (B_phi ./ (r.*B.*sin(theta))) .* dS;
    
        
        r = r - dr;
        theta = theta - dtheta;
        phi = phi - dtheta;



        x_coords3 = [x_coords3; r.*sin(theta).*cos(phi)];
        y_coords3 = [y_coords3; r.*sin(theta).*sin(phi)];       
        z_coords3 = [z_coords3; r.*cos(theta)];


    
    end

    BB3 = B;


end





































function [B_r, B_theta, B_phi] = field1(r, theta, phi, gamma, B_dipole, B_octupole, theta_dip, theta_oct, R_star)




    beta = theta_oct - theta_dip;


    mu_d = (R_star^3 * B_dipole)/2;
    mu_o = (R_star^5 * B_octupole)/4;

    
    
    u_r_d = (mu_d.*sin(theta).*cos(phi).*cos(gamma).*sin(beta)) + (mu_d.*sin(theta).*sin(phi).*sin(gamma).*sin(beta)) + (mu_d.*cos(theta).*cos(beta));
    u_theta_d = (mu_d.*cos(theta).*cos(phi).*cos(gamma).*sin(beta)) + (mu_d.*cos(theta).*sin(phi).*sin(gamma).*sin(beta)) - (mu_d.*sin(theta).*cos(beta));
    u_phi_d = (-mu_d.*sin(phi).*cos(gamma).*sin(beta)) + (mu_d.*cos(phi).*sin(gamma).*sin(beta));

    u_r_o = (mu_o.*sin(theta).*cos(phi).*cos(gamma).*sin(beta)) + (mu_o.*sin(theta).*sin(phi).*sin(gamma).*sin(beta)) + (mu_o.*cos(theta).*cos(beta));
    u_theta_o = (mu_o.*cos(theta).*cos(phi).*cos(gamma).*sin(beta)) + (mu_o.*cos(theta).*sin(phi).*sin(gamma).*sin(beta)) - (mu_o.*sin(theta).*cos(beta));
    u_phi_o = (-mu_o.*sin(phi).*cos(gamma).*sin(beta)) + (mu_o.*cos(phi).*sin(gamma).*sin(beta));




    B_r = (B_dipole .* ((1/r).^3) .* cos(theta) .* cos(theta_dip)) + 0.5 .* (B_octupole .* ((1/r).^5) .* (5 .* (cos(theta).^2) * (cos(0).^2) - 3) .* cos(theta) .* cos(theta_oct));
    B_theta = (0.5 .* B_dipole .* ((1/r).^3) .* sin(theta) .* cos(theta_dip)) + (3/8) .* (B_octupole .* ((1/r).^5) .* (5 .* (cos(theta).^2) .* (cos(theta_oct).^2) - 1) .* sin(theta) .* cos(theta_oct)); 
    B_phi = (((B_dipole./2) * (1./r).^3) .* (-u_phi_d./mu_d)) + (((B_octupole./4) * (1./r).^5) .* (-(15/2).*((u_r_o./mu_o).^2).*(u_phi_o./mu_o) + (1/16).*((u_phi_o./mu_o))));

   
    





end



























function [x_coords, y_coords, z_coords, BB] = dipole(r, theta, phi, gamma, B_dipole, theta_dip, dS, R_star)



    x_coords = [r.*sin(theta).*cos(phi)];
    y_coords = [r.*sin(theta).*sin(phi)];
    z_coords = [r.*cos(theta)];
    

    mu_d = (R_star^3 * B_dipole)/2;
    beta = theta_dip;


    while r > 1
    
    
    
        B_r = (B_dipole * ((1/r)^3) * cos(theta) * cos(theta_dip));
        B_theta = (0.5 * B_dipole * ((1/r)^3) * sin(theta) * cos(theta_dip)); 
        
        u_r_d = (mu_d*sin(theta)*cos(phi)*cos(gamma)*sin(beta)) + (mu_d*sin(theta)*sin(phi)*sin(gamma)*sin(beta)) + (mu_d*cos(theta)*cos(beta));
        u_theta_d = (mu_d*cos(theta)*cos(phi)*cos(gamma)*sin(beta)) + (mu_d*cos(theta)*sin(phi)*sin(gamma)*sin(beta)) - (mu_d*sin(theta)*cos(beta));
        u_phi_d = (-mu_d*sin(phi)*cos(gamma)*sin(beta)) + (mu_d*cos(phi)*sin(gamma)*sin(beta));


        B_phi = (((B_dipole/2) * (1/r)^3) * (-u_phi_d/mu_d));



        B = (B_r.^2 + B_theta.^2 + B_phi.^2).^0.5;
        dr = (B_r ./ B) .* dS;
        dtheta = (B_theta ./ (r.*B)) .* dS;
        dphi = (B_phi ./ (r.*B.*sin(theta))) .* dS;
    
        
        r = r - dr;
        theta = theta - dtheta;
        phi = phi + dtheta;



        x_coords = [x_coords; r.*sin(theta).*cos(phi)];
        y_coords = [y_coords; r.*sin(theta).*sin(phi)];       
        z_coords = [z_coords; r.*cos(theta)];

        


    
    end

    BB = B;

end




















function [varargout] = polarPcolor(R,theta,Z,varargin)
% [h,c] = polarPcolor1(R,theta,Z,varargin) is a pseudocolor plot of matrix
% Z for a vector radius R and a vector angle theta.
% The elements of Z specify the color in each cell of the
% plot. The goal is to apply pcolor function with a polar grid, which
% provides a better visualization than a cartesian grid.
%
%% Syntax
%
% [h,c] = polarPcolor(R,theta,Z)
% [h,c] = polarPcolor(R,theta,Z,'Ncircles',10)
% [h,c] = polarPcolor(R,theta,Z,'Nspokes',5)
% [h,c] = polarPcolor(R,theta,Z,'Nspokes',5,'colBar',0)
% [h,c] = polarPcolor(R,theta,Z,'Nspokes',5,'labelR','r (km)')
%
% INPUT
%	* R :
%        - type: float
%        - size: [1 x Nrr ] where Nrr = numel(R).
%        - dimension: radial distance.
%	* theta :
%        - type: float
%        - size: [1 x Ntheta ] where Ntheta = numel(theta).
%        - dimension: azimuth or elevation angle (deg).
%        - N.B.: The zero is defined with respect to the North.
%	* Z :
%        - type: float
%        - size: [Ntheta x Nrr]
%        - dimension: user's defined .
%	* varargin:
%        - Ncircles: number  of circles for the grid definition.
%        - autoOrigin: 'on' (the first circle of the plar grid has a radius
%          equal to the lowest value of R) or 'off'.
%        - Nspokes: number of spokes for the grid definition.
%        - colBar: display the colorbar or not.
%        - labelR: title for radial axis.
%        - RtickLabel: Tick label for the radial axis.
%        - colormap: Colormap for the pcolor function
%        - ncolor: Number of colors in the colorbar and pcolor
%        - circlesPos: position of the circles with respect to the origin
%        (it overwrites Ncircles if necessary)
%
%
% OUTPUT
% h: returns a handle to a SURFACE object.
% c: returns a handle to a COLORBAR object.
%
%% Examples
% R = linspace(3,10,100);
% theta = linspace(0,180,360);
% Z = linspace(0,10,360)'*linspace(0,10,100);
% figure
% polarPcolor(R,theta,Z,'Ncircles',3)
%
%% Author
% Etienne Cheynet, University of Stavanger, Norway. 23/10/2019
% see also pcolor
%
%%  InputParseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('Ncircles',5);
p.addOptional('autoOrigin','on');
p.addOptional('Nspokes',8);
p.addOptional('labelR','');
p.addOptional('RtickLabel',[]);
p.addOptional('colBar',1);
p.addOptional('Rscale','linear');
p.addOptional('colormap','parula');
p.addOptional('ncolor',[]);
p.addOptional('typeRose','meteo'); % 'meteo' or 'default'
p.addOptional('circlesPos',[]);
p.parse(varargin{:});
Ncircles = p.Results.Ncircles;
Nspokes = p.Results.Nspokes ;
labelR = p.Results.labelR ;
RtickLabel = p.Results.RtickLabel ;
colBar = p.Results.colBar ;
Rscale = p.Results.Rscale ;
autoOrigin = p.Results.autoOrigin ;
myColorMap = p.Results.colormap ;
ncolor = p.Results.ncolor ;
circPos = p.Results.circlesPos ;
typeRose = p.Results.typeRose ;
if ~isempty(circPos)
    Origin = max([min(circPos),min(R)]);
    circPos(circPos<min(R))=[];
    circPos(circPos>max(R))=[];
elseif strcmpi(autoOrigin,'on')
    Origin = min(R);
elseif strcmpi(autoOrigin,'off')
    Origin = 0;
else
    error(' ''autoOrigin'' must be ''on'' or ''of'' ')
end
if Origin==0 && strcmpi(Rscale,'log')
    warning(' The origin cannot be set to 0 if R is expressed on a logarithmic axis. The value ''Rmin'' is used instead')
    Origin = min(R);
end
if isempty(circPos)
    if ~isempty(RtickLabel)
        if numel(RtickLabel)~=Ncircles
            error(' The radial ticklabel must be equal to Ncircles');
        end
        if any(cellfun(@ischar,RtickLabel)==0)
            error(' The radial ticklabel must be a cell array of characters');
        end
    end
end
if ~isempty(circPos)
    circPos = unique([min(R),circPos,max(R)]);
end
%% Preliminary checks
% case where dimension is reversed
Nrr = numel(R);
Noo = numel(theta);
if isequal(size(Z),[Noo,Nrr]) && Noo~=Nrr,
    Z=Z';
end
% case where dimension of Z is not compatible with theta and R
if ~isequal(size(Z),[Nrr,Noo])
    fprintf('\n')
    fprintf([ 'Size of Z is : [',num2str(size(Z)),'] \n']);
    fprintf([ 'Size of R is : [',num2str(size(R)),'] \n']);
    fprintf([ 'Size of theta is : [',num2str(size(theta)),'] \n\n']);
    error(' dimension of Z does not agree with dimension of R and Theta')
end
%% data plot
rMin = min(R);
rMax = max(R);
thetaMin=min(theta);
thetaMax =max(theta);
if strcmpi(typeRose,'meteo')
    theta = theta;
elseif strcmpi(typeRose,'default')
    theta = 90-theta;
else
    error('"type" must be "meteo" or "default" ');
end
% Definition of the mesh
cax = newplot;
Rrange = rMax - rMin; % get the range for the radius
[rNorm] = getRnorm(Rscale,Origin,R,Rrange); % getRnorm is a nested function
YY = (rNorm)'*cosd(theta);
XX = (rNorm)'*sind(theta);
h = pcolor(XX,YY,Z,'parent',cax);
if ~isempty(ncolor)
    cmap = feval(myColorMap,ncolor);
    colormap(gca,cmap);
else
    colormap(gca,myColorMap);
end
% disp([max(R/Rrange),max(rNorm)])
shading flat
set(cax,'dataaspectratio',[1 1 1]);axis off;
if ~ishold(cax);
    % make a radial grid
    hold(cax,'on')
    % Draw circles and spokes
    createSpokes(thetaMin,thetaMax,Ncircles,circPos,Nspokes);
    createCircles(rMin,rMax,thetaMin,thetaMax,Ncircles,circPos,Nspokes)
end
%% PLot colorbar if specified
if colBar==1,
    c =colorbar('location','WestOutside');
    caxis([quantile(Z(:),0.01),quantile(Z(:),0.99)])
else
    c = [];
end
%% Outputs
nargoutchk(0,2)
if nargout==1,
    varargout{1}=h;
elseif nargout==2,
    varargout{1}=h;
    varargout{2}=c;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    
    
    
function createSpokes(thetaMin,thetaMax,Ncircles,circlesPos,Nspokes)
    
    spokeMesh = round(linspace(thetaMin,thetaMax,Nspokes));
    if isempty(circlesPos)
        circleMesh = linspace(rMin,rMax,Ncircles);
    else
        circleMesh  =  circlesPos;
    end
    contourD = abs((circleMesh - circleMesh(1))/Rrange+R(1)/Rrange);
    
    if strcmpi(typeRose,'meteo')
        cost = cosd(90-spokeMesh); % the zero angle is aligned with North
        sint = sind(90-spokeMesh); % the zero angle is aligned with North
    elseif strcmpi(typeRose,'default')
        cost = cosd(spokeMesh); % the zero angle is aligned with east
        sint = sind(spokeMesh); % the zero angle is aligned with east
    else
        error('"type" must be "meteo" or "default" ');
    end
    
    for kk = 1:Nspokes
        
        X = cost(kk)*contourD;
        Y = sint(kk)*contourD;
        
        if  Origin==0
            X(1)=Origin;
            Y(1)=Origin;
        end
        plot(X,Y,'color',[0.5,0.5,0.5],'linewidth',0.75,...
            'handlevisibility','off');
        % plot graduations of angles
        % avoid superimposition of 0 and 360
        if and(thetaMin==0,thetaMax == 360),
            if spokeMesh(kk)<360,
                
                text(1.05.*contourD(end).*cost(kk),...
                    1.05.*contourD(end).*sint(kk),...
                    [num2str(spokeMesh(kk),3),char(176)],...
                    'horiz', 'center', 'vert', 'middle');
            end
        else
            text(1.05.*contourD(end).*cost(kk),...
                1.05.*contourD(end).*sint(kk),...
                [num2str(spokeMesh(kk),3),char(176)],...
                'horiz', 'center', 'vert', 'middle');
        end
        
    end
end
    
    
    
    
    
    
    
    
    
    
function createCircles(rMin,rMax,thetaMin,thetaMax,Ncircles,circlePos,Nspokes)
    
    if isempty(circlePos)
        if Origin ==0 % if the origin is set at rMin
            contourD = linspace(0,1+R(1)/Rrange,Ncircles);
        else % if the origin is automatically centered at 0
            contourD = linspace(0,1,Ncircles)+R(1)/Rrange;
        end
    else
        
        contourD = circlePos-circlePos(1);
        contourD = contourD./max(contourD)*max(R/Rrange);
        contourD =[contourD(1:end-1)./contourD(end),1]+R(1)/Rrange;
    end
    
    if isempty(circlePos)
        if strcmpi(Rscale,'linear')||strcmpi(Rscale,'lin'),
            tickMesh = linspace(rMin,rMax,Ncircles);
        elseif strcmpi(Rscale,'log')||strcmpi(Rscale,'logarithmic'),
            tickMesh  = logspace(log10(rMin),log10(rMax),Ncircles);
        else
            error('''Rscale'' must be ''log'' or ''linear'' ');
        end
    else
        tickMesh  = circlePos;
        Ncircles = numel(tickMesh);
    end
    
    % define the grid in polar coordinates
    
    
    if strcmpi(typeRose,'meteo')
        angleGrid = linspace(90-thetaMin,90-thetaMax,100);
    elseif strcmpi(typeRose,'default')
        angleGrid = linspace(thetaMin,thetaMax,100);
    else
        error('"type" must be "meteo" or "default" ');
    end
    
    xGrid = cosd(angleGrid);
    yGrid = sind(angleGrid);
    spokeMesh = linspace(thetaMin,thetaMax,Nspokes);
    
    % plot circles
    for kk=1:length(contourD)
        X = xGrid*contourD(kk);
        Y = yGrid*contourD(kk);
        plot(X,Y,'color',[0.5,0.5,0.5],'linewidth',1);
    end
    % radius tick label
    
    position = 0.51.*(spokeMesh(min(Nspokes,round(Ncircles/2)))+...
        spokeMesh(min(Nspokes,1+round(Ncircles/2))));
    if strcmpi(typeRose,'meteo'),position = 90-position; end
    if strcmpi(typeRose,'default') && min(90-theta)<5,position = 0; end
    if min(round(theta))==90 && strcmpi(typeRose,'meteo'),  position = 0; end
    if max(round(theta))==90 && strcmpi(typeRose,'meteo'),  position = 0; end
    
    for kk=1:Ncircles
        if isempty(RtickLabel),
            rtick = num2str(tickMesh(kk),2);
        else
            rtick = RtickLabel(kk);
        end
        
        % radial graduations
        t = text(contourD(kk).*cosd(position),...
            (contourD(kk)).*sind(position),...
            rtick,'verticalalignment','BaseLine',...
            'horizontalAlignment', 'right',...
            'handlevisibility','off','parent',cax);
        if min(round(abs(90-theta)))<5 && strcmpi(typeRose,'default'),
            t.Position =  t.Position - [0,0.1,0];
            t.Interpreter = 'latex';
            clear t;
        end
        if min(round(theta))==90 && strcmpi(typeRose,'meteo')
            t.Position =  t.Position + [0,0.02,0];
            t.Interpreter = 'latex';
            clear t;
        elseif max(round(theta))==90 && strcmpi(typeRose,'meteo')
            t.Position =  t.Position - [0,0.05,0];
            t.Interpreter = 'latex';
            clear t;
        end
        
        % annotate spokes
        if max(theta)-min(theta)>180,
            t = text(contourD(end).*1.3.*cosd(position),...
                contourD(end).*1.3.*sind(position),...
                [labelR],'verticalalignment','bottom',...
                'horizontalAlignment', 'right',...
                'handlevisibility','off','parent',cax);
        else
            t = text(contourD(end).*0.6.*cosd(position),...
                contourD(end).*0.6.*sind(position),...
                [labelR],'verticalalignment','bottom',...
                'horizontalAlignment', 'right',...
                'handlevisibility','off','parent',cax);
        end
        
        t.Interpreter = 'latex';
        if min(round(theta))==90 && strcmpi(typeRose,'meteo'),
            t.Position =  t.Position + [0,0.05,0];
            clear t;
        elseif max(round(theta))==90 && strcmpi(typeRose,'meteo'),
            t.Position =  t.Position + [0,0.05,0];
            clear t;
        end
        %                 if min(round(abs(90-theta)))<5 && strcmpi(typeRose,'default'),
        %                     t.Position =  t.Position - [0,0.12,0];
        %                     t.Interpreter = 'latex';
        %                     clear t;
        %                 end
    end
    
end




function [rNorm] = getRnorm(Rscale,Origin,R,Rrange)
        if strcmpi(Rscale,'linear')||strcmpi(Rscale,'lin')
            rNorm = R-R(1)+Origin;
            rNorm = (rNorm)/max(rNorm)*max(R/Rrange);
        elseif strcmpi(Rscale,'log')||strcmpi(Rscale,'logarithmic')
            if rMin<=0
                error(' The radial vector cannot be lower or equal to 0 if the logarithmic scale is used');
            end
            rNorm = log10(R); %normalized radius [0,1]
            rNorm =rNorm-rNorm(1);
            rNorm = (rNorm)/max(rNorm)*max(R/Rrange);
        else
            error('''Rscale'' must be ''log'' or ''linear'' ');
        end
    end
end
