%In this code I will produce the second stage of the multipole graph, i.e.
%   the plots of the dipole + octupole magnetic field components

%Firstly, I will define some constants in cgs units

R_sol = 6.957 * 10^10; 
M_sol = 1.989 * 10^33;
G = 6.674 * 10^-8;
year = 3.1536 * 10^7;
day = 86400;

beta = 0.6;  %Beta Parameter

%Now, I will define values for a typical test stars

M_star = M_sol;
R_star = 2 * R_sol;       
Period = 7 * day;
M_acc_rate = ((10^-8) * M_sol) / year;

B_dipole = 1000;   %Dipole magnetic field strength
B_octupole = 4000;  %Octupole magnetic field strength

u = (B_dipole) * ((R_star^3)/2);   %NOT SURE ABOUT THIS

%Now, I will define the truncation radius for this test star, which will
%   also be the largest possible value of rm

R_t = (beta * (u^(4/7)) * ((2*G*M_star)^(-1/7)) * (M_acc_rate^(-2/7))) / R_sol;

   


%Now, I can start creating the plot of the dipole magnetic field case for
%   different values of rm

% B_r/dr = B_theta/dtheta = B/dS
% B = (B_r^2 + B_theta^2)^0.5

%For the axisymmetric multipole -- r = rm when theta = pi/2
%Hence, theta will range from 0 to pi/2

%The star will be plotted as a unit circle, hence the coordinates x and z
%   in cartesian coordinates will be:
%   x(j) = r(j)*sin(theta(j));            
%   z(j) = r(j)*cos(theta(j));

%r = rm * sin(theta)^2;

dS = 10^-3;
N = 1/dS;

%rm = linspace(1, R_t, N);
%theta = linspace(0, (pi/2), N);

% beta_dipole = 0 degrees

r_null = (3/4 * B_octupole/B_dipole)^0.5;  %In solar radii

rm = linspace (1.5, 6, 8);

col = ['b', 'y', 'g', 'r', 'k'];

%In case of an parallel magnetic field     

% for i = 1:numel(rm)
% 
%     r = rm(i);          %This value is given in R_stars
%     theta = pi/2;
%     
%     x_values = zeros(1, 1);
%     z_values = zeros(1, 1);
% 
% 
%         for j = 1:10000
%     
%             B_r = (B_dipole * ((1/r)^3) * cos(theta) * cos(0)) + 0.5 * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(0)^2) - 3) * cos(theta) * cos(0));
%             B_theta = (0.5 * B_dipole * ((1/r)^3) * sin(theta) * cos(0)) + (3/8) * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(0)^2) - 1) * sin(theta) * cos(0)); 
% 
%             B = (B_r^2 + B_theta^2)^0.5;
% 
%             dr = (B_r / B) * dS;
%             dtheta = (B_theta / (r*B)) * dS;
% 
%             r = r - dr;
%             theta = theta - dtheta;
% 
%             x = r*sin(theta);      %This converts the values to cartesian coordinates
%             z = r*cos(theta);
% 
%             if r > 1
% 
%                 x_values(j) = x;
%                 z_values(j) = z;
%             
%             end
%     
% 
%         end
%             
%     ax = gca;
%     hold on
%     plot(x_values, z_values,'b-', LineWidth=1);
%     ax.XAxisLocation = 'origin';
%     ax.YAxisLocation = 'origin';
%     axis equal
%     grid on
%     xlabel('X axis')
%     ylabel('Z axis')
%     title('Dipole Closed Magnetic Field Lines')
%     xlim(ax, [0 5.5])
%     ylim(ax, [-5.5 5.5])
%     hold on
% 
% 
% end


B_dipole = 1000;
B_octupole = 5000;

r_null = (3/4 * B_octupole/B_dipole)^0.5  %In solar radii

rm = linspace (1.1, 6, 8);
rm2 = linspace (1.1, r_null, 8);

theta = pi/2;
theta_dip = 0;
theta_oct = 0;    %For Parallel magnetic field
dS = 10e-3;

[x,y,z] = sphere;

[x_star, y_star] = circle(0, 0, 1, N);

[x1, z1, x_analytic, z_analytic] = dipole(rm(2), rm(2), theta, B_dipole, dS)
[x2, z2, x_analytic, z_analytic] = dipole(rm(5), rm(5), theta, B_dipole, dS)
[x3, z3, x_analytic, z_analytic] = dipole(rm(3), rm(3), theta, B_dipole, dS)

for i = 1:numel(rm)
    

    [x_coords, z_coords] = multipole(rm(i), theta, B_dipole, B_octupole, theta_dip, theta_oct, dS);
    [x_coords2, z_coords2] = multipole2(rm2(i), B_dipole, B_octupole, theta_dip, theta_oct, dS);
    [x_coords3, z_coords3] = multipole3(rm2(i), B_dipole, B_octupole, theta_dip, theta_oct, dS);

    
    if z_coords(i+1) < 0

        z_coords = z_coords .* -1;

    end
    
%     figure(1)
    subplot(1,2,1)
    ax = gca;
    hold on
    plot(x_coords, z_coords,'k-', LineWidth=1);
    plot(x_coords2, z_coords2,'k-', LineWidth=1);
    plot(x_coords3, z_coords3,'k-', LineWidth=1);

    plot(x1, z1,'r--', LineWidth=1);
    plot(x2, z2,'r--', LineWidth=1);
    plot(x3, z3,'r--', LineWidth=1);

    plot(x_star, y_star, '-', 'Color', [0.91 0.41 0.17], LineWidth=1)
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    axis equal
    grid on
    xlabel('X axis', 'FontSize', 15)
    ylabel('Z axis', 'FontSize', 15)
    ax.FontSize = 15;
    title('Parallel Multipole Closed Magnetic Field Lines', 'FontSize', 18)
    xlim(ax, [0 5.5])
    ylim(ax, [0 5.5])
    hold on



    
    X3D = [x_coords, -x_coords, x_coords, -x_coords];
    Y3D = zeros(numel(x_coords), 4);
    Z3D = [z_coords, z_coords, -z_coords, -z_coords];

    X3D_2 = [x_coords2, -x_coords2, x_coords2, -x_coords2];
    Y3D_2 = zeros(numel(x_coords2), 4);
    Z3D_2 = [z_coords2, z_coords2, -z_coords2, -z_coords2];

    X3D_3 = [x_coords3, -x_coords3, x_coords3, -x_coords3];
    Y3D_3 = zeros(numel(x_coords3), 4);
    Z3D_3 = [z_coords3, z_coords3, -z_coords3, -z_coords3];

    %figure(2)
%     subplot(2,2,3)
%     hs1 = surf(x,y,z);
%     q1 = get(hs1);
%     set(hs1, 'FaceColor', [0.91 0.41 0.17])
%     plot3(X3D, Y3D, Z3D, 'k--', LineWidth=1);
%     plot3(X3D_2, Y3D_2, Z3D_2, 'b--', LineWidth=1);
%     plot3(X3D_3, Y3D_3, Z3D_3, 'r--', LineWidth=1);
% 
%     xlabel('X axis')
%     ylabel('Y axis')
%     zlabel('Z axis')
%     title('Parallel Multipole Closed Magnetic Field Lines 3D')
%     grid on
%     axis equal
%     hold on







end







%To plot the star I will implement a function at the bottom of the file

% [x_star, y_star] = circle(0, 0, 1, N);
% 
% plot(x_star, y_star, '-', 'Color', [0.91 0.41 0.17], LineWidth=1)



% The next two lines are used in case the star needs to be filled
d = [0.91 0.41 0.17];
fill(x_star, y_star, d)



%For anti-parallel magnetic field
theta_oct = 180;

for i = 1:numel(rm)
    

    [x_coords, z_coords] = multipole(rm(i), theta, B_dipole, B_octupole, theta_dip, theta_oct, dS);
    [x_coords2, z_coords2] = multipole2(rm2(i), B_dipole, B_octupole, theta_dip, theta_oct, dS);
    [x_coords3, z_coords3] = multipole3(rm2(i), B_dipole, B_octupole, theta_dip, theta_oct, dS);
    
    if z_coords(i+1) < 0

        z_coords = z_coords .* -1;

    end
    
%     figure(2)
    subplot(1,2,2)
    ax = gca;
    hold on
    plot(x_coords, z_coords,'k-', LineWidth=1);
    plot(x_coords2, z_coords2,'k-', LineWidth=1);
    plot(x_coords3, z_coords3,'k-', LineWidth=1);

    plot(x1, z1,'r--', LineWidth=1);
    plot(x2, z2,'r--', LineWidth=1);
    plot(x3, z3,'r--', LineWidth=1);

    plot(x_star, y_star, '-', 'Color', [0.91 0.41 0.17], LineWidth=1)
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    axis equal
    grid on
    xlabel('X axis', 'FontSize', 15)
    ylabel('Z axis', 'FontSize', 15)
    title('Anti-parallel Multipole Closed Magnetic Field Lines', 'FontSize', 18)
    ax.FontSize = 15;
    xlim(ax, [0 5.5])
    ylim(ax, [0 5.5])
    hold on




    X3D = [x_coords, -x_coords, x_coords, -x_coords];
    Y3D = zeros(numel(x_coords), 4);
    Z3D = [z_coords, z_coords, -z_coords, -z_coords];

    X3D_2 = [x_coords2, -x_coords2, x_coords2, -x_coords2];
    Y3D_2 = zeros(numel(x_coords2), 4);
    Z3D_2 = [z_coords2, z_coords2, -z_coords2, -z_coords2];

    X3D_3 = [x_coords3, -x_coords3, x_coords3, -x_coords3];
    Y3D_3 = zeros(numel(x_coords3), 4);
    Z3D_3 = [z_coords3, z_coords3, -z_coords3, -z_coords3];

    %figure(2)
%     subplot(2,2,4)
%     hs1 = surf(x,y,z);
%     q1 = get(hs1);
%     set(hs1, 'FaceColor', [0.91 0.41 0.17])
%     plot3(X3D, Y3D, Z3D, 'k--', LineWidth=1);
%     plot3(X3D_2, Y3D_2, Z3D_2, 'b--', LineWidth=1);
%     plot3(X3D_3, Y3D_3, Z3D_3, 'r--', LineWidth=1);
% 
%     xlabel('X axis')
%     ylabel('Y axis')
%     zlabel('Z axis')
%     title('Anti-parallel Multipole Closed Magnetic Field Lines 3D')
%     grid on
%     axis equal
%     hold on


end


d = [0.91 0.41 0.17];
fill(x_star, y_star, d)






%To plot the star I will implement a function at the bottom of the file

% [x_star, y_star] = circle(0, 0, 1, N);
% 
% plot(x_star, y_star, '-', 'Color', [0.91 0.41 0.17], LineWidth=1)

hold off


% The next two lines are used in case the star needs to be filled
% d = [0.91 0.41 0.17];
% fill(x_star, y_star, d)









%Function to plot the circle

function [x_units, y_units] = circle(x, y, r, N)
    
    angle = 0:(pi/N):(2*pi);
    x_units = x + (r * cos(angle));
    y_units = y + (r * sin(angle));

end




%This is the function to plot the multipole coordinates, I will try to
%integrate the closed lines given starting from the theta=0 line.

function [x_coords, z_coords] = multipole(r, theta, B_dipole, B_octupole, theta_dip, theta_oct, dS)


    x_coords = [r.*sin(theta)];
    z_coords = [r.*cos(theta)];


    while r > 1
    
    
    
        B_r = (B_dipole * ((1/r)^3) * cos(theta) * cos(theta_dip)) + 0.5 * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(0)^2) - 3) * cos(theta) * cos(theta_oct));
        B_theta = (0.5 * B_dipole * ((1/r)^3) * sin(theta) * cos(theta_dip)) + (3/8) * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(0)^2) - 1) * sin(theta) * cos(theta_oct)); 
        
    
        B = (B_r.^2 + B_theta.^2).^0.5;
        dr = (B_r ./ B) .* dS;
        dtheta = (B_theta ./ (r.*B)) .* dS;
    
        
        r = r - dr;
        theta = theta - dtheta;



        x_coords = [x_coords; r.*sin(theta)];
        z_coords = [z_coords; r.*cos(theta)];
        


    
    end


end







%Below, the two functions to plot the multipole component inside the null
%radius 




function [x_coords2, z_coords2] = multipole2(r, B_dipole, B_octupole, theta_dip, theta_oct, dS)

    theta = acos((3/5 + (2/5) * (B_dipole/B_octupole) * (1/r)^2)^0.5);

    x_coords2 = [r.*sin(theta)];
    z_coords2 = [r.*cos(theta)];


    while r > 1
    
    
    
        B_r = (B_dipole * ((1/r)^3) * cos(theta) * cos(theta_dip)) + 0.5 * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(0)^2) - 3) * cos(theta) * cos(theta_oct));
        B_theta = (0.5 * B_dipole * ((1/r)^3) * sin(theta) * cos(theta_dip)) + (3/8) * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(0)^2) - 1) * sin(theta) * cos(theta_oct)); 
        
    
        B = (B_r.^2 + B_theta.^2).^0.5;
        dr = (B_r ./ B) .* dS;
        dtheta = (B_theta ./ (r.*B)) .* dS;
    
        
        r = r + dr;
        theta = theta + dtheta;



        x_coords2 = [x_coords2; r.*sin(theta)];
        z_coords2 = [z_coords2; r.*cos(theta)];
        


    
    end

    


end







function [x_coords3, z_coords3] = multipole3(r, B_dipole, B_octupole, theta_dip, theta_oct, dS)

    theta = acos((3/5 + (2/5) * (B_dipole/B_octupole) * (1/r)^2)^0.5);

    x_coords3 = [r.*sin(theta)];
    z_coords3 = [r.*cos(theta)];


    while r > 1
    
    
    
        B_r = (B_dipole * ((1/r)^3) * cos(theta) * cos(theta_dip)) + 0.5 * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(0)^2) - 3) * cos(theta) * cos(theta_oct));
        B_theta = (0.5 * B_dipole * ((1/r)^3) * sin(theta) * cos(theta_dip)) + (3/8) * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(0)^2) - 1) * sin(theta) * cos(theta_oct)); 
        
    
        B = (B_r.^2 + B_theta.^2).^0.5;
        dr = (B_r ./ B) .* dS;
        dtheta = (B_theta ./ (r.*B)) .* dS;
    
        
        r = r - dr;
        theta = theta - dtheta;



        x_coords3 = [x_coords3; r.*sin(theta)];
        z_coords3 = [z_coords3; r.*cos(theta)];
        


    
    end

    


end
















function [x_coords, z_coords, x_analytic, z_analytic] = dipole(r, r_check, theta, B_dipole, dS)


    x_coords = [r.*sin(theta)];
    z_coords = [r.*cos(theta)];

    x_analytic = [r_check.*sin(theta)];
    z_analytic = [r_check.*cos(theta)];


    while r > 1
    
    
    
        B_r = B_dipole .* (1./r).^3 .* cos(theta) .* cos(0);
        B_theta = 0.5 .* B_dipole .* (1./r).^3 .* sin(theta) .* cos(0);
        
    
        B = (B_r.^2 + B_theta.^2).^0.5;
        dr = (B_r ./ B) .* dS;
        dtheta = (B_theta ./ (r.*B)) .* dS;
    
        
        r = r - dr;
        theta = theta - dtheta;

        r_analytic = r_check * (sin(theta)^2);



        x_coords = [x_coords; r.*sin(theta)];
        z_coords = [z_coords; r.*cos(theta)];
        
        x_analytic = [x_analytic; r_analytic.*sin(theta)];
        z_analytic = [z_analytic; r_analytic.*cos(theta)];

    
    end


end


