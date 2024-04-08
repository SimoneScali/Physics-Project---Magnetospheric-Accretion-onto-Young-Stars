%In this code I will produce the first stage of the octupole graph, i.e.
%   the plots of the octupole magnetic field components

%Firstly, I will define some constants in cgs units

R_sol = 6.957e10; 
M_sol = 1.989e33;
G = 6.674e-8;
year = 3.1536e7;
day = 86400;

beta = 0.6;  %Beta Parameter

%Now, I will define values for a typical test stars

M_star = M_sol;
R_star = 2 * R_sol;       
Period = 7 * day;
M_acc_rate = ((10^-8) * M_sol) / year;

B_octupole = 4000;  %Octupole magnetic field strength

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





B_octupole = 3000;


rm = linspace (1.1, 6, 8);

theta = pi/2;
theta2 = acos((3/5)^0.5);

x_line = zeros(1,10);
z_line = zeros(1,10);

theta_oct = 0;
dS = 10e-3;

[x,y,z] = sphere;

[x_star, y_star] = circle(0, 0, 1, N);


for i = 1:numel(rm)
    

    [x_coords, z_coords] = octupole(rm(i), rm(i), theta, B_octupole, theta_oct, dS);
    [x_coords2, z_coords2] = octupole2(rm(i), rm(i), theta2, B_octupole, theta_oct, dS);
    [x_coords3, z_coords3] = octupole3(rm(i), rm(i), theta2, B_octupole, theta_oct, dS);

    x_line(i) = x_coords3(1);
    z_line(i) = z_coords3(1);
    
    if z_coords(i+1) < 0

        z_coords = z_coords .* -1;

    end

    figure(1)
%     subplot(1,2,1)
    ax = gca;
    hold on
    plot(x_coords, z_coords,'k-', LineWidth=1);
    plot(x_coords2, z_coords2,'k-', LineWidth=1);
    plot(x_coords3, z_coords3,'k-', LineWidth=1);
    plot(x_line, z_line, 'r-', LineWidth=1);
    plot(x_star, y_star, '-', 'Color', [0.91 0.41 0.17], LineWidth=1)
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    axis equal
    ax.FontSize = 15;
    grid on
    xlabel('X axis')
    ylabel('Z axis')
    title('Octupole Closed Magnetic Field Lines')
    xlim(ax, [0 5.5])
    ylim(ax, [0 5.5])
    hold on
    
    d = [0.91 0.41 0.17];
    fill(x_star, y_star, d)


    X3D = [x_coords, -x_coords, x_coords, -x_coords];
    Y3D = zeros(numel(x_coords), 4);
    Z3D = [z_coords, z_coords, -z_coords, -z_coords];

    X3D_2 = [x_coords2, -x_coords2, x_coords2, -x_coords2];
    Y3D_2 = zeros(numel(x_coords2), 4);
    Z3D_2 = [z_coords2, z_coords2, -z_coords2, -z_coords2];

    X3D_3 = [x_coords3, -x_coords3, x_coords3, -x_coords3];
    Y3D_3 = zeros(numel(x_coords3), 4);
    Z3D_3 = [z_coords3, z_coords3, -z_coords3, -z_coords3];

    figure(2)
%     subplot(1,2,2)

%     hs1 = surf(x,y,z);
%     q1 = get(hs1);
%     set(hs1, 'FaceColor', [0.91 0.41 0.17])

    hs1 = surf(x,y,z,'EdgeColor','none');
    hs1 = colormap("turbo");
    hs1 = colorbar;

    plot3(X3D, Y3D, Z3D, 'k--', LineWidth=1);
    plot3(X3D_2, Y3D_2, Z3D_2, 'b--', LineWidth=1);
    plot3(X3D_3, Y3D_3, Z3D_3, 'r--', LineWidth=1);

    xlabel('X axis')
    ylabel('Y axis')
    zlabel('Z axis')
    title('Octupole Closed Magnetic Field Lines 3D')
    grid on
    axis equal
    hold on

end







%To plot the star I will implement a function at the bottom of the file

% [x_star, y_star] = circle(0, 0, 1, N);
% 
% plot(x_star, y_star, '-', 'Color', [0.91 0.41 0.17], LineWidth=1)


% The next two lines are used in case the star needs to be filled
%d = [0.91 0.41 0.17];
%fill(x_star, y_star, d)







%Function to plot the circle

function [x_units, y_units] = circle(x, y, r, N)
    
    angle = 0:(pi/N):(2*pi);
    x_units = x + (r * cos(angle));
    y_units = y + (r * sin(angle));

end




%This is the function to plot the multipole coordinates, I will try to
%integrate the closed lines given starting from the theta=0 line.

function [x_coords, z_coords] = octupole(r, rm, theta, B_octupole, theta_oct, dS)


    x_coords = [r.*sin(theta)];
    z_coords = [r.*cos(theta)];


    while r > 1
    
    
    
        B_r = 0.5 * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(0)^2) - 3) * cos(theta) * cos(theta_oct));
        B_theta = (3/8) * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(0)^2) - 1) * sin(theta) * cos(theta_oct)); 
        
    
        B = (B_r.^2 + B_theta.^2).^0.5;
        dr = (B_r ./ B) .* dS;
        dtheta = (B_theta ./ (r.*B)) .* dS;
    
        
        r = rm * ((1 - 5*cos(theta)^2) * sin(theta)^2)^(1/3);
        theta = theta - dtheta;



        x_coords = [x_coords; r.*sin(theta)];
        z_coords = [z_coords; r.*cos(theta)];
        


    
    end


end





%Below, the two functions to plot the octupole component at high altitude




function [x_coords2, z_coords2] = octupole2(r, rm, theta, B_octupole, theta_oct, dS)


    x_coords2 = [r.*sin(theta)];
    z_coords2 = [r.*cos(theta)];


    while r > 1
    
    
    
        B_r = 0.5 * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(0)^2) - 3) * cos(theta) * cos(theta_oct));
        B_theta = (3/8) * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(0)^2) - 1) * sin(theta) * cos(theta_oct)); 
        
    
        B = (B_r.^2 + B_theta.^2).^0.5;
        dr = (B_r ./ B) .* dS;
        dtheta = (B_theta ./ (r.*B)) .* dS;
    
        
        r = rm * ((5/4)*(5*cos(theta)^2 - 1) * sin(theta)^2)^(1/3);
        theta = theta + dtheta;



        x_coords2 = [x_coords2; r.*sin(theta)];
        z_coords2 = [z_coords2; r.*cos(theta)];
        


    
    end


end







function [x_coords3, z_coords3] = octupole3(r, rm, theta, B_octupole, theta_oct, dS)


    x_coords3 = [r.*sin(theta)];
    z_coords3 = [r.*cos(theta)];


    while r > 1
    
    
    
        B_r = 0.5 * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(0)^2) - 3) * cos(theta) * cos(theta_oct));
        B_theta = (3/8) * (B_octupole * ((1/r)^5) * (5 * (cos(theta)^2) * (cos(0)^2) - 1) * sin(theta) * cos(theta_oct)); 
        
    
        B = (B_r.^2 + B_theta.^2).^0.5;
        dr = (B_r ./ B) .* dS;
        dtheta = (B_theta ./ (r.*B)) .* dS;
    
        
        r = rm * ((5/4)*(5*cos(theta)^2 - 1) * sin(theta)^2)^(1/3);
        theta = theta - dtheta;



        x_coords3 = [x_coords3; r.*sin(theta)];
        z_coords3 = [z_coords3; r.*cos(theta)];
        


    
    end


end

