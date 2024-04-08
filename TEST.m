%I am going to use this file to try and create a 3D model of my project.
%date this file was started - 19/02/2024


%Defining constants in cgs units

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
M_acc_rate = ((1e-8) * M_sol) / year;

dS = 1e-3;
N = 1/dS;


B_octupole = 3000;
B_dipole = 500;
B = B_octupole/B_dipole;


theta = pi/2;
theta_dip = 0;
theta_oct = 0;    %For Parallel magnetic field

R_dip = [3.5; 5];


[x_coords1, z_coords1, B1] = dipole(R_dip(1), theta, B_dipole, dS, G, R_star, M_star);                %R_dip or R_Mult?
[x_coords2, z_coords2, B2] = dipole(R_dip(2), theta, B_dipole, dS, G, R_star, M_star);




figure(1)

ax = gca;
hold on
plot(x_coords1, z_coords1,'k--', LineWidth=1);
plot(x_coords2, z_coords2,'b--', LineWidth=1);
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
axis equal
grid on
xlabel('X axis')
ylabel('Z axis')
title('Test Dipole')
xlim(ax, [0 5.5])
ylim(ax, [0 5.5])
hold on

[x_star, y_star] = circle(0, 0, 1, N);

plot(x_star, y_star, '-', 'Color', [0.91 0.41 0.17], LineWidth=1)
legend('R1 Dipole', 'R2 Dipole', 'Star', Location='northeast')

hold on




Y1 = zeros(numel(x_coords1), 1);
Y2 = zeros(numel(x_coords2), 1);

x_values1 = [x_coords1, -x_coords1, x_coords1, -x_coords1];
z_values1 = [z_coords1, z_coords1, -z_coords1, -z_coords1]; 
y_values1 = [Y1, Y1, Y1, Y1];

x_values2 = [x_coords2, -x_coords2, x_coords2, -x_coords2];
z_values2 = [z_coords2, z_coords2, -z_coords2, -z_coords2]; 
y_values2 = [Y2, Y2, Y2, Y2];








B1_values = [B1(end), -B1(end), B1(end), -B1(end)] / 1000;
B2_values = [B2(end), B2(end), -B2(end), -B2(end)] / 1000;


B = [B1_values; B2_values];
M = zeros(101,101);






for i = 1:4

    for j = 1:2

        M(i,j) = B(j,i);

    end

end



[x,y,z] = sphere(100);


figure(2)


hs1 = surf(x,y,z,'EdgeColor','none');
hs1 = colormap("turbo");
hs1 = colorbar;

axis equal
hold on
plot3(x_values1, y_values1, z_values1, 'k--', LineWidth=1);
plot3(x_values2, y_values2, z_values2, 'b--', LineWidth=1);



















%Function to plot the circle

function [x_units, y_units] = circle(x, y, r, N)
    
    angle = 0:(pi/N):(2*pi);
    x_units = x + (r * cos(angle));
    y_units = y + (r * sin(angle));

end






function [x_coords, z_coords, B_r_array, theta1, R] = dipole(r, theta, B_dipole, dS,  G, R_star, M_star)


    x_coords = [r.*sin(theta)];
    z_coords = [r.*cos(theta)];
    B_r_array = [];
    R = [];


    while r > 1
    
    
    
        B_r = B_dipole .* (1./r).^3 .* cos(theta) .* cos(0);
        B_theta = 0.5 .* B_dipole .* (1./r).^3 .* sin(theta) .* cos(0);
        
    
        B = (B_r.^2 + B_theta.^2).^0.5;
        dr = (B_r ./ B) .* dS;
        dtheta = (B_theta ./ (r.*B)) .* dS;
    
        
        r = r - dr;
        theta = theta - dtheta;

        x_coords = [x_coords; r.*sin(theta)];
        z_coords = [z_coords; r.*cos(theta)];

        
        B_r_array = [B_r_array; B];

        R = [R; r];

    
    end

    theta1 = theta;

end