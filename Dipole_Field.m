%In this code I will produce the first stage of the multipole graph, i.e.
%   the plots of the dipole magnetic field component

%Firstly, I will define some constants in cgs units

R_sol = 6.957e10; 
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

u = B_dipole * ((R_star^3)/2);

%Now, I will define the truncation radius for this test star, which will
%   also be the largest possible value of rm

R_t = (beta * (u^(4/7)) * ((2*G*M_star)^(-1/7)) * (M_acc_rate^(-2/7))) / R_sol;

   


%Now, I can start creating the plot of the dipole magnetic field case for
%   different values of rm

% B_r/dr = B_theta/dtheta = B/dS
% B = (B_r^2 + B_theta^2)^0.5

%For the axisymmetric dipole -- r = rm when theta = pi/2
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

rm = linspace (1.5, 6, 8);
theta = pi/2;


% for i = 1:numel(rm)
% 
%     r = rm(i);
%     r_check = rm(i);
% 
%     theta = pi/2;
%     
%     x_values = zeros(1, 1);
%     z_values = zeros(1, 1);
% 
%     x_check = zeros(1,1);
%     z_check = zeros(1,1);
% 
% 
%         for j = 1:10000
%     
%             B_r = B_dipole * (R_star/r)^3 * cos(theta) * cos(0);
%             B_theta = 0.5 * B_dipole * ((R_star/r)^3) * sin(theta) * cos(0);
%             %B_r = B_dipole * (1/r)^3 * cos(theta) * cos(0);
%             %B_theta = 0.5 * B_dipole * (1/r)^3 * sin(theta) * cos(0);
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
% 
%             r_check = rm(i) * (sin(theta)^2);
%             xx = r_check * sin(theta);
%             zz = r_check * cos(theta);
% 
% 
%             if r > 1
% 
%                 x_values(j) = x;
%                 z_values(j) = z;
% 
%             end
%     
%             if r_check > 1
% 
%                 x_check(j) = xx;
%                 z_check(j) = zz;
% 
%             end
% 
%         end
%             
%     ax = gca;
%     hold on
%     plot(x_values, z_values,'k--', LineWidth=1);
%     plot(x_check, z_check,'b-', LineWidth=1);
%     ax.XAxisLocation = 'origin';
%     ax.YAxisLocation = 'origin';
%     axis equal
%     grid on
%     xlabel('X axis')
%     ylabel('Z axis')
%     title('Dipole Closed Magnetic Field Lines')
%     xlim(ax, [0 5.5])
%     ylim(ax, [0 5.5])
%     hold on
% 
% 
% end


% Now I will create a while loop version of the above code


[x,y,z] = sphere;

rm = linspace (1.5, 6, 8);
theta = pi/2;
B_dipole = 1000;
dS = 10e-3;


[x_star, y_star] = circle(0, 0, 1, N);

%figure(1)
plot(x_star, y_star, '-', 'Color', [0.91 0.41 0.17], LineWidth=1)

for i = 1:numel(rm)
    

    [x_coords, z_coords, x_analytic, z_analytic] = dipole(rm(i), rm(i), theta, B_dipole, dS);

    figure(1)
%     subplot(1,2,1)
    ax = gca;
    hold on
    plot(x_coords, z_coords,'k-', LineWidth=1);
    plot(x_star, y_star, '-', 'Color', [0.91 0.41 0.17], LineWidth=1)
    %plot(x_analytic, z_analytic,'b-', LineWidth=1);
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    axis equal
    grid on
    ax.FontSize = 15;
    xlabel('X axis')
    ylabel('Z axis')
    title('Dipole Closed Magnetic Field Lines')
    xlim(ax, [0 5.5])
    ylim(ax, [0 3])
    hold on
    d = [0.91 0.41 0.17];
    fill(x_star, y_star, d)

    X3D = [x_coords, -x_coords, x_coords, -x_coords];
    Y3D = zeros(numel(x_coords), 4);
    Z3D = [z_coords, z_coords, -z_coords, -z_coords];
    

%     figure(2)
% %     subplot(1,2,2)
% 
%     hs1 = surf(x,y,z,'EdgeColor','none');
%     hs1 = colormap("turbo");
%     hs1 = colorbar;
% %     hs1 = surf(x,y,z);
% %     q1 = get(hs1);
% %     set(hs1, 'FaceColor', [0.91 0.41 0.17])
% 
%     plot3(X3D, Y3D, Z3D, 'k--', LineWidth=1);
% 
%     xlabel('X axis')
%     ylabel('Y axis')
%     zlabel('Z axis')
%     title('Dipole Closed Magnetic Field Lines 3D')
%     grid on
%     axis equal
%     hold on


end


% %To plot the star I will implement a function at the bottom of the file
% 
% [x_star, y_star] = circle(0, 0, 1, N);
% 
% plot(x_star, y_star, '-', 'Color', [0.91 0.41 0.17], LineWidth=1)


% The next two lines are used in case the star needs to be filled
d = [0.91 0.41 0.17];
fill(x_star, y_star, d)

%legend({['rm = ', num2str(rm(1))], ['rm = ', num2str(rm(2))], ['rm = ', num2str(rm(3))], ['rm = ', num2str(rm(4))], ['rm = ', num2str(rm(end))], 'Star'},'Location','northwest','Orientation','vertical')



%LegendeP (l, cos(theta))
%Rloop = [rloop; r]




%Function to plot the circle

function [x_units, y_units] = circle(x, y, r, N)
    
    angle = 0:(pi/N):(2*pi);
    x_units = x + (r * cos(angle));
    y_units = y + (r * sin(angle));

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


