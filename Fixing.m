%This is gonna be the final coding file where I will merge the previous
%ones to produce various plots such as the fall accretion rate and number
%density


% Firstly, I will calculate A_spot, the area of the sun spots in the star
% caused by different magnetic fields

%Defining constants in cgs units


clear all
clf


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

[x_star, y_star] = circle(0, 0, 1, N); %This will be used to plot the stars as unit circles


B_octupole = 5000;
B_dipole = 500;
B = B_octupole/B_dipole;


%there are 3 possible scenarios 
%1 -- R_in and R_out > r_null
%2 -- R_in < r_null < R_out
%3 -- R_in and R_out < r_null

r_null = (3/4 * B_octupole/B_dipole)^0.5  %In solar radii



%Firstly, this is for a parallel multipole magnetic field

% [y, x, r_a_dip, r_a] = multitrunc(B, B_dipole, R_star, M_star, M_acc_rate, G, R_sol);

% R_T_mult = r_a * beta;
% R_mult = [R_T_mult; R_T_mult + 1];


theta = pi/2;
theta_dip = 0;
theta_oct = [0; 180];    % 0 for Parallel magnetic field and 180 for anti parallel



    
[x_para_RT, z_para_RT, B_RT_para, theta1_para] = multipole(3, theta, B_dipole, B_octupole, theta_dip, theta_oct(1), dS, G, R_star, M_star);
[x_para_r, z_para_r, B_r_para, theta2_para, R_para] = multipole(4, theta, B_dipole, B_octupole, theta_dip, theta_oct(1), dS, G, R_star, M_star);

    
figure(1)
ax = gca;
hold on
plot(x_para_RT, z_para_RT,'r--', LineWidth=1);
plot(x_para_r, z_para_r,'r--', LineWidth=1);
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
axis equal
grid on
xlabel('X axis')
ylabel('Z axis')
title('Accretion Plotting')
xlim(ax, [0 5.5])
ylim(ax, [0 5.5])
hold on

plot(x_star, y_star, '-', 'Color', [0.91 0.41 0.17], LineWidth=1)

hold on


%Now, I will do the same for an anti-parallel multipole magnetic field

    
[x_anti_RT, z_anti_RT, B_RT_anti, theta1_anti] = multipole(3, theta, B_dipole, B_octupole, theta_dip, theta_oct(2), dS, G, R_star, M_star);
[x_anti_r, z_anti_r, B_r_anti, theta2_anti, R_anti] = multipole(4, theta, B_dipole, B_octupole, theta_dip, theta_oct(2), dS, G, R_star, M_star);

    
ax = gca;
hold on
plot(x_anti_RT, z_anti_RT,'b--', LineWidth=1);
plot(x_anti_r, z_anti_r,'b--', LineWidth=1);

hold on




%Now, I will plot the dipole


[x_dip_RT, z_dip_RT, B_RT_dip, theta1_dip] = dipole(3, theta, B_dipole, dS, G, R_star, M_star);                %R_dip or R_Mult?
[x_dip_r, z_dip_r, B_r_dip, theta2_dip, R_dip_final] = dipole(4, theta, B_dipole, dS, G, R_star, M_star);

ax = gca;
hold on
plot(x_dip_RT, z_dip_RT,'k--', LineWidth=1);
plot(x_dip_r, z_dip_r,'k--', LineWidth=1);
legend('R_T para', 'r para', 'Star', 'R_T anti', 'r anti', 'R_T dip', 'r dip', Location='northeast')

hold on




%Now I will plot the null radius curves for the parallel and antiparallel
%fields

[x_out_para, z_out_para, b, theta_null_out, r] = multipole(r_null+0.001, theta, B_dipole, B_octupole, theta_dip, theta_oct(1), dS, G, R_star, M_star);
[x_in_para, z_in_para, b, theta_null_in, r] = multipole(r_null-0.001, theta, B_dipole, B_octupole, theta_dip, theta_oct(1), dS, G, R_star, M_star);

[x_out_anti, z_out_anti] = multipole(r_null+0.01, theta, B_dipole, B_octupole, theta_dip, theta_oct(2), dS, G, R_star, M_star);
[x_in_anti, z_in_anti] = multipole(r_null-0.01, theta, B_dipole, B_octupole, theta_dip, theta_oct(2), dS, G, R_star, M_star);



figure(2)
subplot(1,2,1)
ax = gca;
hold on
plot(x_out_para, z_out_para,'r--', LineWidth=1);
plot(x_in_para, -z_in_para,'b--', LineWidth=1);
plot(x_star, y_star, '-', 'Color', [0.91 0.41 0.17], LineWidth=1)
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
axis equal
grid on
xlabel('X axis')
ylabel('Z axis')
title('Accretion Plotting for Null Radius (Parallel)')
xlim(ax, [0 5.5])
ylim(ax, [0 5.5])
hold on


subplot(1,2,2)
ax = gca;
hold on
plot(x_out_anti, z_out_anti,'r--', LineWidth=1);
plot(x_in_anti, z_in_anti,'b--', LineWidth=1);
plot(x_star, y_star, '-', 'Color', [0.91 0.41 0.17], LineWidth=1)
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
axis equal
grid on
xlabel('X axis')
ylabel('Z axis')
title('Accretion Plotting for Null Radius (Anti-Parallel)')
xlim(ax, [0 5.5])
ylim(ax, [0 5.5])
hold on



%Now I will make a subplot showing the 3 possible scenarios.

B_oct_area = linspace(0, 5000, 30);
B_dipole = 500;


RR1 = [4, 5; 2, 3; 1.5, 2];


for j = 2

    A_para = [];
    A_anti = [];
    A_dip = [];

    R_para = [];
    R_anti = [];
    R_dip = [];

    B_para = [];
    B_anti = [];
    B_dip = [];

    for i = 1:numel(B_oct_area)
    
    
    
        r_null2 = (3/4 * B_oct_area(i)/B_dipole)^0.5;
    
    
    
        if r_null2 < RR1(j,1)
    
            [x1_in_para, z1_in_para, B1_in_para, theta1_para1, R1_in_para] = multipole(RR1(2,1), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
            [x1_out_para, z1_out_para, B1_out_para, theta2_para1, R1_out_para] = multipole(RR1(2,2), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
        
            [x1_in_anti, z1_in_anti, B1_in_anti, theta1_anti1, R1_in_anti] = multipole(RR1(2,1), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(2), dS, G, R_star, M_star);
            [x1_out_anti, z1_out_anti, B1_out_anti, theta2_anti1, R1_out_anti] = multipole(RR1(2,2), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(2), dS, G, R_star, M_star);
        
            [xnull1_out_para, znull1_out_para, b, theta_null_out1, r] = multipole(r_null2+0.001, theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
            [xnull1_in_para, znull1_in_para, b, theta_null_in1, r] = multipole(r_null2-0.001, theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
       
    
    
    
    
            [x_in_para, z_in_para, B_in_para, theta1_para, R_in_para] = multipole(RR1(j,1), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
            [x_out_para, z_out_para, B_out_para, theta2_para, R_out_para] = multipole(RR1(j,2), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
        
            [x_in_anti, z_in_anti, B_in_anti, theta1_anti, R_in_anti] = multipole(RR1(j,1), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(2), dS, G, R_star, M_star);
            [x_out_anti, z_out_anti, B_out_anti, theta2_anti, R_out_anti] = multipole(RR1(j,2), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(2), dS, G, R_star, M_star);
        
            [xnull_out_para, znull_out_para, b, theta_null_out, r] = multipole(r_null2+0.001, theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
            [xnull_in_para, znull_in_para, b, theta_null_in, r] = multipole(r_null2-0.001, theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
    
    
            
            [x_in_dip, z_in_dip, B_in_dip, theta_in_dip, R_in_dip] = dipole(RR1(j,1), theta, B_dipole, dS,  G, R_star, M_star);
            [x_out_dip, z_out_dip, B_out_dip, theta_out_dip, R_out_dip] = dipole(RR1(j,2), theta, B_dipole, dS,  G, R_star, M_star);





            A_para1 = 4*pi*(cos(theta2_para) - cos(theta1_para)) * R_star^2; %% 4*pi*((cos(theta2_out_para) - cos(test1)) + ((cos(test2) - cos(theta2_in_para)))) * R_star^2;
            A_anti1 = 4*pi*(cos(theta2_anti) - cos(theta1_anti)) * R_star^2;
            A_dip1 = 4*pi*(cos(theta_out_dip) - cos(theta_in_dip)) * R_star^2;
    
            A_para = [A_para; A_para1];
            A_anti = [A_anti; A_anti1];
            A_dip = [A_dip; A_dip1];


    
        elseif RR1(j,1) < r_null2 && r_null2 < RR1(j,2)
    
            [x2_in_para, z2_in_para, B2_in_para, theta1_para2, R2_in_para] = multipole(RR1(3,1), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
            [x2_out_para, z2_out_para, B2_out_para, theta2_para2, R2_out_para] = multipole(RR1(2,2), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
        
            [x2_in_anti, z2_in_anti, B2_in_anti, theta1_anti2, R2_in_anti] = multipole(RR1(3,1), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(2), dS, G, R_star, M_star);
            [x2_out_anti, z2_out_anti, B2_out_anti, theta2_anti2, R2_out_anti] = multipole(RR1(2,2), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(2), dS, G, R_star, M_star);
        
            [xnull2_out_para, znull2_out_para, b, theta_null_out2, r] = multipole(r_null2+0.001, theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
            [xnull2_in_para, znull2_in_para, b, theta_null_in2, r] = multipole(r_null2-0.001, theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
    
    
            
    
    
    
            [x_in_para, z_in_para, B_in_para, theta1_para, R_in_para] = multipole(RR1(j,1), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
            [x_out_para, z_out_para, B_out_para, theta2_para, R_out_para] = multipole(RR1(j,2), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
        
            [x_in_anti, z_in_anti, B_in_anti, theta1_anti, R_in_anti] = multipole(RR1(j,1), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(2), dS, G, R_star, M_star);
            [x_out_anti, z_out_anti, B_out_anti, theta2_anti, R_out_anti] = multipole(RR1(j,2), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(2), dS, G, R_star, M_star);
        
            [xnull_out_para, znull_out_para, b, theta_null_out, r] = multipole(r_null2+0.001, theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
            [xnull_in_para, znull_in_para, b, theta_null_in, r] = multipole(r_null2-0.001, theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
    
    

            [x_in_dip, z_in_dip, B_in_dip, theta_in_dip, R_in_dip] = dipole(RR1(j,1), theta, B_dipole, dS,  G, R_star, M_star);
            [x_out_dip, z_out_dip, B_out_dip, theta_out_dip, R_out_dip] = dipole(RR1(j,2), theta, B_dipole, dS,  G, R_star, M_star);


            theta_null_in = pi - theta_null_in;
            theta1_para = pi - theta1_para;


            A_para1 = 4*pi*((cos(theta2_para) - cos(theta_null_out)) + ((cos(theta_null_in) - cos(theta1_para)))) * R_star^2;
            A_anti1 = 4*pi*(cos(theta2_anti) - cos(theta1_anti)) * R_star^2;
            A_dip1 = 4*pi*(cos(theta_out_dip) - cos(theta_in_dip)) * R_star^2;
    
            A_para = [A_para; A_para1];
            A_anti = [A_anti; A_anti1];
            A_dip = [A_dip; A_dip1];
    
    
    
    
        elseif RR1(j,2) < r_null2
    
            [x3_in_para, z3_in_para, B3_in_para, theta1_para3, R3_in_para] = multipole(RR1(3,1), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
            [x3_out_para, z3_out_para, B3_out_para, theta2_para3, R3_out_para] = multipole(RR1(3,2), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
        
            [x3_in_anti, z3_in_anti, B3_in_anti, theta1_anti3, R3_in_anti] = multipole(RR1(3,1), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(2), dS, G, R_star, M_star);
            [x3_out_anti, z3_out_anti, B3_out_anti, theta2_anti3, R3_out_anti] = multipole(RR1(3,2), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(2), dS, G, R_star, M_star);
        
            [xnull3_out_para, znull3_out_para, b, theta_null_out3, r] = multipole(r_null2+0.001, theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
            [xnull3_in_para, znull3_in_para, b, theta_null_in3, r] = multipole(r_null2-0.001, theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
        
    
    
    
    
            [x_in_para, z_in_para, B_in_para, theta1_para, R_in_para] = multipole(RR1(j,1), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
            [x_out_para, z_out_para, B_out_para, theta2_para, R_out_para] = multipole(RR1(j,2), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
        
            [x_in_anti, z_in_anti, B_in_anti, theta1_anti, R_in_anti] = multipole(RR1(j,1), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(2), dS, G, R_star, M_star);
            [x_out_anti, z_out_anti, B_out_anti, theta2_anti, R_out_anti] = multipole(RR1(j,2), theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(2), dS, G, R_star, M_star);
        
            [xnull_out_para, znull_out_para, b, theta_null_out, r] = multipole(r_null2+0.001, theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
            [xnull_in_para, znull_in_para, b, theta_null_in, r] = multipole(r_null2-0.001, theta, B_dipole, B_oct_area(i), theta_dip, theta_oct(1), dS, G, R_star, M_star);
    



            [x_in_dip, z_in_dip, B_in_dip, theta_in_dip, R_in_dip] = dipole(RR1(j,1), theta, B_dipole, dS,  G, R_star, M_star);
            [x_out_dip, z_out_dip, B_out_dip, theta_out_dip, R_out_dip] = dipole(RR1(j,2), theta, B_dipole, dS,  G, R_star, M_star);


    
            theta1_para = pi - theta1_para;
            theta2_para = pi - theta2_para;


            A_para1 = 4*pi*(cos(theta2_para) - cos(theta1_para)) * R_star^2; %% 4*pi*((cos(theta2_out_para) - cos(test1)) + ((cos(test2) - cos(theta2_in_para)))) * R_star^2;
            A_anti1 = 4*pi*(cos(theta2_anti) - cos(theta1_anti)) * R_star^2;
            A_dip1 = 4*pi*(cos(theta_out_dip) - cos(theta_in_dip)) * R_star^2;
    
            A_para = [A_para; A_para1];
            A_anti = [A_anti; A_anti1];
            A_dip = [A_dip; A_dip1];


    
    
        end
   



    end



    B_list = B_oct_area/B_dipole;
    
    figure(3)
    subplot(1,3,j)
    ax = gca;
    hold on
    plot(B_list, (A_para ./ (4*pi*R_star^2)) * 100, 'r--', LineWidth=1);    %A1_spot_para
    plot(B_list, (A_anti ./ (4*pi*R_star^2)) * 100, 'b--', LineWidth=1);    %A1_spot_anti
    grid on
    xlabel('B_{oct} / B_{dip}')
    ylabel('Accretion Filling Factor f_{acc} (%)')
    legend('Parallel Field', 'Antiparallel Field', Location='northeast')
    title(strcat('Multipole Fall Accretion Rate Scenario  ', num2str( j)))
    hold on




    Vff_para_out = ((2*G*M_star)/(R_star)).^(0.5) .* ((1./R_out_para - 1/RR1(j,2)).^0.5);
    Vff_anti_out = ((2*G*M_star)/(R_star)).^(0.5) .* ((1./R_out_anti - 1/RR1(j,2)).^0.5);
    Vff_dip_out = ((2*G*M_star)/(R_star)).^(0.5) .* ((1./R_out_dip - 1/RR1(j,2)).^0.5);

    Vff_para_in = ((2*G*M_star)/(R_star)).^(0.5) .* ((1./R_in_para - 1/RR1(j,1)).^0.5);
    Vff_anti_in = ((2*G*M_star)/(R_star)).^(0.5) .* ((1./R_in_anti - 1/RR1(j,1)).^0.5);
    Vff_dip_in = ((2*G*M_star)/(R_star)).^(0.5) .* ((1./R_in_dip - 1/RR1(j,1)).^0.5);


    
    figure(4)
    ax = gca;
    hold on
    plot(R_out_para(2568:end), Vff_para_out(2568:end)./(100.*1000),'r-', LineWidth=1);
    plot(R_out_anti, Vff_anti_out./(100.*1000),'b-', LineWidth=1);
    plot(R_out_dip, Vff_dip_out./(100.*1000),'k-', LineWidth=1);
    grid on
    ax.FontSize = 15;
    xlabel('r/R_{*}', 'FontSize', 15)
    ylabel('Vff km/s', 'FontSize', 15)
    title('Free falling velocity vs Radius', 'FontSize', 15)
    legend('Multipole Parallel', 'Multipole Antiparallel', 'Dipole')
    hold on 






    p_para_out = (B_out_para ./ B_in_para(end)) .* (M_acc_rate ./ (A_para(end).*(1000000 + Vff_para_out)));
    p_anti_out = (B_out_anti ./ B_in_anti(end)) .* (M_acc_rate ./ (A_anti(end).*(1000000 + Vff_anti_out)));
    p_dip_out = (B_out_dip ./ B_in_dip(end)) .* (M_acc_rate ./ (A_dip(end).*(1000000 + Vff_dip_out)));

%     p_para_in = (B_out_para ./ B_in_para(end)) .* (M_acc_rate ./ (A_para(end).*(1000000 + Vff_para_in)));
%     p_anti_in = (B_out_anti ./ B_in_anti(end)) .* (M_acc_rate ./ (A_anti(end).*(1000000 + Vff_anti_in)));
%     p_dip_in = (B_out_dip ./ B_in_dip(end)) .* (M_acc_rate ./ (A_dip(end).*(1000000 + Vff_dip_in)));



    n_para_out = (p_para_out / (1e-24));
    n_anti_out = (p_anti_out / (1e-24));
    n_dip_out = (p_dip_out / (1e-24));
% 
%     n_para_in = (p_para_in / (1e-24));
%     n_anti_in = (p_anti_in / (1e-24));
%     n_dip_in = (p_dip_in / (1e-24));


    figure(5)
    ax = gca;
    hold on
    plot(R_out_para(2568:end), log10(n_para_out(2568:end)),'r-', LineWidth=1);
    plot(R_out_anti, log10(n_anti_out),'b-', LineWidth=1);
    plot(R_out_dip, log10(n_dip_out),'k-', LineWidth=1);

%     plot(R_in_para, log10(n_para_in),'r--', LineWidth=1);
%     plot(R_in_anti, log10(n_anti_in),'b--', LineWidth=1);
%     plot(R_in_dip, log10(n_dip_in),'k--', LineWidth=1);
 
    grid on
    xlabel('r/R_{*}', 'FontSize', 15)
    ylabel('log(n(r))', 'FontSize', 15)
    ax.FontSize = 15;
    title('Number Density vs Radial Distance', 'FontSize', 18)
    legend('Multipole Parallel', 'Multipole Antiparallel', 'Dipole')
    hold on




end


    



figure(6)
subplot(1,3,1)
ax = gca;
hold on
plot(xnull1_out_para, znull1_out_para,'k-', LineWidth=1);
plot(xnull1_in_para, -znull1_in_para,'k-', LineWidth=1);

plot(x1_in_para, z1_in_para,'b--', LineWidth=1);
plot(x1_out_para, z1_out_para,'r--', LineWidth=1);

plot(x_star, y_star, '-', 'Color', [0.91 0.41 0.17], LineWidth=1)
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
axis equal
grid on
xlabel('X axis')
ylabel('Z axis')
title('Accretion Plotting Scenario 1')
xlim(ax, [0 5.5])
ylim(ax, [0 5.5])
hold on


subplot(1,3,2)
ax = gca;
hold on
plot(xnull2_out_para, znull2_out_para,'k-', LineWidth=1);
plot(xnull2_in_para, -znull2_in_para,'k--', LineWidth=1);

plot(x2_in_para, -z2_in_para,'b--', LineWidth=1);
plot(x2_out_para, z2_out_para,'r--', LineWidth=1);

plot(x_star, y_star, '-', 'Color', [0.91 0.41 0.17], LineWidth=1)
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
axis equal
grid on
xlabel('X axis')
ylabel('Z axis')
title('Accretion Plotting Scenario 2')
xlim(ax, [0 5.5])
ylim(ax, [0 5.5])
hold on



% figure(7)
% plot(R_out_para(1:1360), log10(n_para_out(1:1360)))


% subplot(1,3,3)
% ax = gca;
% hold on
% plot(xnull3_out_para, znull3_out_para,'k-', LineWidth=1);
% plot(xnull3_in_para, -znull3_in_para,'k-', LineWidth=1);
% 
% plot(x3_in_para, -z3_in_para,'b--', LineWidth=1);
% plot(x3_out_para, -z3_out_para,'r--', LineWidth=1);
% 
% plot(x_star, y_star, '-', 'Color', [0.91 0.41 0.17], LineWidth=1)
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% axis equal
% grid on
% xlabel('X axis')
% ylabel('Z axis')
% title('Accretion Plotting Scenario 3')
% xlim(ax, [0 5.5])
% ylim(ax, [0 5.5])
% hold on





% figure(4)
% ax = gca;
% hold on
% plot(R_para, Vff_para./(100.*1000),'r--', LineWidth=1);
% plot(R_anti, Vff_anti./(100.*1000),'b--', LineWidth=1);
% plot(R_dip_final, Vff_dip./(100.*1000),'k--', LineWidth=1);
% grid on
% xlabel('r/R_{*}')
% ylabel('log(n(r))')
% title('Number Density vs Radial Distance')
% legend('Multipole Parallel', 'Multipole Antiparallel', 'Dipole')
% hold on








%Function to plot the circle

function [x_units, y_units] = circle(x, y, r, N)
    
    angle = 0:(pi/N):(2*pi);
    x_units = x + (r * cos(angle));
    y_units = y + (r * sin(angle));

end




%This is the function to plot the multipole coordinates, I will try to
%integrate the closed lines given starting from the theta=0 line.

function [x_coords, z_coords, B_r_array, theta1, R] = multipole(r, theta, B_dipole, B_octupole, theta_dip, theta_oct, dS, G, R_star, M_star)


    x_coords = [r.*sin(theta)];
    z_coords = [r.*cos(theta)];
    B_r_array = [];
    R = [];


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
        B_r_array = [B_r_array; B];
        R = [R; r];


    
    end
    
    theta1 = theta;

end










function [y, x, r_a_dip, r_a] = multitrunc(B, B_dip, R_star, M_star, M_acc, G, R_sol)

    dt = 0.001;

    u_dip = B_dip .* ((R_star.^3)./2);

    r_a_dip = ((u_dip.^(4/7)) .* ((2.*G.*M_star).^(-1/7)) .* (M_acc.^(-2/7))) / R_sol;

 

    x = linspace(r_a_dip-1, r_a_dip+1, (1/dt)+1);

    y = x .* (1 - (3/4).*(B).*(1./(x.^2))).^(-4/7) - r_a_dip;


    r_a = interp1(y,x,0);

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







