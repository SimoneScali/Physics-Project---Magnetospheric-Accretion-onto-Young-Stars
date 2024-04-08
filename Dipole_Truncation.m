%This documents will plot examples graphs for the first tries of the
%project

%The first graph I will want to produce is the magnetic field of a dipole
%against the truncation radius

%Firstly, I will define some constants in cgs units

R_sol = 6.957 * 10^10; 
M_sol = 1.989 * 10^33;
G = 6.674 * 10^-8;
year = 3.1536 * 10^7;
day = 86400;

beta = 0.6;  %Parameter


%Then, I will define values for tests stars

dt = 0.001;      %steps for accuracy

B_dip = linspace(250, 2000, (1/dt)+1);   %values in gauss

M_star = M_sol;      %linspace(M_sol, 2*M_sol, (1/dt)+1);
R_star = 2 * R_sol;       %linspace(1*R_sol, 2*R_sol, (1/dt)+1);
Period = 7 * day;           %linspace(5*day, 10*day, (1/dt)+1);

p = 8:0.2:9;



M_acc_const = ((10.^(-p)).*M_sol) ./ year;

u_dip = B_dip .* ((R_star.^3)./2);

%Now I will calculate the truncation radius in R_sol keeping the accretion rate
%constant for a few different values


R_t_AC = zeros(length(p), (1/dt)+1);

for i = 1:length(p)

    for j = 1:((1/dt) + 1)

        R_t_AC(i, j) = (beta .* (u_dip(j).^(4/7)) .* ((2.*G.*M_star).^(-1/7)) .* (M_acc_const(i).^(-2/7))) / R_sol;

    end
end

subplot(2,2,1)                                   
plot(R_t_AC(1, 1:end), B_dip, 'b--', R_t_AC(2, 1:end), B_dip, 'r--', R_t_AC(3, 1:end), B_dip, 'k--', R_t_AC(4, 1:end), B_dip, 'g--', R_t_AC(5, 1:end), B_dip, 'y--', R_t_AC(6, 1:end), B_dip, 'm--', LineWidth=1.5)
legend('Accretion 1', 'Accretion 2', 'Accretion 3', 'Accretion 4', 'Accretion 5', 'Accretion 6', Location='northwest')
title('Dipole magnetic fields trial 1 - Constant accretion')
xlabel('R_t  [R_{sol}]')
ylabel('B_{dip} [Gauss]')
grid on

hold off


%Now I will try to plot the graph of the mass accretion rate against the
%truncation radius keeping the magnetic dipole field constant for a few
%different values


B_dip_const = linspace(1000, 2000, length(p));

M_acc = linspace(10^-9, 10^-8, (1/dt)+1) .* (M_sol ./ year);


u_dip_BC = zeros(length(p), (1/dt)+1);
R_t_BC = zeros(length(p), (1/dt)+1);

for i = 1:length(p)

    for j = 1:((1/dt) + 1)

        u_dip_BC(i, j) = B_dip_const(i) .* ((R_star.^3)./2);

        R_t_BC(i, j) = (beta .* (u_dip_BC(i, j).^(4/7)) .* ((2.*G.*M_star).^(-1/7)) .* (M_acc(j).^(-2/7))) / R_sol;

    end
end


subplot(2,2,2)                                     
plot(R_t_BC(1, 1:end), (M_acc .* (year/M_sol)), 'b--', R_t_BC(2, 1:end), (M_acc .* (year/M_sol)), 'r--', R_t_BC(3, 1:end), (M_acc .* (year/M_sol)), 'k--', R_t_BC(4, 1:end), (M_acc .* (year/M_sol)), 'g--', R_t_BC(5, 1:end), (M_acc .* (year/M_sol)), 'y--', R_t_BC(6, 1:end), (M_acc .* (year/M_sol)), 'm--', LineWidth=1.5)
legend('B 1', 'B 2', 'B 3', 'B 4', 'B 5', 'B 6', Location='northwest')
title('Dipole magnetic fields trial 1 - Constant Magnetic Field')
xlabel('R_t  [R_{sol}]')
ylabel('Accretion Mass rate [M_o / year]')
grid on

hold off



%Now I will calculate the corotation radius and produce the same two graphs
%above, but the x-axis will be in terms of R_t / R_co

R_co = (((G .* M_star) ./ (((2*pi) ./ Period).^2)) .^ (1/3)) ./ R_sol;

R_t_AC_co = R_t_AC ./ R_co;
R_t_BC_co = R_t_BC ./ R_co;

subplot(2,2,3)                                    
plot(R_t_AC_co(1, 1:end), B_dip, 'b--', R_t_AC_co(2, 1:end), B_dip, 'r--', R_t_AC_co(3, 1:end), B_dip, 'k--', R_t_AC_co(4, 1:end), B_dip, 'g--', R_t_AC_co(5, 1:end), B_dip, 'y--', R_t_AC_co(6, 1:end), B_dip, 'm--', LineWidth=1.5)
legend('Accretion 1', 'Accretion 2', 'Accretion 3', 'Accretion 4', 'Accretion 5', 'Accretion 6', Location='northwest')
title('Dipole magnetic fields trial 2 - Constant accretion')
xlabel('R_t / R_{co}  [unitless]')
ylabel('B_{dip} [Gauss]')
grid on

hold off



subplot(2,2,4)                                    
plot(R_t_BC_co(1, 1:end), (M_acc .* (year/M_sol)), 'b--', R_t_BC_co(2, 1:end), (M_acc .* (year/M_sol)), 'r--', R_t_BC_co(3, 1:end), (M_acc .* (year/M_sol)), 'k--', R_t_BC_co(4, 1:end), (M_acc .* (year/M_sol)), 'g--', R_t_BC_co(5, 1:end), (M_acc .* (year/M_sol)), 'y--', R_t_BC_co(6, 1:end), (M_acc .* (year/M_sol)), 'm--', LineWidth=1.5)
legend('B 1', 'B 2', 'B 3', 'B 4', 'B 5', 'B 6', Location='northwest')
title('Dipole magnetic fields trial 2 - Constant Magnetic Field')
xlabel('R_t / R_{co}  [unitless]')
ylabel('Accretion Mass rate [M_o / year]')
grid on

hold off











