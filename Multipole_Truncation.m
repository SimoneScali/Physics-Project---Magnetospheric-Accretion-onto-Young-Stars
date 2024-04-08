%This documents will plot examples graphs for the multipole components of
%the magnetic field so that i will be able to compare them with the dipole
%ones

%First I have to calculate the truncation radius for a multipole component

%Firstly, I will define some constants in cgs units

R_sol = 6.957 * 10^10; 
M_sol = 1.989 * 10^33;
G = 6.674 * 10^-8;
year = 3.1536 * 10^7;
day = 86400;

beta = 0.6;  %Parameter


%Then, I will define values for tests stars

dt = 0.001;      %steps for accuracy

%B_dip = linspace(1000, 2000, (1/dt)+1);   %values in gauss
%B_oct = linspace(3000, 6000, (1/dt)+1);

M_star = M_sol;      %linspace(M_sol, 2*M_sol, (1/dt)+1);
R_star = 2 * R_sol;       %linspace(1*R_sol, 2*R_sol, (1/dt)+1);
Period = 7 * day;           %linspace(5*day, 10*day, (1/dt)+1);




%Now I Will calculate the truncation radius for the multipole component for
%a single test star

B_oct = 5000;
B_dip = 1000;
B = B_oct./B_dip;
R_star = 2 * R_sol;
M_star = M_sol;
M_acc = ((10e-8) .* M_sol) / year;



[y, x, r_a_dip] = multitrunc(B_oct, B_dip, R_star, M_star, M_acc, G, R_sol);

yy = zeros(1, (1/dt)+1);
xx = x;

r_a = interp1(y,x,0);

figure(1)
plot(x, y, 'b--', xx, yy, 'k-', r_a, 0, 'r+', LineWidth=1.5)
legend('Root Finding', 'y = 0', strcat('r_a =  ', num2str(r_a)), Location='northwest');
title('Findnig R_T for the multipole')
xlabel('r_a / R_*')
ylabel('f(x)')
grid on

hold off

R_T = beta .* r_a;


%Now I will show how the truncation radius changes if the mass accretion
%rate is kept constant but the magnetic field changes

B1 = linspace(2, 6, (1/dt)+1);   %This is B_oct / B_dip

p = 8:0.2:9;
M_acc_const = ((10.^(-p)).*M_sol) ./ year;

beta = 0.6;

B_dip = linspace(250, 2000, (1/dt)+1);   %values in gauss

R_t_AC = zeros(length(p), (1/dt)+1);
R_t_AC_DIP = zeros(length(p), (1/dt)+1);


BDIP = linspace(500, 2000, length(p));   %This is B_dip
BOCT = linspace(0, 10000, (1/dt)+1);   %This is B_oct



for i = 1:length(p)

    for j = 1:((1/dt) + 1)

        [y, x, r_a_dip] = multitrunc(BOCT(j), BDIP(i), R_star, M_star, M_acc_const(i), G, R_sol);

        if 0 == isreal(y);
            y = abs(y);

        end

        r_a = interp1(y,x,0);

        R_t_AC(i, j) = beta .* r_a;

        R_t_AC_DIP(i, j) = beta .* r_a_dip;

    end

    
end




figure(2)
subplot(2,1,1) 
ax = gca;
plot(R_t_AC(1, 1:end), BOCT./BDIP(1), 'b-', R_t_AC(2, 1:end), BOCT./BDIP(2), 'b-', R_t_AC(3, 1:end), BOCT./BDIP(3), 'b-', R_t_AC(4, 1:end), BOCT./BDIP(4), 'b-', R_t_AC(5, 1:end), BOCT./BDIP(5), 'b-', R_t_AC(6, 1:end), BOCT./BDIP(6), 'b-', LineWidth=1.5)
hold on
plot(R_t_AC_DIP(1, 1:end), BOCT./BDIP(1), 'r--', R_t_AC_DIP(2, 1:end), BOCT./BDIP(2), 'r--', R_t_AC_DIP(3, 1:end), BOCT./BDIP(3), 'r--', R_t_AC_DIP(4, 1:end), BOCT./BDIP(4), 'r--', R_t_AC_DIP(5, 1:end), BOCT./BDIP(5), 'r--', R_t_AC_DIP(6, 1:end), BOCT./BDIP(6), 'r--', LineWidth=1.5)
axis equal
ax.FontSize = 15;
legend({'Multipole magnetic field', '','','','','', 'Dipole magnetic field'}, Location='southeast')
% title('Disc truncation variance - Constant accretion', 'FontSize', 14)
xlabel('R_T  [R_*]', 'FontSize', 15)
ylabel('B_{oct}/B_{dip}','FontSize', 15)
xlim(ax, [3 16])
ylim(ax, [0 5])
grid on

hold off


%Now I will try to plot the graph of the mass accretion rate against the
%truncation radius keeping the magnetic multipole field constant for a few
%different values


B2 = linspace(2000, 6000, length(p));

M_acc = linspace(10^-9, 10^-8, (1/dt)+1) .* (M_sol ./ year);
B_dip_const = linspace(1000, 2000, length(p));

R_t_BC = zeros(length(p), (1/dt)+1);
R_t_BC_DIP = zeros(length(p), (1/dt)+1);


for i = 1:length(p)

    for j = 1:((1/dt) + 1)

        [y, x, r_a_dip] = multitrunc(B2(i), B_dip_const(i), R_star, M_star, M_acc(j), G, R_sol);

        r_a = interp1(y,x,0);

        R_t_BC(i, j) = beta .* r_a;

        R_t_BC_DIP(i, j) = beta .* r_a_dip;


    end

end




subplot(2,1,2)
ax = gca;
plot(R_t_BC(1, 1:end), (M_acc .* (year/M_sol)), 'b-', R_t_BC(2, 1:end), (M_acc .* (year/M_sol)), 'b-', R_t_BC(3, 1:end), (M_acc .* (year/M_sol)), 'b-', R_t_BC(4, 1:end), (M_acc .* (year/M_sol)), 'b-', R_t_BC(5, 1:end), (M_acc .* (year/M_sol)), 'b-', R_t_BC(6, 1:end), (M_acc .* (year/M_sol)), 'b-', LineWidth=1.5)
hold on
plot(R_t_BC_DIP(1, 1:end), (M_acc .* (year/M_sol)), 'r--', R_t_BC_DIP(2, 1:end), (M_acc .* (year/M_sol)), 'r--', R_t_BC_DIP(3, 1:end), (M_acc .* (year/M_sol)), 'r--', R_t_BC_DIP(4, 1:end), (M_acc .* (year/M_sol)), 'r--', R_t_BC_DIP(5, 1:end), (M_acc .* (year/M_sol)), 'r--', R_t_BC_DIP(6, 1:end), (M_acc .* (year/M_sol)), 'r--', LineWidth=1.5)
ax.FontSize = 15;
legend({'Multipole magnetic field', '','','','','', 'Dipole magnetic field'}, Location='northeast')
% title('Disc truncation variance - Constant magnetic field magnitude', 'FontSize', 14)
xlabel('R_T  [R_*]', 'FontSize', 15)
ylabel('Mass accretion rate [M_☉ / year]', 'FontSize', 12)
% axis equal
xlim(ax, [5 15])
ylim(ax, [0 10e-9])
grid on

hold off



%Now I will calculate the corotation radius and produce the same two graphs
%above, but the x-axis will be in terms of R_t / R_co

R_co = (((G .* M_star) ./ (((2*pi) ./ Period).^2)) .^ (1/3)) ./ R_sol;

R_t_AC_co = R_t_AC ./ R_co;
R_t_BC_co = R_t_BC ./ R_co;

R_t_AC_co_DIP = R_t_AC_DIP ./ R_co;
R_t_BC_co_DIP = R_t_BC_DIP ./ R_co;

% subplot(2,2,3)
% ax = gca;
% plot(R_t_AC_co(1, 1:end), B1, 'b--', R_t_AC_co(2, 1:end), B1, 'r--', R_t_AC_co(3, 1:end), B1, 'k--', R_t_AC_co(4, 1:end), B1, 'g--', R_t_AC_co(5, 1:end), B1, 'y--', R_t_AC_co(6, 1:end), B1, 'm--', LineWidth=1.5)
% hold on
% plot(R_t_AC_co_DIP(1, 1:end), B1, 'b:', R_t_AC_co_DIP(2, 1:end), B1, 'r:', R_t_AC_co_DIP(3, 1:end), B1, 'k:', R_t_AC_co_DIP(4, 1:end), B1, 'g:', R_t_AC_co_DIP(5, 1:end), B1, 'y:', R_t_AC_co_DIP(6, 1:end), B1, 'm:', LineWidth=1.5)
% %legend('Accretion 1', 'Accretion 2', 'Accretion 3', 'Accretion 4', 'Accretion 5', 'Accretion 6', Location='northwest')
% title('Multipole magnetic fields trial 2 - Constant accretion')
% xlabel('R_t / R_{co}  [unitless]')
% ylabel('B_{oct}/B_{dip} [unitless]')
% %xlim(ax, [5 15])
% %ylim(ax, [-0.1 0.5])
% grid on
% 
% hold off
% 
% 
% 
% subplot(2,2,4)
% ax = gca;
% plot(R_t_BC_co(1, 1:end), (M_acc .* (year/M_sol)), 'b--', R_t_BC_co(2, 1:end), (M_acc .* (year/M_sol)), 'r--', R_t_BC_co(3, 1:end), (M_acc .* (year/M_sol)), 'k--', R_t_BC_co(4, 1:end), (M_acc .* (year/M_sol)), 'g--', R_t_BC_co(5, 1:end), (M_acc .* (year/M_sol)), 'y--', R_t_BC_co(6, 1:end), (M_acc .* (year/M_sol)), 'm--', LineWidth=1.5)
% hold on
% plot(R_t_BC_co_DIP(1, 1:end), (M_acc .* (year/M_sol)), 'b:', R_t_BC_co_DIP(2, 1:end), (M_acc .* (year/M_sol)), 'r:', R_t_BC_co_DIP(3, 1:end), (M_acc .* (year/M_sol)), 'k:', R_t_BC_co_DIP(4, 1:end), (M_acc .* (year/M_sol)), 'g:', R_t_BC_co_DIP(5, 1:end), (M_acc .* (year/M_sol)), 'y:', R_t_BC_co_DIP(6, 1:end), (M_acc .* (year/M_sol)), 'm:', LineWidth=1.5)
% %legend('B 1', 'B 2', 'B 3', 'B 4', 'B 5', 'B 6', Location='northwest')
% title('Multipole magnetic fields trial 2 - Constant Magnetic Field')
% xlabel('R_t / R_{co}  [unitless]')
% ylabel('Accretion Mass rate [M_o / year]')
% %xlim(ax, [5 15])
% %ylim(ax, [-0.1 0.5])
% grid on
% 
% hold off


% B = 6;
% B_dip = 500;
% M_acc = (1e-8 .* M_sol) / year;
% M_star = 0.7 * M_sol;
% 
% 
% [y, x, r_a_dip] = multitrunc(B, B_dip, R_star, M_star, M_acc, G, R_sol);
% 
% r_a = interp1(y,x,0);
% 
% R_T = r_a * 0.6;


function [y, x, r_a_dip] = multitrunc(B_oct, B_dip, R_star, M_star, M_acc, G, R_sol)

    dt = 0.001;

    B = B_oct ./ B_dip;

    u_dip = B_dip .* ((R_star.^3)./2);

    r_a_dip = ((u_dip.^(4/7)) .* ((2.*G.*M_star).^(-1/7)) .* (M_acc.^(-2/7))) / R_sol;

    

    x = linspace(r_a_dip-1, r_a_dip+1, (1/dt)+1);



    y = x .* (1 - (3/4).*(B).*(1./(x.^2))).^(-4/7) - r_a_dip;


end





