%
%   AntiH.m
%
clear all;
close all;

npts=2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   User defined parameters in between the lines

N_L = 3;                        %Number of coils vertically
N_r = 3;                        %Number of coils radially
I = 50;                        %current in amps
t = 0.067 * 2.54;            %diameter of coated wire - 16 gauge wire
r_ave = 2.25 * 2.54;            %interior radius in cm
L = r_ave;          %distance of center coils for Anti-Helmholtz
                                %configuration. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = 2*[-npts:npts]*L/(2*npts);  %Plotting beyond the location of the coils
B=zeros(1,length(z));

for j=1:N_r                     %radial loop
    r_j = r_ave + (j-(N_r+1)/2)*t;
    A = 2 * pi * I * r_j^2 / 10;
    for i=1:N_L                 % Longitudinal loop
        z_i = (L/2) +(i-(N_L+1)/2) * t ;
        C_plus = ((z - z_i).^2 + r_j^2).^(3/2);
        C_minus = ((z + z_i).^2 + r_j^2).^(3/2);
        B = B + A * (1./C_plus - 1./C_minus);           
    end
end

fid1=figure;
plot(z,B);
set(fid1,'Position',[41 578 560 420]);
xlabel('z- cm'); ylabel('B(z)- G');
title('B');

fid2=figure;
plot(z(1:length(z)-1),diff(B)./diff(z))
set(fid2,'Position',[43 84 560 420]);
xlabel('z- cm'); ylabel('dB/dz- G/cm');
title('dB');

fid3=figure;
plot(z(1:length(z)-2),diff(diff(B))./diff(z(1:length(z)-1)).^2)

% fid4=figure;
% plot(z(1:length(z)-3),diff(diff(diff(B)))./diff(z(1:length(z)-2)).^2)
% title('Third derivative')


coil_area = pi * (t/2)^2;
N_t = 2*(N_L * N_r);
pwr = 2 * pi * r_ave * (1.7e-6) * (N_t * I)^2 / (N_t * coil_area);
volts = pwr / I;
display(pwr);
display(volts);

mu = (4*pi)*10e-7;
inductance = N_t^2 * (r_ave/100) * mu * (log(8*r_ave/coil_area) - 2);
display(inductance);

clamp_V = 500;
current = I;
switch_L = inductance;
switching_T = switch_L * current / clamp_V;
disp(strcat('Switching time = ',num2str(switching_T/1e-6),' usec'));
