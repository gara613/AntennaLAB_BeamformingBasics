N = 1:35;
kd = pi;    % half wavelength inter-element separation
alpha = 0;	% broadside radiation
D = zeros(35,1);

for cont = N
    In = ones(1,cont);
    D(cont) = arrayDir(In,kd,alpha);  
end
dir_dB = 10*log10(D);
hpbw = 2./D*180/pi;
hpbw_appr = 0.88*2./N*180/pi;

% Directivity of an antenna array vs number of antennas
figure, 
plot(N,dir_dB);
title('Directivity of an antenna array');
xlabel('Number of antennas'); ylabel('Directivity (dB)'); grid on

% Half-power beamwidth of an antenna array vs number of antennas
figure, 
plot(N,hpbw,N,hpbw_appr)
title('Half power Beamwidth of an antenna array');
xlabel('Number of antennas'); ylabel('HPBW (º)'); grid on

% Required gain to get the same coverage as a monopole at 5.9 GHz
f_base = 5.9e9;
lambda_base = 3e8/f_base;
G_ref = 10^(5.19/10);
P_ref = G_ref*(lambda_base/(4*pi))^2;
freq = 1e9:1e9:100e9;
lambda = 3e8./freq;
Pr = (lambda/(4*pi)).^2;

G_dif = 10*log10(P_ref./Pr);
Grx = 10*log10(G_ref*(lambda_base./lambda).^2);

figure, 
plot(freq,G_dif,freq,Grx);
title('Required gain to get the same coverage as a monopole at 5.9 GHz');
xlabel('frequency (Hz)'); ylabel('Gain (dB)'); grid on
