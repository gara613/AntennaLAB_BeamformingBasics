%% Circular antenna Array for V2I communications
% Based on the expresions in: Antenna theory - Analysis and design (Balanis 2005)
%
% Germán Augusto Ramírez
% UPC AntennaLAB

c = 3e8;
f = 5.9e9;          % Oeration frequency
lambda = c/f;
l = lambda/4;
d = lambda;         % Array diameter
a = d/2; 
k = 2*pi/lambda;

theta = (0:5:180)*pi/180;
phi = (0:5:360)*pi/180;
N_theta = length(theta);
N_phi = length(phi);

N = 16;                             % Number of elements
I = ones(N,1);                      % Feed weights (real)
phi_n = 2*pi/N*(1:N);               
theta_0 = 80*pi/180;                % for pointing toward a BS antenna at 21m from a moving car antenna at 2m
N_Beams = 8;                        % Based on the approximation HPBW = 0.88*lambda/d = 50º, choose 8 beams separated 45º
phi_0 = linspace(0,(N_Beams-1)/N_Beams*2*pi,N_Beams)';
% the phases of the elements (Used in feko simulation) are given by:
alpha_n = -k*a*sin(theta_0)*cos(bsxfun(@minus,phi_0,phi_n))*180/pi; % each row serves one phi direction

%% Array factor with maximum pointing in the (theta_0, phi_0) direction
AF_tot = zeros(N_theta,N_phi,N_Beams);
for contBeams = 1:N_Beams
    AF = zeros(N_theta,N_phi,N);
    rho_0 = a*sqrt( (sin(theta).'*cos(phi) - sin(theta_0)*cos(phi_0(contBeams))).^2 +...
        (sin(theta).'*sin(phi)-sin(theta_0)*sin(phi_0(contBeams))).^2 );
    xi = atan2( (sin(theta).'*sin(phi)-sin(theta_0)*sin(phi_0(contBeams))) ,...
        (sin(theta).'*cos(phi)-sin(theta_0)*cos(phi_0(contBeams))) );

    for cont = 1:N
        AF(:,:,cont) = I(cont).*exp(1i*k*rho_0.*cos(phi_n(cont)-xi));
    end
    AF_tot(:,:,contBeams) = sum(AF,3);
end

%% some plots
beamNumber = [1;4;6];
indTheta = find(theta == theta_0); % may not be exact, look for the closest
[~,indPhi] = intersect(phi,phi_0(beamNumber));

figure, 
for contBeam = 1:length(beamNumber)
    polar(phi,(abs(AF_tot(indTheta,:,beamNumber(contBeam ))))); 
    hold on; 
    polar([phi_0(beamNumber(contBeam)) phi_0(beamNumber(contBeam))],[0 max(abs(AF_tot(:,indPhi(contBeam),beamNumber(contBeam))))], 'r');
end    
title('AF(\phi) for \theta = \theta_0'); 

figure,  
for contBeam = 1:length(beamNumber)
    polar(theta,(abs(AF_tot(:,indPhi(contBeam),beamNumber(contBeam)))).'); 
    hold on;
	polar([theta_0 theta_0],[0 max(abs(AF_tot(:,indPhi(contBeam),beamNumber(contBeam)))).'], 'r');
end
title('AF(\theta) for \phi = \phi_0'); 

figure, imagesc(phi*180/pi,theta*180/pi,abs(AF_tot(:,:,4))); xlabel('\phi'), ylabel('\theta'); colorbar;
% [phi_G,theta_G] = meshgrid(phi,theta);
% [x_G,y_G,AF_Z] = sph2cart(phi_G,theta_G,abs(AF_tot(:,:,1)));
% figure, surf(x_G,y_G,AF_Z);