% Simple script to test the costEval_0x1 fucntion
close all, clear all, clc

inds = rand(50,12)>0.5; % (nInds x nSws) Logical or int matrix, each row is a switch configuration
% Eliminate repeated individuals
inds = unique(inds,'rows');

% Antenna Snp file, NOTE: It is mandatory that fields are saved in the same folder as the S-parameter matrix, using the same name scheme
dataPath.SmatDir = 'D:\Users\german.ramirez\Documents\MATLAB\00 Sims and Meas\Reconfig PL\DualPort_PolDiv_3x3_PL\DualPort_PolDiv_3x3_PL.s14p';
% S2p On Sw parameters
dataPath.SwsDir_On = 'D:\Users\german.ramirez\Documents\MATLAB\00 Sims and Meas\S_Parameters_Diodes\Manufacturer\MA4AGFCP910_10mA_Forward_Bias.s2p';
 % S2p On Sw parameters
dataPath.SwsDir_Off = 'D:\Users\german.ramirez\Documents\MATLAB\00 Sims and Meas\S_Parameters_Diodes\Manufacturer\MA4AGFCP910_0V_Reverse_Bias.s2p'; 

options.Exhaustive = true;              % flag indicating if Exhaustive Search is to be performed 
options.Combination = 'none';           % (convex/none) Output creation (Non)Linear combination or MOGA
options.idealSws = false;               % Debugging option to indicate optimization with idealSws
options.connection = 'Conventional';	% Indicates the way in which loads are connected to the antenna port
options.activePorts = [1;0];%;0,1];     % Indicates the active input ports combinations, that should be served with one configuration. this can be very demanding for some structures
options.maxSval = 0;                    % Maximum value for S parameters range
options.targetS = -15;                  % Optimization target S parameters
options.maxIval = 1;                    % Maximum value for currents range (relative to fed port maximum current)
options.targetI = 0;                    % Optimization target current
options.targetTheta = 30;               % In degrees, Optimization target theta for pattern maximum
options.targetPhi = 0;                  % In degrees, Optimization target phi for pattern maximum
options.maxAngDev = 10;                 % Maximum value for angular deviation from target (In degrees)
options.f_min = 3.3e9;                  % Minimum frequency for S parameters range
options.f_max = 3.8e9;                  % Maximum frequency for S parameters range
options.patFreq = 3.55;                 % Frequency sample(s) for field combination

allIndscost = costEval_0x1(inds,dataPath,options);

figure
scatter3(allIndscost(:,1),allIndscost(:,2),allIndscost(:,3),'r*');
xlabel('Cost_{SL}'); ylabel('Cost_{\theta}'); zlabel('Cost_{Icurr}');

F = scatteredInterpolant(allIndscost(:,1),allIndscost(:,2),allIndscost(:,3),'linear','none');
F.Method = 'linear';

sgr = linspace(min(allIndscost(:,1)),max(allIndscost(:,1)));
ygr = linspace(min(allIndscost(:,2)),max(allIndscost(:,2)));
[XX,YY] = meshgrid(sgr,ygr);
ZZ = F(XX,YY);

figure
surf(XX,YY,ZZ,'LineStyle','none'); % can also be mesh 
hold on
plot3(allIndscost(:,1),allIndscost(:,2),allIndscost(:,3),'r*');
xlabel('Cost_{SL}'); ylabel('Cost_{\theta}'); zlabel('Cost_{Icurr}');

[bestMatched, posMatched] = sort(allIndscost(:,1),1);
[bestPointed, posPointed] = sort(allIndscost(:,2),1);   % This metric is prioritized
[bestCurrs, posCurrents] = sort(allIndscost(:,3),1);

if options.Exhaustive
    options.Exhaustive = false; 
    options.Combination = 'convex';
	best_Inds = dec2bin(posPointed(1:50)-1,12)-'0'; % Leftmost bit corresponds to the lowest index switch
    [a, b] = costEval_0x1(best_Inds,dataPath,options);
else
    best_Inds = inds(pos(1:50),:);
    [a, b] = costEval_0x1(inds(pos(1:50),:),dataPath,options);
end

[best, pos] = sort(a);
bestInd = best_Inds(pos(1),:);
% filteredSet = intersect(posMatched(1:50),posPointed(1:50))

magSpars_dB = 20*log10(abs(b.SL));
figure, 
for cont = 1:6
    subplot(2,3,cont);
	plot(b.freqs,magSpars_dB(pos(cont),:,1,1), b.freqs,magSpars_dB(pos(cont),:,1,2),...
    b.freqs,magSpars_dB(pos(cont),:,2,1), b.freqs,magSpars_dB(pos(cont),:,2,2) ); 
    title(['S-parameters - ' dec2bin(posPointed(pos(cont))-1,12)]);
    axis([min(b.freqs), max(b.freqs), -40,0]); grid on;
    xlabel('frequency'); ylabel('dB'); legend('S_{11}','S_{12}','S_{21}','S_{22}','location','southEast');
end

mag_SwCurr = abs(b.SwCurr);
PtSel = find(options.activePorts);
figure,
for cont = 1:2
    legCurrs{cont} = sprintf('I_{Pt %d}', cont);
end
for cont = 1:12
    legCurrs{cont+2} = sprintf('I_{Sw %d}', cont);
end
for cont = 1:6
    subplot(2,3,cont);
    plot(b.freqs,squeeze(mag_SwCurr(pos(cont),PtSel,:,:))); 
    title(['Currents Pt ', num2str(PtSel), ' - ',  dec2bin(posPointed(pos(cont))-1,12)]);
    axis([min(b.freqs), max(b.freqs), 0,0.1]); grid on;
    xlabel('frequency'); ylabel('dB'); legend(legCurrs,'location','northWest'); 
end

[aa, bb, cc] = costEval_0x1(bestInd,dataPath,options);

figure,
imagesc(squeeze(cc.Pattern(1,:,:)));
title('Combined radiation pattern for the best configuration');
xlabel('\phi'); ylabel('\theta');
set(gca,'XTick',0:10:73); set(gca,'XTickLabel',cc.phi(1,1:10:end)); 
set(gca,'YTick',0:5:37); set(gca,'YTickLabel',cc.theta(1:5:end,1)); 

MagE = 10.^(squeeze(cc.Pattern)/10);
figure,
patternCustom(MagE.', cc.theta(:,1), cc.phi(1,:));

% 10 better individuals for IDEAL switches
% port 1 fed:
    % theta = 30
        % phi = 0: 
        % phi = 90:    
        % phi = 180: 
        % phi = 270: 
% port 2 fed:
    % theta = 30
        % phi = 0:     
        % phi = 90:     
        % phi = 180: 
        % phi = 270:         

% 10 better individuals for MA4AGFCP910 switches
% port 1 fed:
    % theta = 30
        % phi = 0: 
        % phi = 90:  
        % phi = 180: 
        % phi = 270: 
% port 2 fed:
    % theta = 30
        % phi = 0: 
        % phi = 90: 
        % phi = 180: 
        % phi = 270:                