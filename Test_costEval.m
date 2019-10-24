% Simple script to test the costEval_0x1 fucntion
close all, clear all, clc

inds = rand(50,12)>0.5; % (nInds x nSws) Logical or int matrix, each row is a switch configuration
% Eliminate repeated individuals
inds = unique(inds,'rows');

% Antenna Snp file, NOTE: It is mandatory that fields are saved in the same folder as the S-parameter matrix, using the same name scheme
dataPath.SmatDir = 'D:\Users\german.ramirez\Documents\MATLAB\00 Sims and Meas\Reconfig PL\\ReconfigPL\ReconfigPL.s13p';
% S2p On Sw parameters
dataPath.SwsDir_On = 'D:\Users\german.ramirez\Documents\MATLAB\00 Sims and Meas\S_Parameters_Diodes\Manufacturer\MA4AGFCP910_10mA_Forward_Bias.s2p';
 % S2p On Sw parameters
dataPath.SwsDir_Off = 'D:\Users\german.ramirez\Documents\MATLAB\00 Sims and Meas\S_Parameters_Diodes\Manufacturer\MA4AGFCP910_0V_Reverse_Bias.s2p'; 

options.Exhaustive = true;      %flag indicating if Exhaustive Search is to be performed 
options.Combination = 'convex';	% (convex/none) Output creation (Non)Linear combination or MOGA
options.idealSws = true;        % Debugging option to indicate optimization with idealSws
options.connection = 'Conventional'; % Indicates the way in which loads are connected to the antenna port
options.activePorts = 1;        %
options.maxSval = 0;            % Maximum value for S parameters range
options.targetS = -15;          % Optimization target S parameters
options.maxIval = 1;            % Maximum value for currents range (relative to fed port maximum current)
options.targetI = 0;            % Optimization target current
options.targetTheta = 40;       % In degrees, Optimization target theta for pattern maximum
options.targetPhi = 0;          % In degrees, Optimization target phi for pattern maximum
options.maxAngDev = 10;         % Maximum value for angular deviation from target (In degrees)
options.f_min = 3.45e9;         % Minimum frequency for S parameters range
options.f_max = 3.55e9;         % Maximum frequency for S parameters range
options.patFreq = 3.55;         % Frequency sample(s) for field combination

a_cost = costEval_0x1(inds,dataPath,options);

[bestCosts, pos]= sort(a_cost);
% After exhaustive search, best states are: [3828,3012,2884,3908,4036,2564,3588,3716,2692,3860]

if options.Exhaustive
    options.Exhaustive = false; 
    bestInd = dec2bin(pos(1)-1,12)-'0'; % Leftmost bit corresponds to the lowest index switch
    [a, b] = costEval_0x1(dec2bin(pos(1:10)-1,12)-'0',dataPath,options);
else
    bestInd = inds(pos(1),:);
    [a, b] = costEval_0x1(inds(pos(1:10),:),dataPath,options);
end
figure,
plot(b.freqs,20*log10(abs(b.SL(1,:)))); % plotSpars
title('Input port S-parameters of the best configuration')
xlabel('frequency'); ylabel('dB'); grid on;

bestInd
[aa, bb, cc] = costEval_0x1(bestInd,dataPath,options);

figure,
plot(bb.freqs,20*log10(abs(bb.SL))); % plotSpars
title('Input port S-parameters of the best configuration')
xlabel('frequency'); ylabel('dB'); grid on;

figure,
imagesc(squeeze(cc.Pattern(1,:,:)));
title('Combined radiation pattern for the best configuration');
xlabel('\phi'); ylabel('\theta');
set(gca,'XTick',0:10:73); set(gca,'XTickLabel',cc.phi(1,1:10:end)); 
set(gca,'YTick',0:5:37); set(gca,'YTickLabel',cc.theta(1:5:end,1)); 

min(a_cost)
min(a)
b
aa
bb
