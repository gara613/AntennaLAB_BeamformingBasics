% Simple script to test the costEval_0x1 fucntion

inds = rand(50,12)>0.5; % (nInds x nSws) Logical or int matrix, each row is a switch configuration
% Eliminate repeated individuals
inds = unique(inds,'rows');

% Antenna Snp file, NOTE: It is mandatory that fields are saved in the same folder as the S-parameter matrix, using the same name scheme
dataPath.SmatDir = 'D:\Users\german.ramirez\Documents\MATLAB\00 Sims and Meas\Reconfig PL\ReconfigPL_MediumSized_pixels\ReconfigPL.s13p';
% S2p On Sw parameters
dataPath.SwsDir_On = 'D:\Users\german.ramirez\Documents\MATLAB\00 Sims and Meas\S_Parameters_Diodes\Manufacturer\MA4AGFCP910_10mA_Forward_Bias.s2p';
 % S2p On Sw parameters
dataPath.SwsDir_Off = 'D:\Users\german.ramirez\Documents\MATLAB\00 Sims and Meas\S_Parameters_Diodes\Manufacturer\MA4AGFCP910_0V_Reverse_Bias.s2p'; 

options.Exhaustive = true; %flag indicating if Exhaustive Search is to be performed 
options.Combination = 'convex';	% (convex/none) Output creation (Non)Linear combination or MOGA
options.idealSws = true;        % Debugging option to indicate optimization with idealSws
options.connection = 'DiffPts'; % Indicates the way in which loads are connected to the antenna port
options.maxSval = 0;            % Maximum value for S parameters range
options.targetS = -15;          % Optimization target S parameters
options.maxIval = 5e-3;         % Maximum value for currents range
options.targetI = 1e-1;         % Optimization target current
options.targetTheta = 30;       % In degrees, Optimization target theta for pattern maximum
options.targetPhi = 0;          % In degrees, Optimization target phi for pattern maximum
options.maxAngDev = 10;         % Maximum value for angular deviation from target (In degrees)
options.f_min = 3.45e9;         % Minimum frequency for S parameters range
options.f_max = 3.55e9;         % Maximum frequency for S parameters range
options.patFreq = 3.55;         % Frequency sample(s) for field combination

a_cost = costEval_0x1(inds,dataPath,options);

% [bestCosts, pos]= sort(a_cost);
% optimSwState = inds(pos(1),:);  % Leftmost bit corresponds to the lowest index switch
% 
% options.Exhaustive = false;     % 
% [a, b] = costEval_0x1(inds(pos(1:10),:),dataPath,options);
% 
% figure,
% plot(b.freqs,20*log10(abs(b.SL(1,:)))); % plotSpars
% title('Input port S-parameters of the best configuration')
% xlabel('frequency'); ylabel('dB'); grid on;
% 
% [aa, bb, cc] = costEval_0x1(inds(1,:),dataPath,options);
% 
% figure,
% imagesc(cc.Pattern);
% title('Combined radiation pattern for the best configuration');
% xlabel('\theta'); ylabel('\phi');
% 
% min(a_cost)
% min(a)
% b
% aa
% bb
% cc.Gain
