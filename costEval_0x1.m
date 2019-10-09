% function costEval_1x0(inds,dataPath,options)
% Cost evaluator for the optimization of switch states of PL and RA using GA
%
% Inputs: (No default values admitted)
%	* inds -> (nInds x nSws) Logical or int matrix, each row is a switch configuration
%   * dataPath -> Struct with fields:
%       - SmatDir       Antenna Snp file
%       - SwsDir_On     S2p On Sw parameters
%    	- SwsDir_Off    S2p On Sw parameters
%           + NOTE: It is mandatory that fields are saved in the same folder as the S
%             parameter matrix, using the same name shceme: e.g 'ReconfigPL.s13p' ->
%             'ReconfigPL_f_2.45_Pt_n.fss', where 'n' is the excieted port ('2.45' is just an example)
%   * options -> Struct with fields:
%       - Exhaustive -> (logic) indicating if Exhaustive Search is to be performed, WARNING: ignores provided inds
%       - Combination -> (String: Convex/none) Output creation (Non)Linear combination (provide weights and create a more interesting NL mapping) or vector output for MOGA 
%       - idealSws -> (logical) Debugging option to indicate optimization with idealSws
%       - connection -> (string: DiffPts/Conventional) indicates the way in which loads are connected to the antenna port      
%       - maxSval -> (float) Maximum value for S parameters range
%       - targetS -> (float) Optimization target S parameters
%       - maxIval -> (float) Maximum value for currents range
%       - targetI -> (float) Optimization target current
%       - targetTheta -> (float) In degrees, Optimization target theta for pattern maximum
%       - targetPhi -> (float) In degrees, Optimization target phi for pattern maximum
%       - maxAngDev -> (float) In degrees, Maximum value for angular deviation from target
%    	- f_min -> (float) Minimum frequency for S parameters range
%       - f_max -> (float) Maximum frequency for S parameters range
%       - patFreq -> (float) Frequency sample(s) for field combination
%   NOTE: Future release should read options from file rather than from a user-created structure
% Outputs: 
%   * Cost vector
% 
% Germán A. Ramírez, Sept 2019 

function varargout = costEval_0x1(inds,dataPath,options,varargin)
    %% Initialize variables
    [nInds,nSws,dataroute,SW_S_pars_path_On,SW_S_pars_path_Off,...
    ExhaustiveSearch,Combination,idealSws,connection,maxSval,targetS,targetI,maxIval,...
    targetTheta,targetPhi,maxAngDev,f_min,f_max,patFreq] = inputsRead(inds,dataPath,options,varargin);

    if ExhaustiveSearch
        disp('Exhaustive search requested, this can take a long time');
        inds = dec2bin(0:2^nSws-1,nSws)-'0'; %this overrides the individuals provided
        nInds = size(inds,1);
    end

    Cost_S_L = ones(nInds,1);
    Cost_Sw_I = ones(nInds,1);
    Cost_AngDev = ones(nInds,1);
    
    %% Read data - Include a validations to ensure reliable and consistent data is loaded
    [filepath,name,ext] = fileparts(dataroute);
    [Smat, freqsSim, Z_ref] = readSpars(dataroute); % Total S-matrix (N_ports x N_ports x N_freqs) from CST simulation
    totN_ports = size(Smat,1);
    N_A = totN_ports - nSws;    % N_A: Number of antenna feed ports, 
    N_B = nSws;                 % N_B: Number of impedance loaded ports

    fieldFileName = fullfile(filepath,[name, '_f_', regexprep(num2str(patFreq),'\.','_'), '_Pt_']);
    for cont=1:totN_ports 
        [rad_E_field(cont), thetaGrid, phiGrid, P_ant, freqPat] = readCST_FarFieldSrc([fieldFileName num2str(cont) '.ffs']);
    end

    if idealSws 
        freqs_Sw = freqsSim(1:2:end); % Decimation to reduce exhaustive search overhead
        interpSparsSw = 'none';
        Spars_Sw_on = [-ones(1,1,length(freqsSim)) zeros(1,1,length(freqsSim)); zeros(1,1,length(freqsSim)) -ones(1,1,length(freqsSim))];
        Spars_Sw_off = [ones(1,1,length(freqsSim)) zeros(1,1,length(freqsSim)); zeros(1,1,length(freqsSim)) ones(1,1,length(freqsSim))];
    else	% Load the total S matrix (2 x 2 x N_freqs) of the intended switch for its ON/OFF state
        interpSparsSw = 'swRange'; 
        [Spars_Sw_on, freqs_Sw, Zr_Sw]  = readSpars(SW_S_pars_path_On);  
        [Spars_Sw_off, ~, ~]            = readSpars(SW_S_pars_path_Off);
    end

    %% Interpolate Sw S-parameters - Needs revision: take as few sample points as possible to avoid computation burden
    if strcmp(interpSparsSw, 'swRange') % Interpolates over the Switch frequencies
        interpFact = 100;
        Nsamps = interpFact*(length(freqs_Sw)-1)+1;
        for cont1 = 1:2
            for cont2 = 1:2
                [~, ~, ~, S_Sw_on(cont1,cont2,:)] = rationalSpar(freqs_Sw, squeeze(Spars_Sw_on(cont1,cont2,:)), 'canonical',20, Nsamps);
                [~, ~, ~, S_Sw_off(cont1,cont2,:)] = rationalSpar(freqs_Sw, squeeze(Spars_Sw_off(cont1,cont2,:)), 'canonical',20, Nsamps);
            end
        end
        freqs_Sw = linspace(min(freqs_Sw),max(freqs_Sw),Nsamps)'; 
    elseif strcmp(interpSparsSw, 'simRange') % Interpolates in the simulation frequency range (assumed lower than Sw freq range)
        Nsamps = size(Smat,3);
        [freqsComm, ~, indsB] = intersect(freqsSim, freqs_Sw); 
        for cont1 = 1:2
            for cont2 = 1:2
                [~, ~, ~, S_Sw_on(cont1,cont2,:)] = rationalSpar(freqsComm, squeeze(Spars_Sw_on(cont1,cont2,indsB)), 'canonical',20, Nsamps);
                [~, ~, ~, S_Sw_off(cont1,cont2,:)] = rationalSpar(freqsComm, squeeze(Spars_Sw_off(cont1,cont2,indsB)), 'canonical',20, Nsamps); 
            end
        end            
        freqs_Sw = linspace(min(freqs_Sw(indsB)),max(freqs_Sw(indsB)),Nsamps);
    else
        S_Sw_on = Spars_Sw_on;
        S_Sw_off = Spars_Sw_off;
    end

    %% Retain only the common frequencies to perform the required calculations
    [freqsComm, indsA, indsB] = intersect(freqsSim, freqs_Sw);     %freqsComm = freqsSim(indsA) = freqs_Sw(indsB)
    nfreqs = length(freqsComm);

% 	f_range_S = freqsComm(freqsComm>f_min & freqsComm<f_max);
%   f_indx_S = (freqsComm>f_min & freqsComm<f_max);
    
    if idealSws 
        Y_Sw_on = 1e15*ones(nfreqs,1);
        Y_Sw_off = zeros(nfreqs,1);
    else
        Y_Sw_on = myS2Y(S_Sw_on(:,:,indsB), Z_ref);   
        Y_Sw_off = myS2Y(S_Sw_off(:,:,indsB), Z_ref);
    end

    %% Initialize variables for processing 
    % S_L = cell(1,nInds);        % Each cell will contain a N_A x N_A x nFreqs array
    S_L = zeros(nInds,nfreqs,N_A,N_A);
    AngDev = zeros(nInds,2); 
    Icurr = zeros(nInds,totN_ports,nfreqs);
    
    Idmat = eye(N_B);

    % S parameters matrix at the common frequencies:
    Smat = Smat(:,:, indsA);
    Y = zeros(totN_ports,totN_ports,nfreqs);
    for cont = 1:nfreqs
        Y(:,:,cont) = 1/Z_ref* inv(Smat(:,:,cont) + eye(totN_ports))* (eye(totN_ports) - Smat(:,:,cont));     
    end
    % S parameters matrix segmentation (This assumes ascending port numbering, where first N_A = feed ports):
    S_AA = Smat(1:N_A, 1:N_A, :);
    S_AB = Smat(1:N_A, N_A+1:N_A+N_B, :);
    S_BA = Smat(N_A+1:N_A+N_B, 1:N_A, :);
    S_BB = Smat(N_A+1:N_A+N_B, N_A+1:N_A+N_B, :);
    % Y parameters matrix segmentation:
    Y_AA = Y(1:N_A, 1:N_A, :);
    Y_AB = Y(1:N_A, N_A+1:N_A+N_B, :);
    Y_BA = Y(N_A+1:N_A+N_B, 1:N_A, :);
    Y_BB = Y(N_A+1:N_A+N_B, N_A+1:N_A+N_B, :);

    %% Evaluate individuals cost according to switch state, take care when using massive amounts of swithces nSws>12!
    for contInd = 1:nInds
        testInd = inds(contInd,:) > 0.5; 
        Gamma_BB = zeros(N_B,N_B,nfreqs);  % reflection coefficients matrix, depends on the state of the switches
        Ypt = zeros(N_B,N_B,nfreqs);
        V_pt = 1;

        for cont = 1:N_B
            if idealSws   
                if testInd(cont) == true
                    S11_Sw = squeeze(S_Sw_on(1,1,indsB));	
                    Ypt(cont,cont,:) = Y_Sw_on;
                elseif testInd(cont) == false
                    S11_Sw = squeeze(S_Sw_off(1,1,indsB));	
                    Ypt(cont,cont,:) = Y_Sw_off;
                end
                Gamma_BB(cont,cont,:) = S11_Sw;
            else 
                if testInd(cont) == true
                    S11_Sw = squeeze(S_Sw_on(1,1,indsB));
                    S12_Sw = squeeze(S_Sw_on(1,2,indsB));
                    S21_Sw = squeeze(S_Sw_on(2,1,indsB));
                    S22_Sw = squeeze(S_Sw_on(2,2,indsB));
                    Y_Sw = Y_Sw_on;
                elseif testInd(cont) == false
                    S11_Sw = squeeze(S_Sw_on(1,1,indsB));
                    S12_Sw = squeeze(S_Sw_on(1,2,indsB));
                    S21_Sw = squeeze(S_Sw_on(2,1,indsB));
                    S22_Sw = squeeze(S_Sw_on(2,2,indsB));
                    Y_Sw = Y_Sw_off;
                end
                if strcmp(connection,'DiffPts')
                    Ypt(cont,cont,:) = ( Y_Sw(1,1,:).*Y_Sw(2,2,:) - Y_Sw(1,2,:).^2 ) ./ ( Y_Sw(1,1,:) + 2*Y_Sw(1,2,:) + Y_Sw(2,2,:) );
                    Gamma_BB(cont,cont,:) = (1./Ypt(cont,cont,:)-Z_ref)./(1./Ypt(cont,cont,:)+Z_ref);
                elseif strcmp(connection,'Conventional')
                    Gamma_BB(cont,cont,:) = S11_Sw - S12_Sw.*S21_Sw./(1+S22_Sw);
                    Ypt(cont,cont,:) = myS2Y(Gamma_BB);
                end
            end
        end

        Gamma = [zeros(N_A,N_A,nfreqs), zeros(N_A,N_B,nfreqs); zeros(N_B,N_A,nfreqs), Gamma_BB];

        %% calculate the port loaded reduced S matrix and the current across the elements connected at the ports
        V_in = V_pt*ones(N_A,nfreqs); % 1V Voltage excitation connected at feeding (external) port
        V = [V_in; zeros(N_B,nfreqs)];
        Icurr = zeros(totN_ports,nfreqs);
        %S_L{contInd} = zeros(N_A,N_A,nfreqs);

        for cont = 1:nfreqs
            %S_L{contInd}(:,:,cont) = S_AA(:,:,cont) + S_AB(:,:,cont) * ((Idmat - Gamma_BB(:,:,cont)*S_BB(:,:,cont)) \ (Gamma_BB(:,:,cont) * S_BA(:,:,cont)));
            S_L(contInd,cont,:,:) = S_AA(:,:,cont) + S_AB(:,:,cont) * ((Idmat - Gamma_BB(:,:,cont)*S_BB(:,:,cont)) \ (Gamma_BB(:,:,cont) * S_BA(:,:,cont)));
            
            I_A = (Y_AA(:,:,cont) - Y_AB(:,:,cont)* ((Y_BB(:,:,cont) + Ypt(:,:,cont))\ Y_BA(:,:,cont)) )*V_in(:,cont);
            V_B = -(Y_BB(:,:,cont)+Ypt(:,:,cont))\ Y_BA(:,:,cont)* V_in(:,cont);
            % This calculation is unstable when using ideal switches,
            % temporarily disable nearlySingularMatrix warning 
            w = warning('query','last');
            id = w.identifier;
            warning('off',id);
            
            I_B = -Ypt(:,:,cont)*V_B;
            Icurr(:,cont) = [I_A; I_B];
            V(N_A+1:end,cont) = V_B;
        end
%        Cost_S_L(contInd) = max(20*log10(abs(S_L{contInd}(:,:,(freqsComm>f_min & freqsComm<f_max))))); % This considers input port match only
        Cost_S_L(contInd) = max(20*log10(abs(S_L(contInd,(freqsComm>f_min & freqsComm<f_max),:,:)))); % This considers input port match only
        Cost_Sw_I(contInd) = max(max(abs(Icurr(N_A+1:N_A+N_B,(freqsComm>f_min & freqsComm<f_max))))); 

        %% Calculate the combined E field
        a_A = ones(N_A,nfreqs);
        a_B = zeros(N_B,nfreqs);

        for cont = 1:nfreqs
            a_B(:,cont) = (Idmat - Gamma_BB(:,:,cont)*S_BB(:,:,cont))\ Gamma_BB(:,:,cont)*S_BA(:,:,cont)*a_A(cont);
        end

        aA = a_A(1,freqsComm == freqPat);  % Patterns must be exported for a fixed single frequency, this validation must be properly handled 
        comb_radPat.theta = aA*rad_E_field(1).E_Field.E_theta;
        comb_radPat.phi = aA*rad_E_field(1).E_Field.E_phi;

        for cont = 2:N_A
            aA = a_A(cont,freqsComm == freqPat);  
            comb_radPat.theta = comb_radPat.theta + aA*rad_E_field(cont).E_Field.E_theta;
            comb_radPat.phi = comb_radPat.phi + aA*rad_E_field(cont).E_Field.E_phi;
        end
        % Sum the contributions of the loaded ports 
        for cont = 1:N_B
            aB = a_B(cont,freqsComm == freqPat);
            comb_radPat.theta  = comb_radPat.theta + aB*rad_E_field(cont+N_A).E_Field.E_theta;
            comb_radPat.phi = comb_radPat.phi + aB*rad_E_field(cont+N_A).E_Field.E_phi;
        end
        % Calculate the combined radiation pattern 
        %eff = (1-abs(S_L{contInd}(:,:,freqsComm == freqPat)).^2);
        eff = (1-abs(S_L(contInd,freqsComm == freqPat,:,:)).^2);
        mag_comb_radPat = eff*(1/60)*(1/P_ant.Psource)*(comb_radPat.theta.*conj(comb_radPat.theta) + comb_radPat.phi.*conj(comb_radPat.phi));
        mag_comb_radPat_dB = 10*log10(mag_comb_radPat);

        [maxPatternPhi, ind_thetaMax] = max(mag_comb_radPat_dB,[],1);
        [maxPatternVal, ind_phiMax] = max(maxPatternPhi);

        thetaMax = unique(thetaGrid(ind_thetaMax(ind_phiMax),:));
        phiMax = unique(phiGrid(:,ind_phiMax));
        AngDev(contInd,:) = [thetaMax,phiMax];

        Cost_AngDev(contInd) = sqrt((targetTheta-thetaMax).^2 + (targetPhi - phiMax).^2);
        %polDev(contInd)
    end
    % Restore warnings!
    warning('on',id);

    normCost_SL = normCost(Cost_S_L,targetS,maxSval);
    normCost_Sw_I = normCost(Cost_Sw_I,targetI,maxIval);
    normCost_AngDev = normCost(Cost_AngDev,0,maxAngDev);
	
    if strcmp(Combination,'convex') 
        varargout{1} = [normCost_SL, normCost_Sw_I, normCost_AngDev]*[0.3;0;0.7]; % Combining
    elseif strcmp(Combination,'none') 
        varargout{1} = [normCost_SL, normCost_Sw_I, normCost_AngDev];
    end
    
    if nargout == 2 
        varargout{2}.worstS_L = Cost_S_L;
        varargout{2}.AngDev = AngDev*180/pi;
        varargout{2}.SL = S_L;
        varargout{2}.freqs = freqsComm;
        %varargout{2}.SwCurr = SwCurr;
    elseif nargout == 3      
        if nInds == 1
            varargout{2}.worstS_L = Cost_S_L;
            varargout{2}.AngDev = AngDev*180/pi;
            varargout{2}.SL = S_L;
            varargout{2}.freqs = freqsComm;
            varargout{3}.Gain = maxPatternVal;
            varargout{3}.phi = phiGrid*180/pi;
            varargout{3}.theta = thetaGrid*180/pi;
            varargout{3}.Pattern = mag_comb_radPat_dB;
        else 
            error('Field results only available for one individual evaluation');
        end
    end
end

%% Normalized cost calculation
function out = normCost(rawCost,minVal,maxVal)
	rawCost = max(min(rawCost,maxVal),minVal);
    out = (rawCost - minVal)/(maxVal - minVal);
end

% Suited for the case when varargin is implemented
function [nInds,nSws,dataroute,SW_S_pars_path_On,SW_S_pars_path_Off,...
    ExhaustiveSearch,Combination,idealSws,connection,maxSval,targetS,targetI,maxIval,...
    targetTheta,targetPhi,maxAngDev,f_min,f_max,patFreq] =...
    inputsRead(inds,dataPath,options,varargin)

    [nInds, nSws] = size(inds);

    dataroute = dataPath.SmatDir;
	SW_S_pars_path_On = dataPath.SwsDir_On;
    SW_S_pars_path_Off = dataPath.SwsDir_Off; 

    ExhaustiveSearch = options.Exhaustive;
	Combination = options.Combination;

    idealSws = options.idealSws; 
    connection = options.connection; 
    maxSval = options.maxSval; 
    targetS = options.targetS;
    targetI = options.targetI;
    maxIval = options.maxIval;
    targetTheta = options.targetTheta*pi/180;
    targetPhi = options.targetPhi*pi/180;
    maxAngDev = options.maxAngDev*pi/180; 
    f_min = options.f_min; 
    f_max = options.f_max; 
    patFreq = options.patFreq; 
end