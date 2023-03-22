function [proportion]=preprocessing_MIREM(userName, nameEEG, nameScore)

% Pipeline to preprocess EEG data in MIREM project
%   use as preprocessing_MIREM('jb1', '105_NN_Sommeil.edf', '105_NN_Sommeil_scores.csv')
%
% JB Eichenlaub, 2023 || jb.eichenlaub@gmail.com (MATLAB2022a)

%% Parameters
paramsPrepo.highPassFreq  = 0.1;
paramsPrepo.highPassOrder = 4;

paramsPrepo.lowPassFreq  = 3;
paramsPrepo.lowPassOrder = 4;

%% Define main paths and folders
% collect users' and sessions' info
userInfo    = UserSessionInfo_MIREM(userName);
analysisDir = userInfo.analysisDir;
dataDir     = userInfo.dataDir;
% Create, if necessary, analysis folder
if ~exist(analysisDir, 'dir')
    disp(['Making directory ' analysisDir])
    mkdir(analysisDir);
end

%% preprocessing EOG

% open the file in FieldTrip
cfg         = [];
cfg.dataset = [dataDir nameEEG];
data        = ft_preprocessing(cfg);

% keep only EOG
cfg         = [];
cfg.channel = {'eog'};
data_EOG    = ft_selectdata(cfg, data);

% high-pass filter
cfg            = [];
cfg.hpfilter   = 'yes';
cfg.hpfreq     = paramsPrepo.highPassFreq;
cfg.hpfiltord  = paramsPrepo.highPassOrder;
cfg.hpfilttype = 'but';
cfg.hpfiltdir  = 'twopass';
data_EOG       = ft_preprocessing(cfg, data_EOG);

% low-pass filter
cfg            = [];
cfg.lpfilter   = 'yes';
cfg.lpfreq     = paramsPrepo.lowPassFreq;
cfg.lpfiltord  = paramsPrepo.lowPassOrder;
cfg.lpfilttype = 'but';
cfg.lpfiltdir  = 'twopass';
data_EOG       = ft_preprocessing(cfg, data_EOG);

% bipolar
cfg             = [];
cfg.channel     = 'all';
cfg.reref       = 'yes';
cfg.refchannel  = data_EOG.label{2};
data_EOG_bi     = ft_preprocessing(cfg, data_EOG);

% keep only final bipolar
cfg               = [];
cfg.channel       = data_EOG.label{1};
data_EOG_bi       = ft_selectdata(cfg, data_EOG_bi);
data_EOG_bi.label = 'EOG_bi';

% main info
numSamples = size(data_EOG_bi.trial{1}, 2);
fsample    = data_EOG_bi.fsample;

%% 

% open score file
scoreFile               = [dataDir nameScore];
opts                    = detectImportOptions(scoreFile);
opts.DataLines          = 3;
opts.VariableNamesLine  = 2;
opts.VariableNamingRule = 'preserve';
score                   = readtable(scoreFile, opts);
numEpochs               = height(score);

% define REM samples
vectorREM = zeros(1, numSamples);
for ep=1:numEpochs
    if isequal(score.stage{ep}, 'REM')
        idxStart_s                     = score.start(ep) * fsample;
        idxEnd_s                       = score.end(ep) * fsample;
        vectorREM(idxStart_s:idxEnd_s) = 1;
    end
end

%%  
full_data      = data_EOG_bi.trial{1};
full_time      = data_EOG_bi.time{1};
data_REM       = full_data(1, find(vectorREM == 1));
time_REM       = full_time(1, find(vectorREM == 1));


%%
% Step 1: Define amplitude threshold using the Gaussian Mixture Model to
% fit the peaks of the signal.
abs_signal  = abs(data_REM);
pks         = findpeaks(abs_signal);
pks(pks>500)= [];
GMModel     = fitgmdist(transpose(pks),2);
threshold_G = max(GMModel.mu);


%%
% Step 2: Define zero crossing: 

data_thresh_pos  = data_REM - threshold_G;       
data_thresh_neg  = data_REM + threshold_G;
[ zX  , zset]    = detectzerocross(full_data);
[ tpX , tpset]   = detectzerocross(data_thresh_pos);                                 % Applying a zero-crossing detection method to detect where the signal crosses the previously estimated threshold
[ tnX , tnset]   = detectzerocross(data_thresh_neg);
time_REM_zcp     = time_REM(tpX);
time_REM_zcn     = time_REM(tnX);


for i=1:length(zX)
    cross       = [zset(i) zset(i+1)];
    
    if cross(1)==cross(2)
        continue
    end

    switch cross(1)
        
        case 1   % Evenement au dessus de 0
            
            tmp= time_REM_zcp( full_time(zX(i)) < time_REM_zcp & time_REM_zcp < full_time(zX(i+1)) );

            if any(tmp)
                'lolzp'
            end

        case -1  % Evenement en dessous de zero
            
            tmp= time_REM_zcn( full_time(zX(i)) < time_REM_zcn & time_REM_zcn < full_time(zX(i+1)) );
            if any(tmp)
                'lolzn'
            end

                                %Calcul de peak rise/fall slope duration et conditions



    end

end





