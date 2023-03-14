function preprocessing_MIREM(userName, nameEEG, nameScore)

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
cfg.channel       = {'EOG 1'};
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
data_REM    = data_EOG_bi.trial{1}(1, find(vectorREM == 1));
times       = data_EOG_bi.time{1}(1, find(vectorREM == 1));


%%
% Step 1: Define amplitude threshold using the Gaussian Mixture Model to
% fit the peaks of the signal.
abs_signal  = abs(data_REM);
pks         = findpeaks(abs_signal);
pks(pks>500)= [];
GMModel     = fitgmdist(transpose(pks),2);
threshold_G = max(GMModel.mu);


%%
% Step 2: Define zero crossing: zX see https://www.mathworks.com/matlabcentral/answers/267222-easy-way-of-finding-zero-crossing-of-a-function

zci         = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); 
data_thresh = data_REM - threshold_G;
zX          = zci(data_thresh);



candidates=[];
for i=1:length(zX)-1
    crit     = times(zX(i+1))-times(zX(i));
    if crit <=4000
        candidates = [candidates [times(zX(i)); times(zX(i+1))]];  %stores the indexes of the data_REM vector where zero crossings happen 2 by 2

    end

end


%%
% Step3: Define and use the other criterias on the candidates 

candidates2=[];
for i=1:length(candidates)
    slope_dur =  cast( (candidates(2,i) - candidates(1,i))*100 , 'uint32');
    ind       =  cast(  candidates(1,i) , 'uint32' );
    crit      = (data_REM(  ind +slope_dur  )  -  data_REM( ind ))/ cast( slope_dur,'double');
    if crit<0.001
        candidates2 = [candidates2 , candidates(:,i)];
    end
end


total_time=0;
for i=1:length(candidates2)
    total_time=total_time + ( candidates2(2,i) - candidates2(1,i) );
end
total_time=total_time/3600;





time_REM= length(data_REM)/(1000*60*60);
proportion= total_time*100/time_REM;




