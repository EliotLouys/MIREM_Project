function [detected_REM_table] = preprocessing_MIREM(userName, nameEEG, nameScore)

% Pipeline to preprocess EEG data in MIREM project
%   use as [REM_events_ts, REM_events_data, maxSlopes, minSlopes, eventpks]
%   = preprocessing_MIREM('jb1', '105_NN_Sommeil.edf', '105_NN_Sommeil_scores.csv');
%
% JB Eichenlaub, 2023 || jb.eichenlaub@gmail.com (MATLAB2022a)
% 
% input arguments :
% 
% userName  : must be one of the users defined in UserSessionInfo_MIREM.m where all the paths to the data/scoring has been set
% nameEEG   : name of the data edf file of the data and onle the file name; the path has to be set in UserSessionInfo_MIREM.m
% nameScore : name of the scoring csv file of the data and onle the file name; the path has to be set in UserSessionInfo_MIREM.m

% 
% output arguments :
% 
% REM_events_ts   : vector of size 'number of events detected' where the timestamps where REM events were detected is stored. Should be same size as REM_events_data.
% REM_events_data : vector of size 'number of events detected' where the values of the signal corresponding the timestamps in REM_events_ts is stored. Should be same size as REM_events_ts.
% maxSlopes       : vector of size 'number of events detected' where the maximum slope each event is stored.
% minSlopes       : vector of size 'number of events detected' where the mnimum slope each event is stored.
% eventpks        : vector of size 'number of events detected' where the highest value detected in each event is stored.



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
scorDir     = userInfo.scorDir;

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

%resampling at 100Hz
cfg            = [];
cfg.resamplefs = 100;
resampled_data = ft_resampledata(cfg, data_EOG);

% bipolar
cfg             = [];
cfg.channel     = 'all';
cfg.reref       = 'yes';
cfg.refchannel  = resampled_data.label{2};
data_EOG_bi     = ft_preprocessing(cfg, resampled_data);


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
scoreFile               = [scorDir nameScore];
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


%TODO utiliser ismember ou l'autre jsp a regarder, c'est bien aussi
%intersect je crois. 
%TODO sur chaque zero crossing, on regarde les zeros crossings qui sont
%dans les periodes REM. Sur les zeros crossing de periode REM: on va
%regarder si on a au moins une valeur (en abs) qui dépasse le threshold
%voulu. Ensuite on regarde si la durée est <4s et ensuite on calcule pente:
%pour le zero crossing, on va faire un intersect entre le masque et les
%indices de zc trouvés et ducoup on garde juste les indices qui nous
%intéressent. 




dur_crit         = 4  ;
slope_crit       = 10 ;


[ zX  , zset]    = detectzerocross(full_data);            % Applying a zero-crossing detection method to detect where the signal crosses the previously estimated threshold
candidates         = intersect(zX,find(vectorREM));
candidatesset      = zset(ismember(zX,find(vectorREM)));
detected_REM_table = struct('start_index', {0}, 'stop_index', {0}, 'max_slope', {0}, 'min_slope', {0}, 'peak', {0}, 'full_time', full_time, 'full_data', full_data);


for i=1:length(candidates)-1
    cross       = [candidatesset(i) candidatesset(i+1)];
    
    if cross(1)==cross(2)
        continue;
    end
    
    if cross(1)==1
        set='positive crossing';
    
    elseif cross(1)==-1
        set='negative crossing';
    end


    curr_data     = full_data(candidates(i):candidates(i+1));
    thresh_crit   = max( abs( curr_data ) );

    if thresh_crit > threshold_G
        time_crit = full_time(candidates(i+1))-full_time(candidates(i));

        if time_crit < 4

            slope_vector= diff(curr_data);
            maxSlope= max(slope_vector);
            minSlope=min(slope_vector);


            if maxSlope > slope_crit || abs(minSlope) > slope_crit 
                detected_REM_table.max_slope( length(detected_REM_table.max_slope)+1)     = maxSlope;
                detected_REM_table.min_slope( length(detected_REM_table.min_slope)+1)     = minSlope;
                detected_REM_table.start_index( length(detected_REM_table.start_index)+1) = candidates(i);
                detected_REM_table.stop_index( length(detected_REM_table.stop_index)+1)   = candidates(i+1);

                if strcmp(set,'positive crossing')
                    detected_REM_table.peak( length(detected_REM_table.peak)+1)           = max(curr_data);
                elseif strcmp(set,'negative crossing')
                    detected_REM_table.peak( length(detected_REM_table.peak)+1)           = min(curr_data);
                end
         
            end

        end

    end

end   

detected_REM_table.max_slope(1)   = [];
detected_REM_table.min_slope(1)   = [];
detected_REM_table.start_index(1) = [];
detected_REM_table.stop_index(1)  = [];
detected_REM_table.peak(1)        = [];




end


