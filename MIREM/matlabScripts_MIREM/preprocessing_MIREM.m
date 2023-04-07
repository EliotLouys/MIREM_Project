function [detected_REM_table] = preprocessing_MIREM(userName, nameEEG, nameScore)

% Pipeline to preprocess EEG data in MIREM project
%   use as [REM_events_ts, REM_events_data, maxSlopes, minSlopes, eventpks]
%   = preprocessing_MIREM('jb1', '105_NN_Sommeil.edf', '105_NN_Sommeil_scores.csv');
%
% JB Eichenlaub, 2023 || jb.eichenlaub@gmail.com (MATLAB2022a)
% 
% 
% input arguments :
% 
% userName  : must be one of the users defined in UserSessionInfo_MIREM.m where all the paths to the data/scoring has been set
% nameEEG   : name of the data edf file of the data and onle the file name; the path has to be set in UserSessionInfo_MIREM.m
% nameScore : name of the scoring csv file of the data and onle the file name; the path has to be set in UserSessionInfo_MIREM.m
% 
% 
% output arguments :
% 
% detected_REM_table : a structure with 7 fields that returns the following
% fields:
% 
% start_index   : vector of size 'number of events detected' where the starting indexes of each REM are stored.
% stop_index    : vector of size 'number of events detected' where the stopping indexes of each REM are stored.
% max_slope     : vector of size 'number of events detected' where the maximum slope each event is stored. The units is uV/ms.
% min_slope     : vector of size 'number of events detected' where the minimum slope each event is stored. The units is uV/ms
% peak          : vector of size 'number of events detected' where the highest value detected in each event is stored. The unit is uV.
% duration      : vector of size 'number of events detected' where the duration of each event is stored. The unit is seconds.
% full_data     : vector where the original signal in bipolar resampled at 100Hz is stored.
% full_time     : vector where the time vector resampled at 100Hz is stored.



%% Parameters
paramsPrepo.highPassFreq  = 0.1;
paramsPrepo.highPassOrder = 4;

paramsPrepo.lowPassFreq  = 3;
paramsPrepo.lowPassOrder = 4;

downsampling_frequency   = 100;

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
cfg.resamplefs = downsampling_frequency;
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

%%  loading some signals we will need later.
full_data      = data_EOG_bi.trial{1};
full_time      = data_EOG_bi.time{1};
data_REM       = full_data(1, find(vectorREM == 1));
% time_REM       = full_time(1, find(vectorREM == 1));


% A mask of incoherent values in the signal that should not be taken into
% account for the following operations. 
iv_thresh                       = 500;
iv_mask                         = ones( 1, length(full_data) );
iv_mask( full_data > iv_thresh) = 0;


%%
% Step 1: Define the amplitude threshold using the Gaussian Mixture Model to
% fit the peaks of the signal.

abs_signal         = abs(data_REM);
pks                = findpeaks(abs_signal);
pks(pks>iv_thresh) = [];
GMModel            = fitgmdist(transpose(pks),2);
threshold_G        = max(GMModel.mu);


%%
% Step 2: First, we will detect every zero crossing in the signal.
% Analyzing every zero crossing detected that are during REM sleep stages,
% we will apply 3 criteriums on the potential REM events and if the three 
% criteria are met, we will consider the current event as a REM. 



% Defining our criterias values 
dur_crit         = 4  ;   
slope_crit       = 10 ;
 
% Applying a zero-crossing detection method to detect where the signal crosses the previously estimated threshold
[ zX  , zset]    = detectzerocross(full_data);      

% Keeping only the zero-crossings during REM stages of sleep
candidates         = intersect(zX,find(vectorREM));
candidatesset      = zset(ismember(zX,find(vectorREM))); 

% Building the structure where we will return the informations about the
% detected REM events. 
detected_REM_table = struct('start_index', {0}, 'stop_index', {0}, 'max_slope', {0}, 'min_slope', {0}, 'peak', {0}, 'duration', {0}, 'full_time', full_time, 'full_data', full_data);



% Loop over the candidates and applying the criterias to every candidates.
for i=1:length(candidates)-1

    % Checking that the two zero crossing aren't both onsets or offsets. 
    cross       = [candidatesset(i) candidatesset(i+1)];
    if cross(1)==cross(2)
        continue;
    end
    
    if cross(1)==1
        set='positive crossing';
    
    elseif cross(1)==-1
        set='negative crossing';
    end

    % Keeping just the current data to apply our criterias.
    curr_data     = full_data(candidates(i):candidates(i+1));
    curr_iv_mask  = iv_mask(candidates(i):candidates(i+1));

    % Applying the threshold criteria with the threshold defined earlier
    % with the GMM method.
    thresh_crit   = max( abs( curr_data ) );
    if thresh_crit > threshold_G

        % Applying the time criteria, a REM event can't last longer than 4
        % seconds

        duration = full_time(candidates(i+1))-full_time(candidates(i));
        if duration < dur_crit

            % Computing the slopes and applying the slope criteria: the
            % signal must have a maximum or minimum slope above 1uV/ms
            % (absolute value).
            curr_iv_data = curr_data(curr_iv_mask==1);
            slope_vector = diff(curr_iv_data);

            maxSlope     = max(slope_vector);
            minSlope     = min(slope_vector);

            if maxSlope > slope_crit || abs(minSlope) > slope_crit 

                % All criterias are met we store the useful parameters of
                % each event in the results structure. 
                detected_REM_table.max_slope( length(detected_REM_table.max_slope)+1)     = maxSlope;
                detected_REM_table.min_slope( length(detected_REM_table.min_slope)+1)     = minSlope;
                detected_REM_table.start_index( length(detected_REM_table.start_index)+1) = candidates(i);
                detected_REM_table.stop_index( length(detected_REM_table.stop_index)+1)   = candidates(i+1);
                detected_REM_table.duration( length(detected_REM_table.duration)+1)       = duration;

                % Deciding if we must store a maximum or a minimum peak
                % depending on the positive or negative crossing.
                if strcmp(set,'positive crossing')
                    detected_REM_table.peak( length(detected_REM_table.peak)+1)           = max(curr_data);
                elseif strcmp(set,'negative crossing')
                    detected_REM_table.peak( length(detected_REM_table.peak)+1)           = min(curr_data);
                end
         
            end

        end

    end

end   

% Must delete first value of each fields used to create the structure. 
detected_REM_table.max_slope(1)   = [];
detected_REM_table.min_slope(1)   = [];
detected_REM_table.start_index(1) = [];
detected_REM_table.stop_index(1)  = [];
detected_REM_table.peak(1)        = [];
detected_REM_table.duration(1)    = [];

%TODO implanter un masque des valeurs absurdes et enlever ces valeurs
%absurdes au moment de regarder les pentes et les threshold j'ai pas eu
%trop le temps d'y réfléchir pour l'instant oupsi. 

end


