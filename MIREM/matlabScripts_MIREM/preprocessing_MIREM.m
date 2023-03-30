function [REM_events_ts, REM_events_data, maxSlopes, minSlopes, eventpks] = preprocessing_MIREM(userName, nameEEG, nameScore, plot)

% Pipeline to preprocess EEG data in MIREM project
%   use as [REM_events_ts, REM_events_data, maxSlopes, minSlopes, eventpks]
%   = preprocessing_MIREM('jb1', '105_NN_Sommeil.edf', '105_NN_Sommeil_scores.csv','no');
%
% JB Eichenlaub, 2023 || jb.eichenlaub@gmail.com (MATLAB2022a)
% 
% input arguments :
% 
% userName  : must be one of the users defined in UserSessionInfo_MIREM.m where all the paths to the data/scoring has been set
% nameEEG   : name of the data edf file of the data and onle the file name; the path has to be set in UserSessionInfo_MIREM.m
% nameScore : name of the scoring csv file of the data and onle the file name; the path has to be set in UserSessionInfo_MIREM.m
% plot      : must be 'yes' or 'no'. If 'yes', the programm will plot the results.
% 
% output arguments :
% 
% REM_events_ts   : vector of size 'number of events detected' where the timestamps where REM events were detected is stored. Should be same size as REM_events_data.
% REM_events_data : vector of size 'number of events detected' where the values of the signal corresponding the timestamps in REM_events_ts is stored. Should be same size as REM_events_ts.
% maxSlopes       : vector of size 'number of events detected' where the maximum slope each event is stored.
% minSlopes       : vector of size 'number of events detected' where the mnimum slope each event is stored.
% eventpks        : vector of size 'number of events detected' where the highest value detected in each event is stored.

%% Verifying the input argument 
if ~strcmp(plot,'yes') & ~strcmp(plot,'no')
    disp( ['invalid argument plot=' plot '; must be yes or no, please use valid argument'])
    maxSlopes        = [NaN] ;
    minSlopes        = [NaN] ;
    eventpks         = [NaN] ;
    REM_events_ts    = [NaN] ;
    REM_events_data  = [NaN] ;
    return;
end

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

data_thresh_pos  = data_REM - threshold_G;       
data_thresh_neg  = data_REM + threshold_G;
[ zX  , zset]    = detectzerocross(full_data);
[ tpX , tpset]   = detectzerocross(data_thresh_pos);                                 % Applying a zero-crossing detection method to detect where the signal crosses the previously estimated threshold
[ tnX , tnset]   = detectzerocross(data_thresh_neg);
time_REM_zcp     = time_REM(tpX);
time_REM_zcn     = time_REM(tnX);
numEpochSlope    = 20 ;
maxSlopes        = [] ;
minSlopes        = [] ;
eventpks         = [] ;
REM_events_ts    = [] ;
REM_events_data  = [] ;
dur_crit         = 4  ;
slope_crit       = 0.001;

for i=1:length(zX)-1
    cross       = [zset(i) zset(i+1)];
    
    if cross(1)==cross(2)
        continue;
    end

    switch cross(1)
        
        case 1  % Evenement au dessus de 0
            
            thresh_crit = time_REM_zcp( full_time(zX(i)) < time_REM_zcp & time_REM_zcp < full_time(zX(i+1)) );  % un threshold crossing dans le zero crossing 
            

            if any(thresh_crit) %Il y a au moins un threshold crossing dans le zero crossing en question

               if cross(2)-cross(1)< dur_crit                                                                   % critère si le zero-crossing fait moins de 4secondes ?
                   curr_set    = tpX( full_time(zX(i)) < time_REM_zcp & time_REM_zcp < full_time(zX(i+1)) );
                   curr_data   = data_REM(  curr_set(1) :  curr_set(length(curr_set))  );
                   eplength    = fix(length(curr_data)/numEpochSlope);
                   ref_poly    = 1 : eplength ;
                   maxs_tmp    = -1000;
                   mins_tmp    = 1000;

                   for j=1 : eplength : length(curr_data)-eplength                                              % Calcul de pente : sur l'intervalle entre deux zeros-crossing, 
                        est_slope  = polyfit( ref_poly, curr_data( j:j+eplength -1), 1  );                      % on adapte un polynome de degré 1 un nombre numEpochSlope
                        curr_slope = est_slope(1);                                                                      % de fois et on prends le coefficient de cette droite comme pente.
                                                                                                                % On garde ensuite la pente max et min detectée temporairement. 
                        if curr_slope > maxs_tmp
                            maxs_tmp = curr_slope;
                        elseif curr_slope < mins_tmp
                            mins_tmp = curr_slope;
                        end

                   end


                   if maxs_tmp > slope_crit || abs(mins_tmp) > slope_crit                                              % On applique le critère de pente sur les pentes min et max
                       maxSlopes       = [maxSlopes maxs_tmp];                                                         % temporairement stockées. Si le critère est rempli, on stocke alors
                       minSlopes       = [minSlopes mins_tmp];                                                         % ces pentes, le maximum atteint sur l'évènement, les timestamps et 
                       eventpks        = [eventpks max(  data_REM(  curr_set(1) :  curr_set(length(curr_set))  )  )];  % data correspondant à l'event dans des vecteurs dédiés à être retournés par la fonction.
                       REM_events_ts   = [REM_events_ts  full_time(zX(i):zX(i+1)) ];
                       REM_events_data = [REM_events_data full_data(zX(i):zX(i+1))] ;
                   end
               end

            end

        case -1  % Evenement en dessous de zero
            
            thresh_crit = time_REM_zcn( full_time(zX(i)) < time_REM_zcn & time_REM_zcn < full_time(zX(i+1)) ); % un threshold crossing dans le zero crossing ?
            
            if any(thresh_crit)

                if cross(2)-cross(1)< dur_crit % critère si le zero-crossing fait moins de 4secondes ?
                   curr_set    = tnX( full_time(zX(i)) < time_REM_zcn & time_REM_zcn < full_time(zX(i+1)) );
                   curr_data   = data_REM(  curr_set(1) :  curr_set(length(curr_set))  );
                   eplength    = fix(length(curr_data)/numEpochSlope);
                   ref_poly    = 1 : eplength;
                   maxs_tmp    = -1000;
                   mins_tmp    = 1000;
                   
                   for j=1 : eplength : length(curr_data)-eplength                                                    % Calcul de pente : sur l'intervalle entre deux zeros-crossing, 
                        est_slope  = polyfit( ref_poly, curr_data( j:j+eplength -1 ), 1  );                           % on adapte un polynome de degré 1 un nombre numEpochSlope
                        curr_slope = est_slope(1);                                                                            % de fois et on prends le coefficient de cette droite comme pente.
                                                                                                                      % On garde ensuite la pente max et min detectée temporairement. 
                        if curr_slope > maxs_tmp
                            maxs_tmp = curr_slope;
                        elseif curr_slope < mins_tmp
                            mins_tmp = curr_slope;
                        end

                   end


                   if maxs_tmp > slope_crit || abs(mins_tmp) > slope_crit                                             % On applique le critère de pente sur les pentes min et max
                       maxSlopes       = [maxSlopes maxs_tmp];                                                        % temporairement stockées. Si le critère est rempli, on stocke alors
                       minSlopes       = [minSlopes mins_tmp];                                                        % ces pentes, le maximum atteint sur l'évènement, les timestamps et 
                       eventpks        = [eventpks max(  data_REM(  curr_set(1) :  curr_set(length(curr_set))  )  )]; % data correspondant à l'event dans des vecteurs dédiés à être retournés par la fonction.
                       REM_events_ts   = [REM_events_ts  full_time(zX(i):zX(i+1)) ];
                       REM_events_data = [REM_events_data full_data(zX(i):zX(i+1))] ;
                   end
                end
            
                

            end
    end

end


    if strcmp(plot,'yes')                                                              %plotting the results depending on the input argument 'plot'
        figure()
        plot(full_time, full_data)
        hold on 
        scatter(REM_events_ts,REM_events_data,36,'red','filled')
        title('Data REM along with the REM events detected')
        xlabel('Time in sec')
        ylabel('Amplitude in uV')
        hold off
    
    
    end


end


