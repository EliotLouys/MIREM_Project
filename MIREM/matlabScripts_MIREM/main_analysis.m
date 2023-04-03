


userName  = 'el1';
userInfo  = UserSessionInfo_MIREM(userName);

datafilenames = ls(userInfo.dataDir);
datafilenames = datafilenames(3:end,:);
scorfilenames = ls(userInfo.scorDir);
scorfilenames = scorfilenames(3:end,:);



numfiles  = length(datafilenames(:,1));

init= cell(1,numfiles);

results= struct ('filename', init, 'results_table', struct);



for i=1:numfiles
    curr_datafilename        = erase(datafilenames(i,:),' ');
    curr_scorfilename        = erase(scorfilenames(i,:),' ');
    results(i).results_table = preprocessing_MIREM(userName, curr_datafilename(1,:), curr_scorfilename(1,:)) ;
    results(i).filename      = datafilenames(i,:);
end



figure(1)
plot(detection_REM.full_time, detection_REM.full_data)
hold on 
for i=1:length(detection_REM.start_index)
    scatter(detection_REM.full_time(detection_REM.start_index(i):detection_REM.stop_index(i)),detection_REM.full_data(detection_REM.start_index(i):detection_REM.stop_index(i)), 'r')
end
hold off

