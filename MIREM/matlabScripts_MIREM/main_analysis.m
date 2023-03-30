


userName  = 'el1';
userInfo  = UserSessionInfo_MIREM(userName);

datafilenames = ls(userInfo.dataDir);
datafilenames = datafilenames(3:end,:);
scorfilenames = ls(userInfo.scorDir);
scorfilenames = scorfilenames(3:end,:);



numfiles  = length(datafilenames(:,1));

init= cell(1,numfiles);

results= struct ('filename' , init , 'REM_events_timestamps' , init , 'REM_events_data' , init , 'maxSlopes' , init , 'minSlopes' , init , 'eventpeak' , init );



for i=1:numfiles
    curr_datafilename                = erase(datafilenames(i,:),' ');
    curr_scorfilename                = erase(scorfilenames(i,:),' ');
    [REM_events_ts, REM_events_data, maxSlopes, minSlopes, eventpks] = preprocessing_MIREM(userName, curr_datafilename, curr_scorfilename, 'no');
    results(i).REM_events_timestamps = REM_events_ts;
    results(i).REM_events_data       = REM_events_data;
    results(i).maxSlopes             = maxSlopes;
    results(i).minSlopes             = minSlopes;
    results(i).eventpeak             = eventpks;
    results(i).filename              = datafilenames(i,:);
end