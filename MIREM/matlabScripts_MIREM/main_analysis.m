


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


% plotting signals and results in fig 1 and plotting the means of GMM
% alongside the peaks in fig2.
figure(1)
plot(results(1).results_table.full_time, results(1).results_table.full_data)
hold on 
for i=1:length(results(1).results_table.start_index)
    scatter(results(1).results_table.full_time(results(1).results_table.start_index(i):results(1).results_table.stop_index(i)),results(1).results_table.full_data(results(1).results_table.start_index(i):results(1).results_table.stop_index(i)), 'r')
end
title(['Night ' results(1).filename])
xlim([23430 23460])
xlabel('Time (in seconds)')
ylabel('Amplitude (in uV)')
legend('Full signal in bipolar', 'detected REM events')
hold off


[pks, locs ]=findpeaks(full_data);

figure(2)
plot(results(1).results_table.full_time, results(1).results_table.full_data)
hold on
scatter(full_time(locs),full_data(locs),'g')
plot([0 results(1).results_table.full_time(length(results(1).results_table.full_time))], [max(GMModel.mu) max(GMModel.mu)])
plot([0 results(1).results_table.full_time(length(results(1).results_table.full_time))], [min(GMModel.mu) min(GMModel.mu)])
title(['Night ' results(1).filename])
xlim([23430 23460])
xlabel('Time (in seconds)')
ylabel('Amplitude (in uV)')
legend('Full signal in bipolar', 'First mean in the GMM model', 'Second mean in the GMM model')
hold off