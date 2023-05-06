nodesInFull = 149*ones(1,9);
nodesInComp = [81 63 71 68 73 56 49 43 37];
compRate = (nodesInFull - nodesInComp)./nodesInFull*100;

compNodeP = nodesInComp./nodesInFull;
timeComp = [1208.4 758 962.7 883.1 1017.8 598.9 458.5 353.1 261.5];
timeFull = [5310 5231 5302 5286 5297 5120 5118 5093 4903];
timeRatio = timeFull./timeComp;
ff_full = [22.3 12.3 10.3 17 12 9 8 13 15];
ff_comp = [23.2 11.9 9.8 16.9 13.2 8.9 6.7 12.4 13.5];
qualityRatio = ff_full./ff_comp;
[compRate,I] = sort(compRate);
timeRatio = timeRatio(I);
qualityRatio = qualityRatio(I);
figure(1)
subplot(2,1,1);
plot(compRate,timeRatio,'r*--');
title("time Ratio");
xlabel('Compress Rate');
ylabel("Time Ratio(T/T')");
axis([min(compRate),max(compRate),0, 20]);
subplot(2,1,2);
plot(compRate,qualityRatio,'o--');
title("Quality Ratio(ff/ff')");
xlabel('Compress Rate');
ylabel("Quality Ratio(ff/ff')");
axis([min(compRate),max(compRate),0, 2]);

disp("Elapsed time is 123.23213 seconds.")