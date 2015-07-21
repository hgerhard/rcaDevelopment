function plotSNR(mainRcaData,mainNoiseData,rcaSettings,comparisonRcaData,comparisonNoiseData)

% ### add functionality for condition separation
% ### for x-axis labels
% ### for frequency names
% ### for comparison name (if requested)

if nargin<2, error('You must provide both signal and noise data.'); end
if nargin<3, error('You must provide the rcaSettings struct created when your rca data were created.'); end
if nargin<4, plotComparison = false; end
if nargin==5, plotComparison = true; end

nFreqs = length(rcaSettings.freqsToUse);

avgRcaData = aggregateData(mainRcaData,rcaSettings);
avgNoise1Data = aggregateData(mainNoiseData.lowerSideBand,rcaSettings);
avgNoise2Data = aggregateData(mainNoiseData.higherSideBand,rcaSettings);

snrMain = computeSnr(avgRcaData,avgNoise1Data,avgNoise2Data);

if plotComparison
    avgCompData = aggregateData(comparisonRcaData,rcaSettings);
    avgCompNoise1Data = aggregateData(comparisonNoiseData.lowerSideBand,rcaSettings);
    avgCompNoise2Data = aggregateData(comparisonNoiseData.higherSideBand,rcaSettings);
    
    snrComparison = computeSnr(avgCompData,avgCompNoise1Data,avgCompNoise2Data);
end

for rc=1:rcaSettings.nComp
    figure;
    set(gca,'Color','w');
    for f=1:nFreqs
        subplot(nFreqs,1,f); hold on
        plot(snrMain(:,f,rc),'--ok','MarkerFaceColor','k');
        dataLabels = sprintf('RC%d',rc);
        if plotComparison
            plot(snrComparison(:,f,1),'--or','MarkerFaceColor','r'); % ### make function flexible enough to handle comparison to single RC and multiple RC
            dataLabels = {dataLabels,'Comparison'};
        end
        title(['Frequency ' num2str(rcaSettings.freqsToUse(f))']);
        ylabel('SNR');
        if f==1, hlg=legend(dataLabels); end
    end
    set(hlg,'box','off');
end