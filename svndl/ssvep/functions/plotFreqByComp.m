function plotFreqByComp(mainData,comparisonData)

figure
for f=1:nFreqs
    for c=1:nComp
        %subplot(nFreqs,nComp,(f-1)*nComp+c); hold on
        subplot(nFreqs,nComp+1,(f-1)*(nComp+1)+c); hold on
        plot(ampBins(:,f,c),'--ok','MarkerFaceColor','k');
        plot(noise1AmpBins(:,f,c),'--or','MarkerFaceColor','r');
        plot(noise2AmpBins(:,f,c),'--og','MarkerFaceColor','g');
        
        if f==1, title(['Component ' num2str(c)']); end
        if c==1, ylabel(['Frequency ' num2str(freqsToUse(f))']); end
        if f==1 && c==1, hlg=legend('signal','noise1','noise2'); end
        
        subplot(nFreqs,nComp+1,f*(nComp+1)); hold on
        plot(ozAmpBins(:,f,c),'--ok','MarkerFaceColor','k');
        plot(ozNoise1AmpBins(:,f,c),'--or','MarkerFaceColor','r');
        plot(ozNoise2AmpBins(:,f,c),'--og','MarkerFaceColor','g');
        if f==1, title('Electrode Oz'); end
        
        
    end
end
set(hlg,'box','off');





