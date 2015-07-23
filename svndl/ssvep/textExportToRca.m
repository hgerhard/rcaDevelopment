function [cellData,indF,indB,noiseCell1,noiseCell2,freqsAnalyzed,binLevels]=textExportToRca(pathname,binsToUse,freqsToUse,channelsToUse,dataType,conditionNumbers)

if nargin<6, selectedConditions = false; elseif nargin==6 selectedConditions = true; end
if nargin<5, dataType='RLS'; end
if nargin<4 || (nargin>3 && isempty(channelsToUse)), channelsToUse=1:128; end;
if nargin<3, freqsToUse=1; fprintf('Using frequency 1 only.'); end;
if nargin<2, binsToUse=0; fprintf('Using average across bins \n'); end;
if nargin<1, error('At least one input argument required'); end

nFreqs=numel(freqsToUse);
nChannels=numel(channelsToUse);
nBins=numel(binsToUse);

if ~strcmp(pathname(end),'/'), pathname=cat(2,pathname,'/'); end

switch dataType
    case {'DFT','RLS'}
        filenames=dir([pathname sprintf('%s*.txt',dataType)]);
    otherwise
        error('dataType must be ''RLS'' or ''DFT''');
end
nFilenames=numel(filenames);
if nFilenames<1
    error('No %s files were found in %s',dataType,pathname);
end

useFile = ones(1,nFilenames);
if selectedConditions
    for f = 1:nFilenames
        if ~ismember(str2double(filenames(f).name(6:8)),conditionNumbers)
            useFile(f) = 0;
        end
    end
    nFilenames = sum(useFile);
end

cellData=cell(nFilenames,1);
noiseCell1=cell(nFilenames,1);
noiseCell2=cell(nFilenames,1);

fcnt = 1;
for f=1:length(filenames)
    if useFile(f)
        
        [~,freqsAnalyzed,binLevels,data]=getSweepDataFlex(fullfile(pathname,filenames(f).name),channelsToUse);
        
        if isempty(data)
            error('No data found in %s',filenames(f).name);
        end
        
        %% check congruency of trials
        trials=data(:,1);
        trialInds=unique(trials(:));
        nTrials=numel(trialInds);
        nEntriesPerTrialInd=zeros(nTrials,1);
        for tr=1:nTrials
            trialData=data(trials==tr,:);
            nEntriesPerTrialInd(tr)=numel(trialData(:));
        end
        
        % only keep trials in which the number of entries is congruent
        trialIndsToKeep=trialInds(nEntriesPerTrialInd==mode(nEntriesPerTrialInd));
        nTrialsToKeep=numel(trialIndsToKeep);
        
        %% organize in data
        eeg=zeros(nFreqs*nBins*2,nChannels,nTrialsToKeep);
        noise1=zeros(nFreqs*nBins*2,nChannels,nTrialsToKeep);
        noise2=zeros(nFreqs*nBins*2,nChannels,nTrialsToKeep);
        for tr=1:nTrialsToKeep
            
            thisTrialInd=trialIndsToKeep(tr);
            trialData=data(trials==thisTrialInd,:);
            trialChannels=trialData(:,2);
            trialFreqs=trialData(:,3);
            trialBins=trialData(:,4);
            
            % freqs x bins
            for ch=1:nChannels
                theseReals=trialData(trialChannels==ch & ismember(trialFreqs,freqsToUse) & ismember(trialBins,binsToUse),5);
                theseImags=trialData(trialChannels==ch & ismember(trialFreqs,freqsToUse) & ismember(trialBins,binsToUse),6);
                
                % noise 1
                theseNoiseReals1=trialData(trialChannels==ch & ismember(trialFreqs,freqsToUse) & ismember(trialBins,binsToUse),7);
                theseNoiseImags1=trialData(trialChannels==ch & ismember(trialFreqs,freqsToUse) & ismember(trialBins,binsToUse),8);
                
                % noise 2
                theseNoiseReals2=trialData(trialChannels==ch & ismember(trialFreqs,freqsToUse) & ismember(trialBins,binsToUse),9);
                theseNoiseImags2=trialData(trialChannels==ch & ismember(trialFreqs,freqsToUse) & ismember(trialBins,binsToUse),10);
                
                eeg(:,ch,tr)=[theseReals; theseImags];
                noise1(:,ch,tr)=[theseNoiseReals1; theseNoiseImags1];
                noise2(:,ch,tr)=[theseNoiseReals2; theseNoiseImags2];
                
            end
            
        end
        
        cellData{fcnt}=eeg;
        noiseCell1{fcnt}=noise1;
        noiseCell2{fcnt}=noise2;
        fcnt = fcnt + 1;
    end
end

binLevels = binLevels(binsToUse);
freqsAnalyzed = freqsAnalyzed(freqsToUse);

% figure out feature vector indices (this is being done on the last trial
% FIXTHIS) ### 
indF=trialFreqs(trialChannels==1 & ismember(trialFreqs,freqsToUse) & ismember(trialBins,binsToUse));
indB=trialBins(trialChannels==1 & ismember(trialFreqs,freqsToUse) & ismember(trialBins,binsToUse));

end
%%

function [colHdr, freqsAnalyzed, binLevels, dataMatrix]=getSweepDataFlex(datafile, chanToSave)


%% Imports Sweep Data from a Text File
%
% [colHdr, freqsAnalyzed, binLevels, dataMatrix]=GetSweepDataFlex(datafile, chanToSave)
%
%% Inputs:
% datafile     string containing the data file 
% chanToSave   the requested electrode(s)
%
%% Outputs:
% colHdr          is a string with column fields
% freqsAnalyzed   are the individual frequencies of the VEPs ('1F1' etc.)
% binLevels       are the contrasts or stimulus values used
% dataMatrix      matrix containing the desired data

fname=datafile;
%  SetUp import columns and data formats, from JA
%  Asterisk (*) after % exclude selected columns
%

hdrFields = {
    'iSess'         '%s\t'      0 %1
    'iCond'         '%f\t'      1 %2
    'iTrial'        '%f\t'      1 %3
    'iCh'           '%s\t'      1 %4 This becomes %f later in the function
    'iFr'           '%f\t'      1 %5
    'AF'            '%f\t'      1 %6
    'xF1'           '%f\t'      1 %7
    'xF2'           '%f\t'      1 %8
    'Harm'          '%s\t'      2 %9 
    'FK_Cond'       '%f\t'      1 %10
    'iBin'          '%f\t'      1 %11
    'SweepVal'      '%f\t'      1 %12
    'Sr'            '%f\t'      1 %13
    'Si'            '%f\t'      1 %14
    'N1r'           '%f\t'      1 %15
    'N1i'           '%f\t'      1 %16
    'N2r'           '%f\t'      1 %17
    'N2i'           '%f\t'      1 %18
    'Signal'        '%f\t'      1 %19
    'Phase'         '%f\t'      1 %20
    'Noise'         '%f\t'      1 %21
    'StdErr'        '%f\t'      1 %22
    'PVal'          '%f\t'      1 %23
    'SNR'           '%f\t'     2 %24
    'LSB'           '%f\t'     2 %25
    'RSB'           '%f\t'     2 %26
    'UserSc'        '%s\t'     2 %27
    'Thresh'        '%f\t'     2 %28
    'ThrBin'        '%f\t'     2 %29
    'Slope'         '%f\t'     2 %30
    'ThrInRange'    '%s\t'     2 %31
    'MaxSNR'        '%f\t'     2 };%32


channelIx = 4;
harmIx = 9;
freqIx = 5;

fid=fopen(fname);
if fid == -1
    error('File %s could not be opened.',fname);
end

tline=fgetl(fid);
dati=textscan(fid, [hdrFields{:,2}], 'delimiter', '\t', 'EmptyValue', nan);
% Convert the chanToSave string into digit only
for i=1:size(dati{1,4})
    chan{1,i}=sscanf(dati{1, 4}{i}, 'hc%d');
end

dati{1,channelIx}=chan';
usCols=[3 4 5 11 13 14 15 16 17 18];


% Fill in essential matrix
for s=1:length(usCols)
    o=usCols(s);
    if o ~= channelIx
        dataMatrix(:,s)=(dati{1, o}(:));
    else
        if isempty(cell2mat((dati{1, o}(:))))
            colHdr = {};
            freqsAnalyzed ={};
            binLevels= nan;
            dataMatrix = nan;
            fprintf('ERROR! rawdata is empty..\n')
            return;
        else
            dataMatrix(:,s)=cell2mat((dati{1, o}(:)));
        end
    end
    
end

binIndices = unique(dataMatrix(:, 4)); % this will always include 0
binLevels=(dati{1, 12}(2:length(binIndices)));
binLevels=binLevels';

dataMatrix=dataMatrix(ismember(int16(dataMatrix(:, 3)),int16(chanToSave(:))),:); % Restricts the running matrix to the selected electrodes

dataMatrix=dataMatrix(dataMatrix(:, 4)>0, :); % Selects all bins but the 0th one (i.e. the average bin)

dataMatrix=dataMatrix(dataMatrix(:, 1)>0, :); % Selects all trials but the 0th one (i.e. the average trial)

[freqsAnalyzed,tmpIx]=unique(dati{1, harmIx}(:));
freqNum = nan(1,length(freqsAnalyzed));
for f = 1:length(freqsAnalyzed)
    freqNum(f) = dati{1,freqIx}(tmpIx(f));
end
[tmp,tmpIx] = sort(freqNum);
freqsAnalyzed = freqsAnalyzed(tmpIx);

for m = 1:length(usCols)
    colHdr{m} = hdrFields{usCols(m),1};
end
colHdr{end+1} = 'ampl';
end



