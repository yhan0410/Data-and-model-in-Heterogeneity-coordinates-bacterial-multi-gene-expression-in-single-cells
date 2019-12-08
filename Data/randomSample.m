sampleSizeSet = [10, 20, 50, 100, 200, 500, 1000];
numSampleSize = size(sampleSizeSet, 2);
numGroups = 7;
numSamplingReplicates = 100;
temp = 14;
PCCOverall = zeros(numSampleSize, numSamplingReplicates);
for sampleSizeIndex = 1:numSampleSize
    sampleSize = sampleSizeSet(sampleSizeIndex);
    for replicateIndex = 1:numSamplingReplicates
        CARmCherry = zeros(1, numGroups * sampleSize);
        GFP = zeros(1, numGroups * sampleSize);
        for groupIndex = 1:numGroups
            numCell = size(Fluorescence{groupIndex+temp},1);
            randomSampleIndex = randsample(numCell, sampleSize);
            startIndex = (groupIndex - 1) * sampleSize + 1;
            endIndex = groupIndex * sampleSize;
            CARmCherry(startIndex : endIndex) = Fluorescence{groupIndex+temp}(randomSampleIndex,2);
            GFP(startIndex : endIndex) = Fluorescence{groupIndex+temp}(randomSampleIndex,1);
        end
        PCCOverall(sampleSizeIndex,replicateIndex) = corr(CARmCherry',GFP');
    end
end
meanPCC = mean(PCCOverall,2);
stdPCC =  std(PCCOverall,0,2);
