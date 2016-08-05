%Analysis of Proliferative samples July 2016
clear; clc; close all

%% Input
%Must modify this section for each assay:
filename='D:\OneDrive\Research\0C Triculture Christi\Exp7\ProlifPrimary_CellLines_Compare_revised.xlsx';
standardLot = 64020782;
standardDilutions = 1/3;

%% Import Data
beadCounts=readtable(filename,'Sheet','BeadCounts','ReadRowNames',false);
sampleMFI = readtable(filename,'Sheet','SampleMFI','ReadRowNames',false);
sampleMFI{1:end,2:end}(beadCounts{1:end,2:end}<50)=NaN; %Remove MFIs with low bead counts

standardLabels = {'S0','S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14','S15','S16'}; %Possible sample names for standards
blankLabels = {'Blank','blank','BLANK'};
standardMFI = sampleMFI(ismember(sampleMFI{:,1},standardLabels),:);
[~, idx] = sort(standardMFI{:,1});
standardMFI = standardMFI(idx,:);
blankMFI = sampleMFI(ismember(sampleMFI{:,1},blankLabels),:);
sampleMFIsorted = sampleMFI(~ismember(sampleMFI{:,1},[standardLabels,blankLabels]),:);
[~, idx] = sort(sampleMFIsorted{:,1});
sampleMFIsorted = sampleMFIsorted(idx,:);
[sorted, sortIdx] = sort(sampleMFIsorted.Properties.VariableNames(2:end));
sampleMFIsorted = sampleMFIsorted(:,[1,sortIdx+1]);

%Check for outliers in standard MFIs
%*********************************************************************************
%To Do
%Grubbs test a nice balance of easy+appropriate, but supposed to be for >6
%replicates; see also http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-123
%**********************************************************************************

%Read in standard concentrations
standardTable = readtable('D:\OneDrive\Research\Protocols\Luminex\LuminexStandards.xlsx');
standardConc=standardTable(standardTable{:,1}==standardLot,:);
[sorted, sortIdx] = sort(standardConc.Properties.VariableNames(2:end));
standardConcSorted = standardConc(:,[1,sortIdx+1]);
%Calculate concentrations of standards
currDilution = 1; standards = unique(standardMFI{:,1});
dilutionFactor=table(standardMFI{:,1},ones(length(standardMFI{:,1}),1));
for n = 1:length(standards)
    dilutionFactor{strcmp(dilutionFactor{:,1},standards(n)),2}=currDilution;
    currDilution = currDilution*standardDilutions;
end
allStandardsConc = dilutionFactor{:,2}*standardConcSorted{1,2:end};
allStandardsConc = [array2table(dilutionFactor{:,1}),array2table(allStandardsConc)];
allStandardsConc.Properties.VariableNames = ['Standard',standardTable.Properties.VariableNames(2:end)];
allStandardsConc.Name = dilutionFactor(:,1);
%% Method 1: Fit 5-parameter logistic curves to standards, omitting blanks
sampleConc_1 = sampleMFIsorted; sampleConc_1{1:end,2:end}=0;
A=zeros(1,size(sampleConc_1,2)-1); B=A; C=A; D=A; E=A; %initialize
flag=zeros(size(sampleConc_1,1),size(sampleConc_1,2)-1);
sampleMFIadj=sampleMFIsorted;
for c=1:size(sampleConc_1,2)-1
    [cf{c} gf{c}] = L5P(allStandardsConc{1:end,c+1},standardMFI{1:end,c+1}); 
    A(c)=cf{c}.A; B(c)=cf{c}.B; C(c)=cf{c}.C; D(c)=cf{c}.D; E(c)=cf{c}.E;
    flag(sampleMFIsorted{:,c+1}<A(c),c)=-1; sampleMFIadj{sampleMFIsorted{:,c+1}<A(c),c+1}=A(c); 
    flag(sampleMFIsorted{:,c+1}>D(c),c)=1; sampleMFIadj{sampleMFIsorted{:,c+1}>D(c),c+1}=0.99*D(c); 
    sampleConc_1{:,c+1} = L5Pinv(cf{c},sampleMFIadj{:,c+1});
    rsquares(c) = gf{c}.rsquare;
end
%Average replicates - assumes that sample replicates have identical names,
%omit NaNs
% **********************************************************************
% To Do: check for bad wells (e.g. no sample added? all MFIs lower than
% blank?)
% **********************************************************************
uniqueSampleNames = unique(sampleMFIsorted{:,1});
meanSampleConc_1 = zeros(length(uniqueSampleNames),size(sampleConc_1,2)-1);
for sampleIdx=1:length(uniqueSampleNames)
    meanSampleConc_1(sampleIdx,:) = mean(sampleMFIsorted{strcmp(sampleMFIsorted{:,1},uniqueSampleNames(sampleIdx)),2:end},1,'omitnan');
end
meanSampleConcTab_1 = sampleConc_1(1:length(uniqueSampleNames),:);
meanSampleConcTab_1{:,1} = uniqueSampleNames;
meanSampleConcTab_1{:,2:end} = meanSampleConc_1;
