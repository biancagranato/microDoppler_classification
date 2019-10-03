%% Introduction

% This script was used to import and train alexnet and implement transfer learning.
% Requires: 
%	file/app: export_fig (available at https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig)
%				this function can be replaced 
%	path: expects image categories to be on the 8th level (i.e. a sequence of subfolders with at least
%		8 entries). This can be changed in line 130. It also expects paths to contain backslashes '\',
%		which is the default output of the function on Windows. (i.e. 'C:Users\MATLAB\a\b\c\d\e\f\g')
%	toolboxes: distrib_computing_toolbox, neural_network_toolbox, phased_array_system_toolbox, 
%				signal_blocks, signal_toolbox, statistics_toolbox
%	data: images of different categories organized in folders, whose name is the name of the category. 
%			i.e. images of cats in a folder named cat. Label source can be modified in line 57
%	additional file: at the moment, it also requires a csv file that describes image subclasses (see line 38)
%					the file should be named 'cond.csv' and must contain 4 columns
%
%
% Written and tested on MATLAB 2019b (though it should work on 2017b and above).
%
% Available under MIT license. Copyright (c) 2019 Bianca Granato


%% Make Confmat directory
folders = getFolders; %Get folders with each of the conditions 
mkdir('Confmats') %Make folder where images are saved

%% Run ML code
scores = struct('Accuracy',[],'MCC',[]);
scores(numel(folders)).Accuracy = 0;

net = alexnet;
layersTransfer = net.Layers(1:end-3);
tic
for i = 1:numel(folders)
    tmp = ML(folders{i},layersTransfer);
    scores(i).Accuracy = tmp.acc;
    scores(i).MCC = tmp.mcc;
end
toc
disp('done')

%% Analyze
cond = readtable('cond.csv'); 
%In the future, I will add a field where the different conditions
%are automatically extracted. For now, this will do.
cond.Properties.VariableNames = {'window-func','window-size','no-colors','overlap'};

for i = 1:size(cond,2)
    cats = table2array(cond(:,i));
    figure(1); subplot(2,2,i)
    boxplot([scores.Accuracy],categorical(cats))    
    ylabel('Accuracy'); xlabel(cond.Properties.VariableNames{i})
    set(gca, 'linewidth', 1.25,'TickDir', 'out'); box off
    figure(2); subplot(2,2,i)
    boxplot([scores.MCC],categorical(cond{:,i}))    
    ylabel('MCC'); xlabel(cond.Properties.VariableNames(i))
    set(gca, 'linewidth', 1.25,'TickDir', 'out'); box off
end

function score = ML(f,layersTransfer)
%get name of conditions
slash = find(f == '\');
dirnam = f(slash(end-3)+1:end);
filenam = erase(dirnam,'\');

%make image datastore and get labels from subfolders
imds = imageDatastore(f,'IncludeSubfolders',true, 'LabelSource','foldernames');

%split each label so there are 100 images in the training set and 30 in the validation test
[imdsTrain,imdsValidation,imdsTest] = splitEachLabel(imds,5,10,'randomized');

%modify training set so images match Alexnet's required input size (227x227x3)
augTrain = augmentedImageDatastore([227 227 3],imdsTrain);
augVal = augmentedImageDatastore([227 227 3],imdsValidation);
augTest = augmentedImageDatastore([227 227 3],imdsTest);

%remove last layers of AlexNet

numClasses = numel(categories(imdsTrain.Labels));
layers = [
    layersTransfer
    fullyConnectedLayer(numClasses,'WeightLearnRateFactor',20,'BiasLearnRateFactor',20)
    softmaxLayer
    classificationLayer];

%Set options
options = trainingOptions('sgdm', ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.2, ...
    'LearnRateDropPeriod',5, ...
    'MaxEpochs',4, ...
    'InitialLearnRate',1e-4, ...    
    'ValidationData',augVal, ...
    'ValidationFrequency',3, ...
    'Verbose',false, ...
    'ExecutionEnvironment','gpu');


%Train
netTransfer = trainNetwork(augTrain,layers,options);

%Classify
[testPred,~] = classify(netTransfer,augTest);

%Plot confusion matrix and save
confusionchart(imdsTest.Labels,testPred,...
    'ColumnSummary','column-normalized', ...
    'RowSummary','row-normalized',...
    'Title',filenam)
savename = sprintf('Confmats/Confmat%s.png',filenam);
export_fig(savename) %and this one

%Return accuracy 
testAccuracy = mean(testPred == imdsTest.Labels);

%Get true positive/false positive
numericTest = grp2idx(imdsTest.Labels);
numericPred = grp2idx(testPred);
confIdx = sub2ind([numClasses numClasses], numericTest, numericPred);
mat = zeros(numClasses);

for n = 1:numClasses*numClasses
    mat(n) = sum(confIdx == n);
end

score.mcc = getEmcc(mat);
score.acc = testAccuracy;
end

function f = getFolders
f = dir('**/');
f = {f([f(:).isdir]).folder}';
f(ismember(f,{'.','..'})) = [];
f = unique(f);
f = f(cellfun(@(x) length(find(x == '\'))==7,f));
end

function emcc = getEmcc(c)
% Extended Matthew's Correlation Coefficient
% Gorodkin, J. (2004). doi:10.1016/j.compbiolchem.2004.09.006
N = 0;
ckl = 0;
ckct = 0;
ctcl = 0;
ct = c';
for k = 1:size(c,1)
    for l = 1:size(c,2)
            N = N+c(k,l);
            ckl = ckl + c(k,:)*c(:,l);
            ckct = ckct + c(k,:)*ct(:,l);
            ctcl = ctcl + ct(k,:)*c(:,l);
    end
end

top = N*trace(c) - ckl;
bottom = sqrt(N^2 - ckct)*sqrt(N^2-ctcl);
emcc = top/bottom;
end

%Copyright (c) 2019 Bianca Granato
%
%Permission is hereby granted, free of charge, to any person obtaining a copy
%of this software and associated documentation files (the "Software"), to deal
%in the Software without restriction, including without limitation the rights
%to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%copies of the Software, and to permit persons to whom the Software is
%furnished to do so, subject to the following conditions:
%
%The above copyright notice and this permission notice shall be included in all
%copies or substantial portions of the Software.
%
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%SOFTWARE.