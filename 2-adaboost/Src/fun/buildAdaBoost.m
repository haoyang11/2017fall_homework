function abClassifier = buildAdaBoost(trnX, trnY, iter, tstX, tstY)
if nargin < 4
    tstX = [];
    tstY = [];
end
abClassifier = initAdaBoost(iter);

N = size(trnX, 1); % 训练样本个数
sampleWeight = repmat(1/N, N, 1);

for i = 1:iter
    weakClassifier = buildStump(trnX, trnY, sampleWeight);% 选择一个维度
    abClassifier.WeakClas{i} = weakClassifier;% 第一个分类器
    abClassifier.nWC = i;
    % 分类器权值
    abClassifier.Weight(i) = 0.5*log((1-weakClassifier.error)/weakClassifier.error); % 分类器权值
%     abClassifier.Weight(i)=0.1;
    % 样本权值
    label = predStump(trnX, weakClassifier);% 分类器进行预测
    tmpSampleWeight = -1*abClassifier.Weight(i)*(trnY.*label); % N x 1 %得到alpha
    tmpSampleWeight = sampleWeight.*exp(tmpSampleWeight); % N x 1
    sampleWeight = tmpSampleWeight./sum(tmpSampleWeight); % Normalized
%     sampleWeight = tmpSampleWeight; % Normalized
    
    % Predict on training data
    [ttt, abClassifier.trnErr(i)] = predAdaBoost(abClassifier, trnX, trnY);%得到训练误差
    % Predict on test data
    if ~isempty(tstY)
        abClassifier.hasTestData = true;
        [ttt, abClassifier.tstErr(i)] = predAdaBoost(abClassifier, tstX, tstY); %得到测试误差
    end
    % fprintf('\tIteration %d, Training error %f\n', i, abClassifier.trnErr(i));
end
end
