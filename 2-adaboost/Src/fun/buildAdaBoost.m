function abClassifier = buildAdaBoost(trnX, trnY, iter, tstX, tstY)
if nargin < 4
    tstX = [];
    tstY = [];
end
abClassifier = initAdaBoost(iter);

N = size(trnX, 1); % ѵ����������
sampleWeight = repmat(1/N, N, 1);

for i = 1:iter
    weakClassifier = buildStump(trnX, trnY, sampleWeight);% ѡ��һ��ά��
    abClassifier.WeakClas{i} = weakClassifier;% ��һ��������
    abClassifier.nWC = i;
    % ������Ȩֵ
    abClassifier.Weight(i) = 0.5*log((1-weakClassifier.error)/weakClassifier.error); % ������Ȩֵ
%     abClassifier.Weight(i)=0.1;
    % ����Ȩֵ
    label = predStump(trnX, weakClassifier);% ����������Ԥ��
    tmpSampleWeight = -1*abClassifier.Weight(i)*(trnY.*label); % N x 1 %�õ�alpha
    tmpSampleWeight = sampleWeight.*exp(tmpSampleWeight); % N x 1
    sampleWeight = tmpSampleWeight./sum(tmpSampleWeight); % Normalized
%     sampleWeight = tmpSampleWeight; % Normalized
    
    % Predict on training data
    [ttt, abClassifier.trnErr(i)] = predAdaBoost(abClassifier, trnX, trnY);%�õ�ѵ�����
    % Predict on test data
    if ~isempty(tstY)
        abClassifier.hasTestData = true;
        [ttt, abClassifier.tstErr(i)] = predAdaBoost(abClassifier, tstX, tstY); %�õ��������
    end
    % fprintf('\tIteration %d, Training error %f\n', i, abClassifier.trnErr(i));
end
end
