close all;
clear all;
clc;

addpath 'fun';
% load('FHR_data.mat')
load('data_wineq.mat')

nfold = 5;
iter = 30;
tstError = zeros(nfold, iter);
trnError = zeros(nfold, iter);
[trnM, tstM] = buildCVMatrix(size(X, 1), nfold);% 构建交叉验证矩阵
for n = 1:nfold
    fprintf('\tFold %d\n', n);
    idx_trn = logical(trnM(:, n) == 1);
    trnX = X(idx_trn, :);
    tstX = X(~idx_trn, :);
    trnY = Y(idx_trn);
    tstY = Y(~idx_trn);
    abClassifier = buildAdaBoost(trnX, trnY, iter, tstX, tstY);
    trnError(n, :) = abClassifier.trnErr;
    tstError(n, :) = abClassifier.tstErr;
end

plot(1:iter, mean(trnError, 1));
hold on;
plot(1:iter, mean(tstError, 1));
legend('平均训练误差','平均测试误差') 

