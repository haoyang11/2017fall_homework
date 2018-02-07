close all;
clear all;
clc;

%% 实现基于用户的协同过滤
load testMitrx.mat;
load trainMitrx.mat

t1=clock;
%计算相似度矩阵
norm=sqrt(sum(trainM.^2,2)); % 求训练集范数
normM=norm*norm'; %  范数矩阵
sim=(trainM*trainM')./normM;% x*y/|x||y|
% 预测分数，计算RMSE
% normS=repmat(sum(sim,2),1,10000);

score=((sim*trainM)./(sim*trainA)).*testA; % 预测并做归一化,选择有分数的测试点
% normS=repmat(sum(sim,2),1,10000);
% score=((sim*trainM)./normS).*testA; % 预测并做归一化,选择有分数的测试点
rmse=RMSE(score,testM);
t2=clock;
disp(['运行时间：',num2str(etime(t2,t1)),'s'])

%% 计算RMSE
function [rmse]= RMSE(A,B)
diff=A-B;
tem=double(diff~=0);
n=sum(sum(tem));
rmse=sqrt(sum(sum(diff.*diff))/n);
end

