close all;
clear all;
clc;

%% 实现基于用户的协同过滤
load testMitrx.mat;
load trainMitrx.mat
load testT.mat;
load trainT.mat

t1=clock;
%计算相似度矩阵
norm=sqrt(sum(trainM.^2,2)); % 求训练集范数
normM=norm*norm'; %  范数矩阵
sim=(trainM*trainM')./normM;% x*y/|x||y|
%增加修正项
bias=zeros(10000,10000);
for i=1:10000
    buffer=repmat(trainT(i,:),10000,1);
    dif=abs(trainT-buffer);%求解时间差
    dif(dif>10)=0;
    A=double(dif~=0);
    dif=A.*(1-dif/10)/1000;% 1-时间差/10,由于矩阵稀疏，考虑a=1/1000，使总和不超过1
    cor=sum(dif,2);
    cor(i)=0;
    bias(i,:)=cor';
    bias(:,i)=cor;
end
bias=bias+eye(10000,10000);
%进行相似性修正
sim=sim+bias;

% 预测分数，计算RMSErmse
score=((sim*trainM)./(sim*trainA)).*testA; % 预测并做归一化,选择有分数的测试点
rmse=RMSE(score,testM);
t2=clock;
disp(['运行时间：',num2str(etime(t2,t1)),'s'])

%% 计算RMSE
function [rmse]= RMSE(A,B)
diff=A-B;
tem=double(diff~=0);
n=sum(sum(tem))
rmse=sqrt(sum(sum(diff.*diff))/n);
end

