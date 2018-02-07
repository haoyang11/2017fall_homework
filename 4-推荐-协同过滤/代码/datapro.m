close all;
clear all;
clc;

%构建用户ID数据表
username = load('users.txt');
[m,n]=size(username);
for i =1:m
    users(username(i)) = i;
end
save users.mat users;
clear username;
%构建训练集行为矩阵
tic
train = load('netflix_train.txt');
train = train(:,1:3);
[m,n]=size(train);
for i =1:m
    trainM(users(train(i,1)),train(i,2))=train(i,3);
end
trainA=double(trainM~=0);
toc
save trainMitrx.mat trainM trainA;
clear trainM trainA
t1=clock;
%构建测试集行为矩阵
test = load('netflix_test.txt');
test = test(:,1:3);
[m,n]=size(test);
for i= 1:m
    testM(users(test(i,1)),test(i,2))=test(i,3);
end
testA=double(testM~=0);
t2=clock;
etime(t2,t1);
save testMitrx.mat testM testA;
clear all;

[movieN,~, movieM] = xlsread('movie.xlsx');
save movieM.mat movieN movieM
clear all;