close all;
clear all;
clc;

%�����û�ID���ݱ�
username = load('users.txt');
[m,n]=size(username);
for i =1:m
    users(username(i)) = i;
end
save users.mat users;
clear username;
%����ѵ������Ϊ����
tic
train = load('netflix_train.txt');
train = train(:,[1,2,4]);
[m,n]=size(train);
for i =1:m
    trainT(users(train(i,1)),train(i,2))=train(i,3);
end
trainTA=double(trainT~=0);
toc
save trainT.mat trainT trainTA;
clear trainT trainTA
t1=clock;
%�������Լ���Ϊ����
test = load('netflix_test.txt');
test = test(:,[1,2,4]);
[m,n]=size(test);
for i= 1:m
    testT(users(test(i,1)),test(i,2))=test(i,3);
end
testTA=double(testT~=0);
t2=clock;
etime(t2,t1);
save testT.mat testT testTA;
clear all;