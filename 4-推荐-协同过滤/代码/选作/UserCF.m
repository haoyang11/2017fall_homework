close all;
clear all;
clc;

%% ʵ�ֻ����û���Эͬ����
load testMitrx.mat;
load trainMitrx.mat
load testT.mat;
load trainT.mat

t1=clock;
%�������ƶȾ���
norm=sqrt(sum(trainM.^2,2)); % ��ѵ��������
normM=norm*norm'; %  ��������
sim=(trainM*trainM')./normM;% x*y/|x||y|
%����������
bias=zeros(10000,10000);
for i=1:10000
    buffer=repmat(trainT(i,:),10000,1);
    dif=abs(trainT-buffer);%���ʱ���
    dif(dif>10)=0;
    A=double(dif~=0);
    dif=A.*(1-dif/10)/1000;% 1-ʱ���/10,���ھ���ϡ�裬����a=1/1000��ʹ�ܺͲ�����1
    cor=sum(dif,2);
    cor(i)=0;
    bias(i,:)=cor';
    bias(:,i)=cor;
end
bias=bias+eye(10000,10000);
%��������������
sim=sim+bias;

% Ԥ�����������RMSErmse
score=((sim*trainM)./(sim*trainA)).*testA; % Ԥ�Ⲣ����һ��,ѡ���з����Ĳ��Ե�
rmse=RMSE(score,testM);
t2=clock;
disp(['����ʱ�䣺',num2str(etime(t2,t1)),'s'])

%% ����RMSE
function [rmse]= RMSE(A,B)
diff=A-B;
tem=double(diff~=0);
n=sum(sum(tem))
rmse=sqrt(sum(sum(diff.*diff))/n);
end

