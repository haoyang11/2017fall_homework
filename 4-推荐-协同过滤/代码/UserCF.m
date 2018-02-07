close all;
clear all;
clc;

%% ʵ�ֻ����û���Эͬ����
load testMitrx.mat;
load trainMitrx.mat

t1=clock;
%�������ƶȾ���
norm=sqrt(sum(trainM.^2,2)); % ��ѵ��������
normM=norm*norm'; %  ��������
sim=(trainM*trainM')./normM;% x*y/|x||y|
% Ԥ�����������RMSE
% normS=repmat(sum(sim,2),1,10000);

score=((sim*trainM)./(sim*trainA)).*testA; % Ԥ�Ⲣ����һ��,ѡ���з����Ĳ��Ե�
% normS=repmat(sum(sim,2),1,10000);
% score=((sim*trainM)./normS).*testA; % Ԥ�Ⲣ����һ��,ѡ���з����Ĳ��Ե�
rmse=RMSE(score,testM);
t2=clock;
disp(['����ʱ�䣺',num2str(etime(t2,t1)),'s'])

%% ����RMSE
function [rmse]= RMSE(A,B)
diff=A-B;
tem=double(diff~=0);
n=sum(sum(tem));
rmse=sqrt(sum(sum(diff.*diff))/n);
end

