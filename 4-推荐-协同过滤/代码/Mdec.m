close all
clear all
clc

load testMitrx.mat;
load trainMitrx.mat


lr=0.0001;%�趨����ֵ
weight=0.01;%�趨�˵�ֵ
k=100;%�趨k��ֵ
 U=0.01*rand(10000,k);%��ʼ��U���� 
 V=0.01*rand(10000,k);%��ʼ��V����

 thresld=10E-6;
 pos=1;
 loss=0;
 RM=0;
 t1=clock;
while( pos>0 )
    diff=trainA.*(U*V'-trainM);
    dU=diff*V+2*weight*U;%J��U��ƫ��
    dV=diff'*U+2*weight*V;%J��V��ƫ��
    U=U-lr*dU;%����U
    V=V-lr*dV;%����V
    [J,R]= floss(U,V,weight,trainA,trainM,testA,testM)%����loss�Լ�RMSE
    loss=[loss,J]; % �洢loss
    RM=[RM,R];
    if(pos>100) %����ֹͣ����
        break;
    end
%     if( abs(J-loss(pos))<thresld)%ѭ���������ж�������R��ֵ������ֵ
%         break;
%     end
    pos=pos+1;
    t3=clock;
end
t2=clock;
disp(['����ʱ�䣺',num2str(etime(t2,t1)),'s'])
%%  ��ͼ
loss=loss(2:end);
RM=RM(2:end);
subplot(2,1,1)
plot(loss)
title('�Ż�����ֵ�仯���')
subplot(2,1,2)
plot(RM)
title('RMSE�仯���')


%% ����loss

function [J,R]= floss(U,V,weight,trainA,trainM,testA,testM)
    XX=U*V';%�������ƾ���
    diff=trainA.*(XX-trainM);% ���ƾ�����������ֵ֮�����  
    J=sum(sum(diff.*diff))*0.5+sum(sum(U.*U))*weight+sum(sum(V.*V))*weight;%����Ŀ�꺯��J��ֵ
    testX=testA.*XX;
    R=RMSE(testX,testM);
end



%% ����RMSE
function [rmse]= RMSE(A,B)
diff=A-B;
tem=double(diff~=0);
n=sum(sum(tem));
rmse=sqrt(sum(sum(diff.*diff))/n);
end



