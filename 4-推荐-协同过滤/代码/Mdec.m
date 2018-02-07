close all
clear all
clc

load testMitrx.mat;
load trainMitrx.mat


lr=0.0001;%设定α的值
weight=0.01;%设定λ的值
k=100;%设定k的值
 U=0.01*rand(10000,k);%初始化U矩阵 
 V=0.01*rand(10000,k);%初始化V矩阵

 thresld=10E-6;
 pos=1;
 loss=0;
 RM=0;
 t1=clock;
while( pos>0 )
    diff=trainA.*(U*V'-trainM);
    dU=diff*V+2*weight*U;%J对U的偏导
    dV=diff'*U+2*weight*V;%J对V的偏导
    U=U-lr*dU;%更新U
    V=V-lr*dV;%更新V
    [J,R]= floss(U,V,weight,trainA,trainM,testA,testM)%计算loss以及RMSE
    loss=[loss,J]; % 存储loss
    RM=[RM,R];
    if(pos>100) %迭代停止条件
        break;
    end
%     if( abs(J-loss(pos))<thresld)%循环结束的判定条件，R的值低于阈值
%         break;
%     end
    pos=pos+1;
    t3=clock;
end
t2=clock;
disp(['运行时间：',num2str(etime(t2,t1)),'s'])
%%  作图
loss=loss(2:end);
RM=RM(2:end);
subplot(2,1,1)
plot(loss)
title('优化函数值变化情况')
subplot(2,1,2)
plot(RM)
title('RMSE变化情况')


%% 计算loss

function [J,R]= floss(U,V,weight,trainA,trainM,testA,testM)
    XX=U*V';%产生估计矩阵
    diff=trainA.*(XX-trainM);% 估计矩阵结果和已有值之间的误差。  
    J=sum(sum(diff.*diff))*0.5+sum(sum(U.*U))*weight+sum(sum(V.*V))*weight;%计算目标函数J的值
    testX=testA.*XX;
    R=RMSE(testX,testM);
end



%% 计算RMSE
function [rmse]= RMSE(A,B)
diff=A-B;
tem=double(diff~=0);
n=sum(sum(tem));
rmse=sqrt(sum(sum(diff.*diff))/n);
end



