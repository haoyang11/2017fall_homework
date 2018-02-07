close all;
clear all;
clc;

%% 设置基本参数
N=500;
SNR=25;
K=11;
mu=0.5;
delay=7;
%% 单次实验结果
% 信号产生
[SigM,Sy,Sd] = Sgen(N,SNR,K,delay);
%NLMS算法估计
[wlms,elms] = LMS(Sd,K,mu,SigM);
elms=elms.^2;
%RLS算法估计
lambda=0.8;
[wrls,erls] = RLS(Sd,Sy,SigM,K,lambda );
erls=erls.^2;


%% 多次实验
wlms_mean=zeros(K,1);
wrls_mean=zeros(K,1);
erls_mean=zeros(N,1);
elms_mean=zeros(N,1);
iter=20;



for i=1:iter
    % 信号产生
    [SigM,Sy,Sd] = Sgen(N,SNR,K,delay);
    %NLMS算法估计
    [wlms_iter,elms_iter] = LMS(Sd,K,mu,SigM);
    wlms_mean=wlms_mean+wlms_iter;
    elms_mean=elms_mean+elms_iter.^2;
    %RLS算法估计
    [wrls_iter,erls_iter] = RLS(Sd,Sy,SigM,K,lambda );
    wrls_mean=wrls_mean+wrls_iter;
    erls_mean=erls_mean+erls_iter.^2;
end
elms_mean=elms_mean/iter;
erls_mean=erls_mean/iter;
wlms_mean=wlms_mean/iter;
wrls_mean=wrls_mean/iter;


%% 实验结果显示

figure;
subplot(1,2,1);
cmap = colormap;
color = cmap(int8(64*2/8), :);
plot(elms(delay+1:end),'o','MarkerSize',3,'MarkerFaceColor',color,'MarkerEdgeColor',color,'linewidth',2,'LineStyle','-');
xlabel('迭代次数'); ylabel('LMS误差');
title('单次LMS误差平方曲线')
hold on; 
subplot(1,2,2);
color = cmap(int8(64*3/8), :);
plot(elms_mean(delay+1:end),'o','MarkerSize',3,'MarkerFaceColor',color,'MarkerEdgeColor',color,'linewidth',2,'LineStyle','-');
xlabel('迭代次数'); ylabel('LMS误差');
title('多次LMS平均误差平方曲线')

figure;
subplot(1,2,1);
cmap = colormap;
color = cmap(int8(64*4/8), :);
plot(erls(delay+1:end),'o','MarkerSize',3,'MarkerFaceColor',color,'MarkerEdgeColor',color,'linewidth',2,'LineStyle','-');
xlabel('迭代次数'); ylabel('LMS误差');
title('单次RLS误差平方曲线')
hold on; 
subplot(1,2,2);
color = cmap(int8(64*5/8), :);
plot(erls_mean(delay+1:end),'o','MarkerSize',3,'MarkerFaceColor',color,'MarkerEdgeColor',color,'linewidth',2,'LineStyle','-');
xlabel('迭代次数'); ylabel('LMS误差');
title('多次RLS平均误差平方曲线')




