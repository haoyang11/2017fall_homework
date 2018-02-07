close all;
clear all;
clc;


%%   ����������
N=300;
noise=0;
[data]=data_gen(N,noise);
% ������չʾ
plot(data(:,1),data(:,2),'bx')

%% ���в������Ʋ��ȶԽ��
disp('--------------------ԭʼ����ֵ------------------');
mu1 = [10,3]';%��ѧ����
mu2 = [1,1]';
mu3 = [5,4]';
sigma1=[1,0;0,1];%Э�������
sigma2=[1.5,0;0,1.5];
sigma3=[2,0;0,2];
mu_ori=[mu1,mu2,mu3]
sigma1=[1,0;0,1];%Э�������
sigma2=[3,0;0,3];
sigma3=[2,0;0,2];
sigma_ori=[sigma1,sigma2,sigma3]
w=[1,1,1]./3

disp('--------------------�ֶ���ʼ��------------------');
epson =1e-10;
% �ֶ���ʼ��
mu1 = [3,5]';%��ѧ����
mu2 = [2,1]';
mu3 = [0,3]';
mu=[mu1,mu2,mu3]
sigma1=[1,0;0,1];%Э�������
sigma2=[3,0;0,3];
sigma3=[2,0;0,2];
sigma=[sigma1,sigma2,sigma3]
phi=[0.5,0.4,0.1]
disp('---------------------����ֵ------------------')
[L,mu_1,sigma_1,weight_1]=EM_GMM(data,mu,sigma,phi,epson);

mu_1
sigma_1
weight_1
figure
plot(L,'R*')
title('EM�㷨����GMM����-��Ȼ��������');





disp('---------------------��ȫ�޲��ĳ�ʼֵ����------------------')
mu1 = [0,0]';%��ѧ����
mu2 = [0,0]';
mu3 = [0,0]';
mu=[mu1,mu2,mu3]
sigma1=[1,0;0,1];%Э�������
sigma2=[1,0;0,1];
sigma3=[1,0;0,1];
sigma=[sigma1,sigma2,sigma3]
phi=[1/3,1/3,1/3]


disp('---------------------�޲���ʼֵʱ�Ĺ���ֵ------------------')
[L,mu_1,sigma_1,weight_1]=EM_GMM(data,mu,sigma,phi,epson);

mu_1
sigma_1
weight_1
figure
plot(L,'R')
title('EM�㷨����GMM����-�޲���ʼֵ-��Ȼ��������');



%%  �����ʼ���������ȶԽ��

% ��һ��
[mu,sigma,weight]=param_gen();
[L1,mu_1,sigma_1,weight_1]=EM_GMM(data,mu,sigma,phi,epson);
mu_1
sigma_1
weight_1
% �ڶ���
[mu,sigma,weight]=param_gen()
[L2,mu_2,sigma_2,weight_2]=EM_GMM(data,mu,sigma,phi,epson);
mu_2
sigma_2
weight_2
% ������
[mu,sigma,weight]=param_gen()
[L3,mu_3,sigma_3,weight_3]=EM_GMM(data,mu,sigma,phi,epson);
mu_3
sigma_3
weight_3

% �ԱȲ�ͬ���
figure
plot(L1','*');
hold on
plot(L2','s');
plot(L3','x');

legend('��һ�������ȻֵL1','�ڶ��������ȻֵL2','�����������ȻֵL3');

title('EM�㷨����GMM������Ȼ�����仯');



%%  �����������������������
disp('---------------------�����Թ��Ƶ�Ӱ��------------------')
% ���ó�ʼ����ֵ
mu1 = [3,5]';%��ѧ����
mu2 = [2,1]';
mu3 = [0,3]';
mu=[mu1,mu2,mu3];
sigma1=[1,0;0,1];%Э�������
sigma2=[3,0;0,3];
sigma3=[2,0;0,2];
sigma=[sigma1,sigma2,sigma3];
phi=[0.5,0.4,0.1];
% ԭʼ������
data;
[L1,~,~,~]=EM_GMM(data,mu,sigma,phi,epson);

% ��������һ��
[m,n]=size(data);
noise=0.01;
ndata=noise*rand(m,n);
data=data+ndata;
[L2,~,~,~]=EM_GMM(data,mu,sigma,phi,epson);
% �������ڶ���
 noise=0.1;
ndata=noise*rand(m,n);
data=data+ndata;
[L3,~,~,~]=EM_GMM(data,mu,sigma,phi,epson);

% �������ڶ���
 noise=2;
ndata=noise*rand(m,n);
data=data+ndata;
[L4,~,~,~]=EM_GMM(data,mu,sigma,phi,epson);




figure
plot(L1','*');
hold on
plot(L2','s');
plot(L3','x');
plot(L4','>');

legend('������L1','0.01��������L2','0.1��������L3','2��������L4');

title('�����Թ��ƹ��̵�Ӱ��-��Ȼ��������');



%%  �������ݺ���

function   [data]=data_gen(N,noise)
mu1 = [10,3]';%��ѧ����
mu2 = [1,1]';
mu3 = [5,4]';
sigma1=[1,0;0,1];%Э�������
sigma2=[1.5,0;0,1.5];
sigma3=[2,0;0,2];
w=[1,1,1]./3;

r1=mvnrnd(mu1,sigma1,N);
r2=mvnrnd(mu2,sigma2,N);
r3=mvnrnd(mu3,sigma3,N);
rawdata=[r1,r2,r3];
% ���������
rate = floor(rand(N,1)*3)+1;
data=zeros(N,2);
% ȡ����
    for i=1:3
        idx=find(rate==i);
        data(idx,:)=rawdata(idx,i*2-1:i*2);
    end
% ��������
[m,n]=size(data);
ndata=noise*rand(m,n);
data=data+ndata;

end
%%  �����ʼ������
function   [mu,sigma,weight]=param_gen()
% ������ֵ
mu=floor(rand(2,3)*10);
%Э�������
s1=ceil(rand(1)*10)*eye(2);
s2=ceil(rand(1)*10)*eye(2);
s3=ceil(rand(1)*10)*eye(2);
sigma=[s1,s2,s3];
% Ȩֵ
w_i=rand(1,3);
weight=w_i/sum(w_i);
end


%%  ʹ��EM�㷨��GMM�������й���
function   [L,mu_s,sigma_s,weight_s]=EM_GMM(data,mu,sigma,phi,epson)
N=length(data);
T=5000;
w = zeros(N,3);
muarr=[];
sigmarr=[];
phiarr=[];
error=100;
L(1)=1;
pos=2;

% while(true)
for y=1:1000
    % Expectation
    for k = 1 : 3
        w(:,k) = phi(k)*mvnpdf(data,mu(:,k)',sigma(:,k*2-1:k*2));%  ����ÿһ���������Բ���w���й���
    end
     % �м䴩�������Ȼ����
    L(pos)=sum(log(sum(w,2)))/N;
    err=L(pos)-L(pos-1);
    if(abs(err)<epson)
        break;
    end
    pos=pos+1;
    
    w = w./repmat(sum(w,2),[1 3]);

    % Maximization
    for k = 1 : 3
        mu(:,k) = (w(:,k)'*data / sum(w(:,k)))';
%         sigma(:,k*2-1:k*2) = w(:,k)'*((data-mu(:,k)')'*(data-mu(:,k)')) / sum(w(:,k));
%         temp= kron(w(:,k),((data-mu(:,k)')'*(data-mu(:,k)'))) / sum(w(:,k));% �Է������
        % ��һ��������͹���
%         sigma(:,k*2-1:k*2)= sqrt(reshape(sum(reshape(temp',4,N),2),2,2));

%  ֱ�ӹ������
        cov_mat=diag(w(:,k)'*((data-mu(:,k)').*(data-mu(:,k)')))/sum(w(:,k));
        sigma(:,k*2-1:k*2)=cov_mat;
        phi(k) = sum(w(:,k)) / N;
    end   
    muarr = [muarr;reshape(mu,1,6)];
    sigmarr= [sigmarr;reshape(sigma,1,12)];
    phiarr= [phiarr;phi];
    y=y+1;
end

mu_s=mu;
sigma_s=sigma;
weight_s=phi;
L=L(2:end);
end

