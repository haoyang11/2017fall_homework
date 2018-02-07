
function [ w,e ] = LMS(d,M,mu,sigM)
% NLMS�㷨
% d��������Ӧ��y��ʵ����Ӧmu�Ƕ�Ӧ�ĸ�������,SigM�ǹ������������
w=zeros(M,1);
e = zeros(length(d),1); 
for i = 1:length(d) 
    e(i,1) = d(i,1) - sigM(i,:)*w;
    mu2 = mu/(norm(sigM(i,:),2)^2);
    w = w+mu2*sigM(i,:)'*e(i,1); 
end 
end
