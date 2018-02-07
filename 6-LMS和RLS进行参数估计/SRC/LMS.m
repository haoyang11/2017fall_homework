
function [ w,e ] = LMS(d,M,mu,sigM)
% NLMS算法
% d是期望响应，y是实际响应mu是对应的更新速率,SigM是构建的输入矩阵
w=zeros(M,1);
e = zeros(length(d),1); 
for i = 1:length(d) 
    e(i,1) = d(i,1) - sigM(i,:)*w;
    mu2 = mu/(norm(sigM(i,:),2)^2);
    w = w+mu2*sigM(i,:)'*e(i,1); 
end 
end
