function [ w,e ] = RLS(d,y,sigM,M,lambda )
%RLS 此处显示有关此函数的摘要
%初始条件
w=zeros(M,1);
delta = 0.01*norm(y,2)^2;
P = delta^-1 * eye ( M, M );
e = zeros(length(d),1);
%迭代更新并给出结果
    for i = 1:length(d)
        k = P*sigM(i,:)'/(lambda + sigM(i,:)*P*sigM(i,:)');
        e(i,1) = d(i,1) -sigM(i,:)*w;
        w = w + k*e(i,1);
        P = (P - k*sigM(i,:)*P)/lambda;
    end
end

