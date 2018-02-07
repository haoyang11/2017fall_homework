function [ w,e ] = RLS(d,y,sigM,M,lambda )
%RLS �˴���ʾ�йش˺�����ժҪ
%��ʼ����
w=zeros(M,1);
delta = 0.01*norm(y,2)^2;
P = delta^-1 * eye ( M, M );
e = zeros(length(d),1);
%�������²��������
    for i = 1:length(d)
        k = P*sigM(i,:)'/(lambda + sigM(i,:)*P*sigM(i,:)');
        e(i,1) = d(i,1) -sigM(i,:)*w;
        w = w + k*e(i,1);
        P = (P - k*sigM(i,:)*P)/lambda;
    end
end

