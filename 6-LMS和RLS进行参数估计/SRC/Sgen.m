function [SigM,Sy,Sd] = Sgen(N,SNR,K,delay)
%SGEN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
sdata=rand(N,1);
sratio=0.5;
sig=(sdata>sratio)*2-1;
% �ŵ������������
h=[0.3,0.9,0.3];
% ʵ�����
Sy = awgn(filter(h, 1, sig), SNR);

% �������
Sd=[zeros(delay, 1);sig(1:N-delay)];
Sy_add=[zeros(K-1, 1);Sy];


SigM = [];
for i = 1:K 
    SigM = [SigM,Sy_add(K+1-i:end+1-i,1)]; 
end

end

