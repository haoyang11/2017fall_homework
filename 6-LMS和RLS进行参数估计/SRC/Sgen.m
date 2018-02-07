function [SigM,Sy,Sd] = Sgen(N,SNR,K,delay)
%SGEN 此处显示有关此函数的摘要
%   此处显示详细说明
sdata=rand(N,1);
sratio=0.5;
sig=(sdata>sratio)*2-1;
% 信道输出增加噪声
h=[0.3,0.9,0.3];
% 实际输出
Sy = awgn(filter(h, 1, sig), SNR);

% 期望输出
Sd=[zeros(delay, 1);sig(1:N-delay)];
Sy_add=[zeros(K-1, 1);Sy];


SigM = [];
for i = 1:K 
    SigM = [SigM,Sy_add(K+1-i:end+1-i,1)]; 
end

end

