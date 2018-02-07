close all;
clear all;
clc;


list=dir('./data');
mkdir('matdata');% �½�һ���ļ��д洢mat
% ȫ�������mat�ļ�������ŷʽ����
[m,n]=size(list);
cmap=colormap;
for i=3:m     
  filename=[list(i).folder,'\',list(i).name];
  rawdata=load(filename);
  %% ��ÿһ����л�ͼ
  nc=length(unique(rawdata(:,3)));
  str=[list(i).name,'���ݼ�'];
  figure;
  for k = 1:nc
    color = cmap(int8(64*(nc-k+1)/nc), :);
    idx=find(rawdata(:,3)==k);
    hold on
    plot(rawdata(idx,1),rawdata(idx,2), 'o','MarkerSize',1.5,'MarkerFaceColor',color,'MarkerEdgeColor',color);
    title(str);
  end
  
  [L,w,h]=size(rawdata);
  dist=pdist2(rawdata(:,1:2),rawdata(:,1:2));
  len=L*(L-1)/2;
  data=zeros(len,3);
  pos=1;
  for t=1:L
      for j=t+1:L
         data(pos,:)=[t,j,dist(t,j)];
         pos=pos+1;
      end
  end
  %�洢����
 
  datapath=['matdata','\',list(i).name,'.mat'];
  save(datapath,'data');
end

