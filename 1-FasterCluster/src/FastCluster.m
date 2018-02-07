clear all;
close all; 
clc;
%% ��������
rawdata = load('example_distances.dat');
% rawdata = load('matdata\Aggregation0.8.txt.mat');
label=length(rawdata)
if label~=1
    data=rawdata;
else
    data=rawdata.data;
end
%��������
% total_n=length(data);
% sample=rand(total_n,1);
% ratio=0.90;
% tidx=find(sample<ratio);
% data=data(tidx,:);
% 

%% ��������
np = max(max(data)); % �����нϴ�ģ���Ϊ�ܵ���
[N,~] = size(data); % ��������С

%% ѡ��dc
percent=3.0;
fprintf('average percentage of neighbours (hard coded): %5.6f\n', percent);
position=round(N*percent/100);
sorted_data = sort(data(:,3));% ��������
dc=sorted_data (position);  % ȷ��dc

%% �����������
dist = zeros(np,np); % dist matrix
dist(data(:,1)+(data(:,2)-1)*np) = data(:,3); % �����������
dist=dist+dist'; % �������Գƣ�����һ��


%% �ܶȺ;����������
rho = zeros(np,1);% �洢�ܶ�
delta = zeros(np,1);% �洢����
neigh = zeros(np,1);% �洢�ܶȽϴ����Ͻ��ھ�

%% ��˹kernel
dist_temp=exp(-(dist./dc).^2); % �ֱ����˹ָ��
rho=sum(dist_temp,2)-1; % �ܶȼӺͣ���1��Ϊ��ȥ���Խǣ��Լ����Լ������룩
%% cut-off kernel
% dist_temp=(dist<dc); % �ض�
% rho=sum(dist_temp,2); % 

%% �������
max_dist = max(data(:,3));
[~,idx] = sort(rho, 'descend');  % ���ܶ�����׼�������

for i = 2:np
    delta(idx(i)) = max_dist;
    [min_dist,nidx]=min(dist(idx(i),idx(1:i-1)));  % ѡ��������ǰ�����Сֵ
    if min_dist < delta(idx(i)) % ��Ĭ��ֵ�Ƚϲ�����ֵ
        delta(idx(i)) = min_dist; 
        neigh(idx(i)) = idx(nidx);  
    end
end
delta(idx(1)) = max(delta); % conventional choice for point with highest density
%% ���˻�ͼ
trho=rho/max(rho);
tdelta=delta/max(delta);
tproduct=trho.*tdelta;
tsort=sort(tproduct,'descend');
figure;
plot(tsort(1:100), 'o','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k');
title('product graph');
xlabel('\gamma'); ylabel('n');

%% ��������ͼ
disp('decision graph generated');
figure; subplot(2,2,1);
plot(rho,delta, 'o','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k');
title('decision graph');
xlabel('\rho'); ylabel('\delta');






%%  ѡ���������

% ѡ����ֵ
rect = getrect
rho_thres=rect(1);
delta_thres=rect(2);
% rho_thres = 0.2 * max(rho); % coefficient may be revised
% delta_thres = 0.5 * max(delta);

% �Ҿ�������
cluster = zeros(np,1);% �洢���
nc = 0;
ridx=(rho > rho_thres);% �����ܶȹ���
didx=(delta > delta_thres);%������빻��
sidx=ridx+didx;
ncl=find(sidx==2);% ���߶��㹻��
nc=length(ncl);% ���ĸ���
cluster(ncl)=(1:nc);%���
cluster_idx=ncl;%�洢����

fprintf('NUMBER OF CLUSTERS: %i \n', nc);



%% ���ಢ�ж�����core����halo
for i = 1:np
    if cluster(idx(i)) == 0 % ���࣬�ܶȽϴ�Ŀ�ʼ
        cluster(idx(i)) = cluster(neigh(idx(i)));
    end
end
rho_b = zeros(nc,1); % ����ƽ���ֲ��ܶ��Ͻ磬ȡƽ����Ȼ��ȡ����Ǹ�
for i = 1:np-1
    for j = i+1:np
        if cluster(i)~=cluster(j) && dist(i,j)<=dc % find the border region
            local_rho = (rho(i)+rho(j)) / 2;
            if local_rho > rho_b(cluster(i))
                rho_b(cluster(i)) = local_rho; % update the highest density
            end
            if local_rho > rho_b(cluster(j))
                rho_b(cluster(j)) = local_rho; % update the highest density
            end
        end
    end
end
halo = cluster; % �ж��ǲ���halo
for i = 1:np
    if rho(i) < rho_b(cluster(i))
        halo(i) = 0; % update halo data
    end
end

%% ��ͼ
subplot(2,2,2); % show the cluster centers found
plot(rho,delta, 'o','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k');
hold on; title('cluster centers found');
xlabel('\rho'); ylabel('\delta');
cmap = colormap;
for i = 1:nc
    color = cmap(int8(64*i/nc), :);
    plot(rho(cluster_idx(i)),delta(cluster_idx(i)),'o','MarkerSize',3,'MarkerFaceColor',color,'MarkerEdgeColor',color);
end

subplot(2,2,3); % show the result without halo
points = mdscale(dist, 2, 'criterion','metricstress');
hold on; title('clustering result (without halo)');
xlabel('x'); ylabel('y');
for i = 1:nc
    color = cmap(int8(64*i/nc), :);
    plot(points(cluster==i,1),points(cluster==i,2), 'o','MarkerSize',1.5,'MarkerFaceColor',color,'MarkerEdgeColor',color);
end

subplot(2,2,4); % show the result with halo
plot(points(:,1),points(:,2), 'o','MarkerSize',1.5,'MarkerFaceColor','k','MarkerEdgeColor','k');
hold on; title('clustering result (with halo)');
xlabel('x'); ylabel('y');
for i = 1:nc
    color = cmap(int8(64*i/nc), :);
    plot(points(halo==i,1),points(halo==i,2), 'o','MarkerSize',1.5,'MarkerFaceColor',color,'MarkerEdgeColor',color);
end
