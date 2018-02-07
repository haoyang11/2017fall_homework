clear all;
close all; 
clc;
%% 载入数据
rawdata = load('example_distances.dat');
% rawdata = load('matdata\Aggregation0.8.txt.mat');
label=length(rawdata)
if label~=1
    data=rawdata;
else
    data=rawdata.data;
end
%采样数据
% total_n=length(data);
% sample=rand(total_n,1);
% ratio=0.90;
% tidx=find(sample<ratio);
% data=data(tidx,:);
% 

%% 基本参数
np = max(max(data)); % 两列中较大的，作为总点数
[N,~] = size(data); % 距离矩阵大小

%% 选择dc
percent=3.0;
fprintf('average percentage of neighbours (hard coded): %5.6f\n', percent);
position=round(N*percent/100);
sorted_data = sort(data(:,3));% 距离排序
dc=sorted_data (position);  % 确定dc

%% 构建距离矩阵
dist = zeros(np,np); % dist matrix
dist(data(:,1)+(data(:,2)-1)*np) = data(:,3); % 构建距离矩阵
dist=dist+dist'; % 距离矩阵对称，加另一半


%% 密度和距离基本参数
rho = zeros(np,1);% 存储密度
delta = zeros(np,1);% 存储距离
neigh = zeros(np,1);% 存储密度较大距离较近邻居

%% 高斯kernel
dist_temp=exp(-(dist./dc).^2); % 分别求高斯指数
rho=sum(dist_temp,2)-1; % 密度加和，减1是为了去掉对角（自己和自己做距离）
%% cut-off kernel
% dist_temp=(dist<dc); % 截断
% rho=sum(dist_temp,2); % 

%% 计算距离
max_dist = max(data(:,3));
[~,idx] = sort(rho, 'descend');  % 对密度排序，准备求距离

for i = 2:np
    delta(idx(i)) = max_dist;
    [min_dist,nidx]=min(dist(idx(i),idx(1:i-1)));  % 选出排在最前面的最小值
    if min_dist < delta(idx(i)) % 和默认值比较并更新值
        delta(idx(i)) = min_dist; 
        neigh(idx(i)) = idx(nidx);  
    end
end
delta(idx(1)) = max(delta); % conventional choice for point with highest density
%% 做乘积图
trho=rho/max(rho);
tdelta=delta/max(delta);
tproduct=trho.*tdelta;
tsort=sort(tproduct,'descend');
figure;
plot(tsort(1:100), 'o','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k');
title('product graph');
xlabel('\gamma'); ylabel('n');

%% 产生决策图
disp('decision graph generated');
figure; subplot(2,2,1);
plot(rho,delta, 'o','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k');
title('decision graph');
xlabel('\rho'); ylabel('\delta');






%%  选择距离中心

% 选择阈值
rect = getrect
rho_thres=rect(1);
delta_thres=rect(2);
% rho_thres = 0.2 * max(rho); % coefficient may be revised
% delta_thres = 0.5 * max(delta);

% 找聚类中心
cluster = zeros(np,1);% 存储类别
nc = 0;
ridx=(rho > rho_thres);% 满足密度够大
didx=(delta > delta_thres);%满足距离够大
sidx=ridx+didx;
ncl=find(sidx==2);% 两者都足够大
nc=length(ncl);% 中心个数
cluster(ncl)=(1:nc);%标号
cluster_idx=ncl;%存储中心

fprintf('NUMBER OF CLUSTERS: %i \n', nc);



%% 分类并判断属于core还是halo
for i = 1:np
    if cluster(idx(i)) == 0 % 归类，密度较大的开始
        cluster(idx(i)) = cluster(neigh(idx(i)));
    end
end
rho_b = zeros(nc,1); % 估计平均局部密度上界，取平均，然后取最大那个
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
halo = cluster; % 判断是不是halo
for i = 1:np
    if rho(i) < rho_b(cluster(i))
        halo(i) = 0; % update halo data
    end
end

%% 绘图
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
