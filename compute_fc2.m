 clear all; close all; clc;
%申明预处理完的数据所在的路径 注意修改
D = 'D:\Scientific research\Ana_date\shuju\A_Homer\A_25Hz 10min oxy\ROI\4ROI\D';
J = 'D:\Scientific research\Ana_date\shuju\A_Homer\A_25Hz 10min oxy\ROI\4ROI\J';
P = 'D:\Scientific research\Ana_date\shuju\A_Homer\A_25Hz 10min oxy\ROI\4ROI\P_10';
%利用dir函数寻找数据路径下面所有以.mat结尾的文件，为后面for循环服务
D_Files = dir(fullfile(D,'*.mat')); %不同组的话Condition可以改为group
J_Files = dir(fullfile(J,'*.mat')); 
P_Files = dir(fullfile(P,'*.mat'));
%获取每个被试具体的文件名存到胞元
FileNamesD = {D_Files.name};
FileNamesJ = {J_Files.name};
FileNamesP = {P_Files.name};

%% compute the measures
fs = 25;  %采样率  注意修改

%计算组D 的功能连接指标，以oxy实例
for sub = 1:length(FileNamesD) %沿着每个被试循环
    %导入第sub个被试的数据
    load(strcat(D, '\', FileNamesD{1,sub}));
    
    % correlation
    % corr_oxy_con1 三维数组 通道*通道*被试
    corr_oxy_D(:,:,sub) = corr(oxy);%汇总所有被试计算得到的功能连接矩阵
  %  r = corr(oxy_down);
    %fisher Z  变换，使其服从正太分布
   % corr_oxy_D_z(:,:,sub) = 0.5*log((1+r)./(1-r));
    
    % coherence 相干
    %for channel1 = 1:size(oxy,2)  %第一个通道  size(oxy,2)这里的2是第2个维度
      %  for channel2 = 1:size(oxy,2) %第二个通道
            %coh_oxy_con1 四维数组 频率*通道*通道*被试  Freq 频率点分布  提取对应频段上的相干指标
            %[coh_oxy_D(:,channel1,channel2,sub),Freq] = mscohere(oxy(:,channel1),oxy(:,channel2),[],[],[],fs);
     %   end
 %   end     
    
    % phase locking value  相位锁值
%    for channel1 = 1:size(oxy,2)  %沿着第一个通道
   %     for channel2 = 1:size(oxy,2) %沿着第二个通道
            %提取1号通道的数据，进行希尔伯特变换，提取瞬时相位
 %^           channel1_oxy_phase = angle(hilbert(oxy(:,channel1)));
            %提取2号通道的数据，进行希尔伯特变换，提取瞬时相位
 %           channel2_oxy_phase = angle(hilbert(oxy(:,channel2)));
            %计算相位差
%            rp_oxy = channel1_oxy_phase - channel2_oxy_phase;
            %计算PLV 通道*通道*被试
%            plv_oxy_D(channel1,channel2,sub) = abs(sum(exp(1i*rp_oxy))/length(rp_oxy));
%        end
 %   end
    %进度条，反映数据处理的进度
    waitbar(sub/length(FileNamesD))
end
%提取所有通道所有被试感兴趣频段上的COH值[0.01 -- 0.1HZ]  coh_oxy_con1  三维数组  通道*通道*被试
%coh_oxy_D = squeeze(mean(coh_oxy_D(Freq > 0.01 & Freq < 0.1,:,:,:),1));
%将条件1的功能连接指标保存到本地磁盘  FCmaxtrix_con1.mat
save FCmaxtrix_D.mat  corr_oxy_D %corr_oxy_D_z coh_oxy_D  plv_oxy_D

%% 

fs = 25;
%计算 组J 的功能连接指标（以oxy作为示例）
for sub = 1:length(FileNamesD) %沿着每个被试循环
    %导入第sub个被试的数据
    load(strcat(D, '\', FileNamesD{1,sub}));
    
    % correlation
    % corr_oxy_con1 三维数组 通道*通道*被试
    corr_roi_D(:,:,sub) = corr(roi);%汇总所有被试计算得到的功能连接矩阵
 %   r = corr(oxy_down);
    %fisher Z  变换，使其服从正太分布
%    
    % coherence 相干
%    for channel1 = 1:size(oxy,2)  %第一个通道  size(oxy,2)这里的2是第2个维度
 %       for channel2 = 1:size(oxy,2) %第二个通道
%          %coh_oxy_con1 四维数组 频率*通道*通道*被试  Freq 频率点分布  提取对应频段上的相干指标
%            [coh_oxy_J(:,channel1,channel2,sub),Freq] = mscohere(oxy(:,channel1),oxy(:,channel2),[],[],[],fs);
 %       end
 %   end     
    
    % phase locking value  相位锁值
 %   for channel1 = 1:size(oxy,2)  %沿着第一个通道
 %       for channel2 = 1:size(oxy,2) %沿着第二个通道
            %提取1号通道的数据，进行希尔伯特变换，提取瞬时相位
  %          channel1_oxy_phase = angle(hilbert(oxy(:,channel1)));
            %提取2号通道的数据，进行希尔伯特变换，提取瞬时相位
  %          channel2_oxy_phase = angle(hilbert(oxy(:,channel2)));
            %计算相位差
 %           rp_oxy = channel1_oxy_phase - channel2_oxy_phase;
            %计算PLV 通道*通道*被试
   %         plv_oxy_J(channel1,channel2,sub) = abs(sum(exp(1i*rp_oxy))/length(rp_oxy));
  %      end
  %  end
    %进度条，反映数据处理的进度
    waitbar(sub/length(FileNamesD))
end
%提取所有通道所有被试感兴趣频段上的COH值[0.01 -- 0.1HZ]  coh_oxy_con1  三维数组  通道*通道*被试
%coh_oxy_J = squeeze(mean(coh_oxy_J(Freq > 0.01 & Freq < 0.1,:,:,:),1));
%将条件1的功能连接指标保存到本地磁盘  FCmaxtrix_con1.mat
save FCmaxtrix_D.mat  corr_roi_D %corr_oxy_J_z coh_oxy_J  plv_oxy_J

%% 

fs = 10;
%计算组P的功能连接指标（以oxy作为示例）
for sub = 1:length(FileNamesP) %沿着每个被试循环
    %导入第sub个被试的数据
    load(strcat(P, '\', FileNamesP{1,sub}));
    
    % correlation
    % corr_oxy_con1 三维数组 通道*通道*被试
    corr_roi_P(:,:,sub) = corr(roi);%汇总所有被试计算得到的功能连接矩阵
 %   r = corr(oxy_down);
    %fisher Z  变换，使其服从正太分布
%    corr_oxy_P_z(:,:,sub) = 0.5*log((1+r)./(1-r));
    
    % coherence 相干
 %   for channel1 = 1:size(oxy,2)  %第一个通道  size(oxy,2)这里的2是第2个维度
 %       for channel2 = 1:size(oxy,2) %第二个通道
            %coh_oxy_con1 四维数组 频率*通道*通道*被试  Freq 频率点分布  提取对应频段上的相干指标
 %           [coh_oxy_P(:,channel1,channel2,sub),Freq] = mscohere(oxy(:,channel1),oxy(:,channel2),[],[],[],fs);
  %      end
  %  end     
    
    % phase locking value  相位锁值
%    for channel1 = 1:size(oxy,2)  %沿着第一个通道
 %       for channel2 = 1:size(oxy,2) %沿着第二个通道
  %          %提取1号通道的数据，进行希尔伯特变换，提取瞬时相位
 %           channel1_oxy_phase = angle(hilbert(oxy(:,channel1)));
            %提取2号通道的数据，进行希尔伯特变换，提取瞬时相位
  %          channel2_oxy_phase = angle(hilbert(oxy(:,channel2)));
            %计算相位差
 %           rp_oxy = channel1_oxy_phase - channel2_oxy_phase;
            %计算PLV 通道*通道*被试
   %         plv_oxy_P(channel1,channel2,sub) = abs(sum(exp(1i*rp_oxy))/length(rp_oxy));
 %       end
%   end
    %进度条，反映数据处理的进度
    waitbar(sub/length(FileNamesP))
end
%提取所有通道所有被试感兴趣频段上的COH值[0.01 -- 0.1HZ]  coh_oxy_con1  三维数组  通道*通道*被试
% coh_oxy_P = squeeze(mean(coh_oxy_P(Freq > 0.01 & Freq < 0.1,:,:,:),1));
%将条件1的功能连接指标保存到本地磁盘  FCmaxtrix_con1.mat
save FCmaxtrix_P_10.mat  corr_roi_P % corr_oxy_P_z coh_oxy_P  plv_oxy_P

%% 绘图组水平

clc;clear all;close all
%导入条件1和条件2的功能连接指标
load('FCmaxtrix_D.mat')
load('FCmaxtrix_J.mat')
load('FCmaxtrix_P.mat')
%对被试维度平均，计算组平均的功能连接矩 阵
%corr  皮尔逊相关系数
corr_oxy_D_avg = mean(corr_oxy_D,3);   corr_oxy_J_avg = mean(corr_oxy_J,3);  corr_oxy_P_avg = mean(corr_oxy_P,3);%通道*通道
%corr fisher-z fisherZ后的皮尔逊相关系数
%corr_oxy_con1_z_avg = mean(corr_oxy_con1_z,3);   corr_oxy_con2_z_avg = mean(corr_oxy_con2_z,3);
%coh 相干
%coh_oxy_con1_avg = mean(coh_oxy_con1,3);   coh_oxy_con2_avg = mean(coh_oxy_con2,3);
%plv 相位锁值
%plv_oxy_con1_avg = mean(plv_oxy_con1,3);   plv_oxy_con2_avg = mean(plv_oxy_con2,3);

%绘制相关的组平均矩阵
figure; %生成空画布
%将画布分割成1行2列，在第一列绘制条件1的相关矩阵 将colorbar范围定为0.1-0.5 标题为con1
imagesc(1:27, 1:27, corr_oxy_D_avg); caxis([0 1]); title('corr_D');%subplot(1,2,1)1，2划分一行两列个画布，后面的1代表在第一个位置画
%将画布分割成1行2列，在第二列绘制条件2的相关矩阵 将colorbar范围定为0.1-0.5 标题为con2
print(gcf,'corr_D.tif','-dtiff','-r600')
figure;
imagesc(1:27, 1:27, corr_oxy_J_avg); caxis([0 1]); title('corr_J');
print(gcf,'corr_J.tif','-dtiff','-r600')
figure;
imagesc(1:27, 1:27, corr_oxy_P_avg); caxis([0 1]); title('corr_P');
print(gcf,'corr_P.tif','-dtiff','-r600')



%% 

%绘制fisherZ后的相关矩阵图
figure;
subplot(1,2,1); imagesc(1:24, 1:24, corr_oxy_con1_z_avg); caxis([0.1 0.5]);title('con1');colorbar
subplot(1,2,2); imagesc(1:24, 1:24, corr_oxy_con2_z_avg); caxis([0.1 0.5]); title('con2');colorbar
suptitle('corr fisher-Z')
print(gcf,'corr fisher-Z.tif','-dtiff','-r600')
%绘制相干的组平均矩阵
figure;
subplot(1,2,1); imagesc(1:24, 1:24, coh_oxy_con1_avg); caxis([0.1 0.5]);title('con1');colorbar
subplot(1,2,2); imagesc(1:24, 1:24, coh_oxy_con2_avg); caxis([0.1 0.5]);title('con2');colorbar  
suptitle('coh')
print(gcf,'COH.tif','-dtiff','-r600')%保存图面，（'COH.tif','-dtiff','-r600'保存的名称，保存的类型，图片分辨率）

%绘制相位锁值的组平均矩阵
figure;
subplot(1,2,1); imagesc(1:24, 1:24, plv_oxy_con1_avg); caxis([0.1 0.5]);title('con1');
subplot(1,2,2); imagesc(1:24, 1:24, plv_oxy_con2_avg); caxis([0.1 0.5]); title('con2'); 
suptitle('plv')
print(gcf,'plv.tif','-dtiff','-r600')


%% 绘制单个的FC矩阵
clc;clear all;close all
datadir = 'D:\Scientific research\Ana_date\shuju\Homer\NO10HZ\FC';%预处理后数据所在的文件夹
dirinfo = dir([datadir, filesep, '*.txt']);%筛选出文件夹下所有.txt结尾的文件
file_names = {dirinfo.name};%取出文件名放入胞元
for j = 1:length(file_names)
    fig = load(strcat(datadir, filesep,file_names{1,j}));
    figure;
    imagesc(1:27, 1:27,fig); caxis([0 1]); title('P_corr');colorbar

end 


