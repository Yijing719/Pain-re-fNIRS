%%注意！两次滤波会对数据造成严重的影响，所以使用homer进行预处理后
%如果还想使用nirs_spm进行激活分析，谨慎滤波
%% nirs文件转换成mat文件
clear;clc;close all
datadir = 'D:\Scientific research\Ana_date\shuju\Homer\preprocess.nirs';
cd(datadir)
nirs = dir('*.nirs');
fn = {nirs.name};
savepath =  'D:\Scientific research\Ana_date\shuju\Homer\preprocess.mat';
mkdir(savepath)
for i = 1:length(fn)
    load(fn{i},'-mat');
    nirs_data.oxyData = squeeze(procResult.dc(:,1,:));
    nirs_data.dxyData = squeeze(procResult.dc(:,2,:));
    nirs_data.tHbData = squeeze(procResult.dc(:,3,:));
    nirs_data.fs = 1/(t(2) - t(1));%采样率
    nirs_data.nch = size(procResult.dc,3);%通道数
    nirs_name = [fn{i}(1:end-4),'mat'];
    save([savepath,'\',nirs_name],'nirs_data');
end





