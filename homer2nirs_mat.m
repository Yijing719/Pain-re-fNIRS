%%ע�⣡�����˲��������������ص�Ӱ�죬����ʹ��homer����Ԥ�����
%�������ʹ��nirs_spm���м�������������˲�
%% nirs�ļ�ת����mat�ļ�
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
    nirs_data.fs = 1/(t(2) - t(1));%������
    nirs_data.nch = size(procResult.dc,3);%ͨ����
    nirs_name = [fn{i}(1:end-4),'mat'];
    save([savepath,'\',nirs_name],'nirs_data');
end





