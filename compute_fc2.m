 clear all; close all; clc;
%����Ԥ��������������ڵ�·�� ע���޸�
D = 'D:\Scientific research\Ana_date\shuju\A_Homer\A_25Hz 10min oxy\ROI\4ROI\D';
J = 'D:\Scientific research\Ana_date\shuju\A_Homer\A_25Hz 10min oxy\ROI\4ROI\J';
P = 'D:\Scientific research\Ana_date\shuju\A_Homer\A_25Hz 10min oxy\ROI\4ROI\P_10';
%����dir����Ѱ������·������������.mat��β���ļ���Ϊ����forѭ������
D_Files = dir(fullfile(D,'*.mat')); %��ͬ��Ļ�Condition���Ը�Ϊgroup
J_Files = dir(fullfile(J,'*.mat')); 
P_Files = dir(fullfile(P,'*.mat'));
%��ȡÿ�����Ծ�����ļ����浽��Ԫ
FileNamesD = {D_Files.name};
FileNamesJ = {J_Files.name};
FileNamesP = {P_Files.name};

%% compute the measures
fs = 25;  %������  ע���޸�

%������D �Ĺ�������ָ�꣬��oxyʵ��
for sub = 1:length(FileNamesD) %����ÿ������ѭ��
    %�����sub�����Ե�����
    load(strcat(D, '\', FileNamesD{1,sub}));
    
    % correlation
    % corr_oxy_con1 ��ά���� ͨ��*ͨ��*����
    corr_oxy_D(:,:,sub) = corr(oxy);%�������б��Լ���õ��Ĺ������Ӿ���
  %  r = corr(oxy_down);
    %fisher Z  �任��ʹ�������̫�ֲ�
   % corr_oxy_D_z(:,:,sub) = 0.5*log((1+r)./(1-r));
    
    % coherence ���
    %for channel1 = 1:size(oxy,2)  %��һ��ͨ��  size(oxy,2)�����2�ǵ�2��ά��
      %  for channel2 = 1:size(oxy,2) %�ڶ���ͨ��
            %coh_oxy_con1 ��ά���� Ƶ��*ͨ��*ͨ��*����  Freq Ƶ�ʵ�ֲ�  ��ȡ��ӦƵ���ϵ����ָ��
            %[coh_oxy_D(:,channel1,channel2,sub),Freq] = mscohere(oxy(:,channel1),oxy(:,channel2),[],[],[],fs);
     %   end
 %   end     
    
    % phase locking value  ��λ��ֵ
%    for channel1 = 1:size(oxy,2)  %���ŵ�һ��ͨ��
   %     for channel2 = 1:size(oxy,2) %���ŵڶ���ͨ��
            %��ȡ1��ͨ�������ݣ�����ϣ�����ر任����ȡ˲ʱ��λ
 %^           channel1_oxy_phase = angle(hilbert(oxy(:,channel1)));
            %��ȡ2��ͨ�������ݣ�����ϣ�����ر任����ȡ˲ʱ��λ
 %           channel2_oxy_phase = angle(hilbert(oxy(:,channel2)));
            %������λ��
%            rp_oxy = channel1_oxy_phase - channel2_oxy_phase;
            %����PLV ͨ��*ͨ��*����
%            plv_oxy_D(channel1,channel2,sub) = abs(sum(exp(1i*rp_oxy))/length(rp_oxy));
%        end
 %   end
    %����������ӳ���ݴ���Ľ���
    waitbar(sub/length(FileNamesD))
end
%��ȡ����ͨ�����б��Ը���ȤƵ���ϵ�COHֵ[0.01 -- 0.1HZ]  coh_oxy_con1  ��ά����  ͨ��*ͨ��*����
%coh_oxy_D = squeeze(mean(coh_oxy_D(Freq > 0.01 & Freq < 0.1,:,:,:),1));
%������1�Ĺ�������ָ�걣�浽���ش���  FCmaxtrix_con1.mat
save FCmaxtrix_D.mat  corr_oxy_D %corr_oxy_D_z coh_oxy_D  plv_oxy_D

%% 

fs = 25;
%���� ��J �Ĺ�������ָ�꣨��oxy��Ϊʾ����
for sub = 1:length(FileNamesD) %����ÿ������ѭ��
    %�����sub�����Ե�����
    load(strcat(D, '\', FileNamesD{1,sub}));
    
    % correlation
    % corr_oxy_con1 ��ά���� ͨ��*ͨ��*����
    corr_roi_D(:,:,sub) = corr(roi);%�������б��Լ���õ��Ĺ������Ӿ���
 %   r = corr(oxy_down);
    %fisher Z  �任��ʹ�������̫�ֲ�
%    
    % coherence ���
%    for channel1 = 1:size(oxy,2)  %��һ��ͨ��  size(oxy,2)�����2�ǵ�2��ά��
 %       for channel2 = 1:size(oxy,2) %�ڶ���ͨ��
%          %coh_oxy_con1 ��ά���� Ƶ��*ͨ��*ͨ��*����  Freq Ƶ�ʵ�ֲ�  ��ȡ��ӦƵ���ϵ����ָ��
%            [coh_oxy_J(:,channel1,channel2,sub),Freq] = mscohere(oxy(:,channel1),oxy(:,channel2),[],[],[],fs);
 %       end
 %   end     
    
    % phase locking value  ��λ��ֵ
 %   for channel1 = 1:size(oxy,2)  %���ŵ�һ��ͨ��
 %       for channel2 = 1:size(oxy,2) %���ŵڶ���ͨ��
            %��ȡ1��ͨ�������ݣ�����ϣ�����ر任����ȡ˲ʱ��λ
  %          channel1_oxy_phase = angle(hilbert(oxy(:,channel1)));
            %��ȡ2��ͨ�������ݣ�����ϣ�����ر任����ȡ˲ʱ��λ
  %          channel2_oxy_phase = angle(hilbert(oxy(:,channel2)));
            %������λ��
 %           rp_oxy = channel1_oxy_phase - channel2_oxy_phase;
            %����PLV ͨ��*ͨ��*����
   %         plv_oxy_J(channel1,channel2,sub) = abs(sum(exp(1i*rp_oxy))/length(rp_oxy));
  %      end
  %  end
    %����������ӳ���ݴ���Ľ���
    waitbar(sub/length(FileNamesD))
end
%��ȡ����ͨ�����б��Ը���ȤƵ���ϵ�COHֵ[0.01 -- 0.1HZ]  coh_oxy_con1  ��ά����  ͨ��*ͨ��*����
%coh_oxy_J = squeeze(mean(coh_oxy_J(Freq > 0.01 & Freq < 0.1,:,:,:),1));
%������1�Ĺ�������ָ�걣�浽���ش���  FCmaxtrix_con1.mat
save FCmaxtrix_D.mat  corr_roi_D %corr_oxy_J_z coh_oxy_J  plv_oxy_J

%% 

fs = 10;
%������P�Ĺ�������ָ�꣨��oxy��Ϊʾ����
for sub = 1:length(FileNamesP) %����ÿ������ѭ��
    %�����sub�����Ե�����
    load(strcat(P, '\', FileNamesP{1,sub}));
    
    % correlation
    % corr_oxy_con1 ��ά���� ͨ��*ͨ��*����
    corr_roi_P(:,:,sub) = corr(roi);%�������б��Լ���õ��Ĺ������Ӿ���
 %   r = corr(oxy_down);
    %fisher Z  �任��ʹ�������̫�ֲ�
%    corr_oxy_P_z(:,:,sub) = 0.5*log((1+r)./(1-r));
    
    % coherence ���
 %   for channel1 = 1:size(oxy,2)  %��һ��ͨ��  size(oxy,2)�����2�ǵ�2��ά��
 %       for channel2 = 1:size(oxy,2) %�ڶ���ͨ��
            %coh_oxy_con1 ��ά���� Ƶ��*ͨ��*ͨ��*����  Freq Ƶ�ʵ�ֲ�  ��ȡ��ӦƵ���ϵ����ָ��
 %           [coh_oxy_P(:,channel1,channel2,sub),Freq] = mscohere(oxy(:,channel1),oxy(:,channel2),[],[],[],fs);
  %      end
  %  end     
    
    % phase locking value  ��λ��ֵ
%    for channel1 = 1:size(oxy,2)  %���ŵ�һ��ͨ��
 %       for channel2 = 1:size(oxy,2) %���ŵڶ���ͨ��
  %          %��ȡ1��ͨ�������ݣ�����ϣ�����ر任����ȡ˲ʱ��λ
 %           channel1_oxy_phase = angle(hilbert(oxy(:,channel1)));
            %��ȡ2��ͨ�������ݣ�����ϣ�����ر任����ȡ˲ʱ��λ
  %          channel2_oxy_phase = angle(hilbert(oxy(:,channel2)));
            %������λ��
 %           rp_oxy = channel1_oxy_phase - channel2_oxy_phase;
            %����PLV ͨ��*ͨ��*����
   %         plv_oxy_P(channel1,channel2,sub) = abs(sum(exp(1i*rp_oxy))/length(rp_oxy));
 %       end
%   end
    %����������ӳ���ݴ���Ľ���
    waitbar(sub/length(FileNamesP))
end
%��ȡ����ͨ�����б��Ը���ȤƵ���ϵ�COHֵ[0.01 -- 0.1HZ]  coh_oxy_con1  ��ά����  ͨ��*ͨ��*����
% coh_oxy_P = squeeze(mean(coh_oxy_P(Freq > 0.01 & Freq < 0.1,:,:,:),1));
%������1�Ĺ�������ָ�걣�浽���ش���  FCmaxtrix_con1.mat
save FCmaxtrix_P_10.mat  corr_roi_P % corr_oxy_P_z coh_oxy_P  plv_oxy_P

%% ��ͼ��ˮƽ

clc;clear all;close all
%��������1������2�Ĺ�������ָ��
load('FCmaxtrix_D.mat')
load('FCmaxtrix_J.mat')
load('FCmaxtrix_P.mat')
%�Ա���ά��ƽ����������ƽ���Ĺ������Ӿ� ��
%corr  Ƥ��ѷ���ϵ��
corr_oxy_D_avg = mean(corr_oxy_D,3);   corr_oxy_J_avg = mean(corr_oxy_J,3);  corr_oxy_P_avg = mean(corr_oxy_P,3);%ͨ��*ͨ��
%corr fisher-z fisherZ���Ƥ��ѷ���ϵ��
%corr_oxy_con1_z_avg = mean(corr_oxy_con1_z,3);   corr_oxy_con2_z_avg = mean(corr_oxy_con2_z,3);
%coh ���
%coh_oxy_con1_avg = mean(coh_oxy_con1,3);   coh_oxy_con2_avg = mean(coh_oxy_con2,3);
%plv ��λ��ֵ
%plv_oxy_con1_avg = mean(plv_oxy_con1,3);   plv_oxy_con2_avg = mean(plv_oxy_con2,3);

%������ص���ƽ������
figure; %���ɿջ���
%�������ָ��1��2�У��ڵ�һ�л�������1����ؾ��� ��colorbar��Χ��Ϊ0.1-0.5 ����Ϊcon1
imagesc(1:27, 1:27, corr_oxy_D_avg); caxis([0 1]); title('corr_D');%subplot(1,2,1)1��2����һ�����и������������1�����ڵ�һ��λ�û�
%�������ָ��1��2�У��ڵڶ��л�������2����ؾ��� ��colorbar��Χ��Ϊ0.1-0.5 ����Ϊcon2
print(gcf,'corr_D.tif','-dtiff','-r600')
figure;
imagesc(1:27, 1:27, corr_oxy_J_avg); caxis([0 1]); title('corr_J');
print(gcf,'corr_J.tif','-dtiff','-r600')
figure;
imagesc(1:27, 1:27, corr_oxy_P_avg); caxis([0 1]); title('corr_P');
print(gcf,'corr_P.tif','-dtiff','-r600')



%% 

%����fisherZ�����ؾ���ͼ
figure;
subplot(1,2,1); imagesc(1:24, 1:24, corr_oxy_con1_z_avg); caxis([0.1 0.5]);title('con1');colorbar
subplot(1,2,2); imagesc(1:24, 1:24, corr_oxy_con2_z_avg); caxis([0.1 0.5]); title('con2');colorbar
suptitle('corr fisher-Z')
print(gcf,'corr fisher-Z.tif','-dtiff','-r600')
%������ɵ���ƽ������
figure;
subplot(1,2,1); imagesc(1:24, 1:24, coh_oxy_con1_avg); caxis([0.1 0.5]);title('con1');colorbar
subplot(1,2,2); imagesc(1:24, 1:24, coh_oxy_con2_avg); caxis([0.1 0.5]);title('con2');colorbar  
suptitle('coh')
print(gcf,'COH.tif','-dtiff','-r600')%����ͼ�棬��'COH.tif','-dtiff','-r600'��������ƣ���������ͣ�ͼƬ�ֱ��ʣ�

%������λ��ֵ����ƽ������
figure;
subplot(1,2,1); imagesc(1:24, 1:24, plv_oxy_con1_avg); caxis([0.1 0.5]);title('con1');
subplot(1,2,2); imagesc(1:24, 1:24, plv_oxy_con2_avg); caxis([0.1 0.5]); title('con2'); 
suptitle('plv')
print(gcf,'plv.tif','-dtiff','-r600')


%% ���Ƶ�����FC����
clc;clear all;close all
datadir = 'D:\Scientific research\Ana_date\shuju\Homer\NO10HZ\FC';%Ԥ������������ڵ��ļ���
dirinfo = dir([datadir, filesep, '*.txt']);%ɸѡ���ļ���������.txt��β���ļ�
file_names = {dirinfo.name};%ȡ���ļ��������Ԫ
for j = 1:length(file_names)
    fig = load(strcat(datadir, filesep,file_names{1,j}));
    figure;
    imagesc(1:27, 1:27,fig); caxis([0 1]); title('P_corr');colorbar

end 


