% cd ('/Users/clairezhao/Documents/data/MRI/pv6/CEST/glucoCEST-coop-heme/Untitled Folder/20191008_050744_ANEST_71983_1_1/10/');
clear all;close all; clc;
% 
load('AcqAndLsqFit.mat');
load('T1Map_lsq.mat');

% load('T2AcqDataAndFit');%T2

%debug plot
h = figure,imagesc(squeeze(T1Map_lsq(:,:,1))),colormap jet;
colormap jet;title('T1Map');
cb = colorbar;
title(cb, 'ms');

i=0;
while i==0
clear ii jj
[jj,ii]=getpts(h);

ii=int8(ii);
jj=int8(jj);

% fun=@(A,X)A(1)*(1-2*A(3)*exp(-X./A(2)));%T1Fit
fun=@(A,X)A(1).*(1-2.*exp(-X./A(2)));%brukerfunction
% fun=@(A,X)A(1)*(1-2*exp(-X./A(2)))+A(3);%T1Fit for 71669_1015
% fun = @(A,X)A(1)*exp(-X./A(2))+A(3);%T2Fit

%plot for single fitting
% figure,plot(EffectiveTE,T2AcqData(:,1,ii,jj),'o');%T2
figure,plot(TI',squeeze(AcqData(:,1,ii,jj)),'o');

hold on
% fitResultY=fun(T2fitResults(:,1,ii,jj),linspace(EffectiveTE(1),EffectiveTE(end),100));%T2
fitResultY=fun(fitResults(:,1,ii,jj),linspace(TI(1),TI(end),100));
% T2Value=T2fitResults(2,1,ii,jj);
T1Value=fitResults(2,1,ii,jj);
% plot(linspace(EffectiveTE(1),EffectiveTE(end),100),fitResultY,'r-');%T2

plot(linspace(TI(1),TI(end),100),fitResultY,'-');
% title(sprintf('wt#71670_1_0_0_8Heme @(%d,%d), T2=%1f',ii,jj,T2Value));%T2
title(sprintf('ko#71670_1_0_0_8Heme @(%d,%d), T1=%3f',ii,jj,T1Value));

xlabel('TI (ms)');
ylabel('Intensity (a.u.)')
%  axis([0 80 0 30]);% 71669_1015Base axis
end

