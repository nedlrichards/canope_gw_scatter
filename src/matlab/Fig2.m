clear
close all
load('xpos.Canape2016-T1-11.0.mat')
load('DVLAmeanSoundSpeed.mat')
% lat and long
lat=xpos.xlat;
lon=xpos.xlon;
zt=xpos.z;
Re=6356.75*1e3;
xt=Re*cos(75.3633*pi/180)*((lon+145.0522)*pi/180);
yt=Re*((lat-75.3633)*(pi/180));

%load the interragator and hydprohne data
load('nav75_25130011_5458.mat')


intensity=intens11';
title('OW/IF')
nc=0;
fc=11;
imax=max(max(intensity));

for i=1:length(yearday)
    nc=nc+1;
    hold on
    imagesc(taxi,yearday(i),10*log10(abs(intensity(i,:))/imax))
    hold on
end
xlim([5.6 7.2])
ylim([yearday(1) yearday(end)])
cbh=colorbar
colorbar off
colormap('jet')
set(gca, 'clim', [-50 0]);
axis('ij')
ylabel(' Yearday, t','FontSize',17)
%xlabel('Arrival Time (s) ','FontSize',20)
%title([num2str(fc) ' (kHz), T1'],'FontSize',25)
ylim([250 300])
%ylim([251 609])
xlim([6 7.2])
%set(cbh,'XTickLabel',[0:-10:-50])
set(gca,'FontSize',17);
set(gca,'YTick',[250:10:300],'FontSize', 17)
set(gca,'XTick',[6:0.2:7.2],'FontSize', 17)
hold on
plot(tt2(1,:),yearday,'k*','MarkerSize',0.5)
set(gca,'FontSize',17);


