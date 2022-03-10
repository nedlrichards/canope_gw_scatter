%  this code plots the ping intensity and saves specific time period
% HM75 on mooring T1. 156-m from the surface
close all
clear
hr1=0:23;
nc=0;
yrday1=251;yrday2=300;
intens11=zeros(62501,1200);
intens115=zeros(62501,1200);
intens12=zeros(62501,1200);
tt1=zeros(3,1200);
tt2=zeros(3,1200);
direct_intens=zeros(3,1200);
reflected_intens=zeros(3,1200);

for yrday=yrday1:yrday2
for hr=0:23
nc=nc+1
if hr < 10
    fname=['../../data/raw/nav_' num2str(yrday) '0' num2str(hr) '5458.nc'];
else
    fname=['../../data/raw/nav_' num2str(yrday) num2str(hr) '5458.nc'];
end
varargin=1;

% load one file to get some basic data
[rcvx, tax, jdn ]=HMreadnvx(fname, varargin);
dt=tax(2)-tax(1);

% A few bad files
if yrday == 371 & hr == 18
    tax=tax+0.026;
end

if yrday == 437 & hr == 9
    tax=tax+0.026;
end

if yrday == 518 & hr == 0
    tax=tax+0.026;
end

fs=1/dt;
yrday0=jdn-datenum(2016,1,1,0,0,0)+1;
% rcvx vector of acoustic data in uPa
% uPa= 1e6 Pa
% tax  time axis relative to jdn
% jdn  matlab date number 

yearday(:,nc)=yrday0;

% demodulation frequency in kHz
fc=11;
% demodulate and filter
[mf1]=navpro( rcvx, fc);


fc=11.5;
% demodulate and filter
[mf2]=navpro( rcvx, fc);


fc=12;
% demodulate and filter
[mf3]=navpro( rcvx, fc);


% choose the time duration
taxi=5.6:2.5600e-05:7.2; %time 

% compute the intensity
tmp1=abs(interp1(tax-1,mf1,taxi,'linear')).^2; % filtered amplitude
tmp2=abs(interp1(tax-1,mf2,taxi,'linear')).^2; % filtered amplitude
tmp3=abs(interp1(tax-1,mf3,taxi,'linear')).^2; % filtered amplitude
adsa

% filter the intensity
fN=2*pi/(2*dt);
fcc=2*pi/(200*dt);
N=4;
Wn=fcc/fN;
[B,A] = butter(N,Wn);


tmp1=abs(filtfilt(B, A, tmp1)); % 11 kHz
tmp2=abs(filtfilt(B, A, tmp2)); % 11.5 kHz
tmp3=abs(filtfilt(B, A, tmp3)); % 12 kHz


intens11(:,nc)=tmp1;   % all recording as a function of yearday for 11 kHz 
intens115(:,nc)=tmp2;  % all recording as a function of yearday for 11.5 kHz 
intens12(:,nc)=tmp3;   % all recording as a function of yearday for 12 kHz 

intens(1,:)=tmp1;  % 11 kHz
intens(2,:)=tmp2;  % 11.5 kHz
intens(3,:)=tmp3;  % 12 kHz


% find the peaks for each frequency seperately
for i=1:3   
inten=intens(i,:);
% find the peaks
[pks,locs]=findpeaks(inten,'sortstr','descend');
good1=find((inten)==pks(1));
good2=find((inten)==pks(2));
good3=find((inten)==pks(3));
good4=find((inten)==pks(4));

a=taxi(good1);
b=taxi(good2);
c=taxi(good3);
d=taxi(good4);


%find the peaks
if taxi(good1)<taxi(good2) && abs(a-b) >0.05
tt1(i,nc)=taxi(good1);
tt2(i,nc)=taxi(good2);
direct_intens(i,nc)=inten(good1);
reflected_intens(i,nc)=inten(good2);

elseif taxi(good2)<taxi(good1) && abs(a-b) >0.05
tt1(i,nc)=taxi(good2);
tt2(i,nc)=taxi(good1);
direct_intens(i,nc)=inten(good2);
reflected_intens(nc,i)=inten(good1);

elseif taxi(good1)<taxi(good3) && abs(a-c) >0.05
tt1(i,nc)=taxi(good1);
tt2(i,nc)=taxi(good3);
direct_intens(i,nc)=inten(good1);
reflected_intens(i,nc)=inten(good3);


elseif taxi(good3)<taxi(good1) && abs(c-a) >0.05
tt1(i,nc)=taxi(good3);
tt2(i,nc)=taxi(good1);
direct_intens(i,nc)=inten(good3);
reflected_intens(i,nc)=inten(good1);

elseif taxi(good1)<taxi(good4) && abs(a-d) >0.05
tt1(i,nc)=taxi(good1);
tt2(i,nc)=taxi(good4);
direct_intens(i,nc)=inten(good1);
reflected_intens(i,nc)=inten(good4);

elseif taxi(good4)<taxi(good1) && abs(a-d) >0.05
tt1(i,nc)=taxi(good4);
tt2(i,nc)=taxi(good1);
direct_intens(i,nc)=inten(good4);
reflected_intens(i,nc)=inten(good1);
end

end
end
end

%tt1(1,:) & tt2(1,:) arrival time for direct and reflected path for 11 kHz
%tt1(2,:) & tt2(2,:) arrival time for direct and reflected path for 11.5 kHz
%tt1(3,:) & tt2(3,:) arrival time for direct and reflected path for 12 kHz


save nav75_25130011_5458.mat intens11 taxi yearday tt1 tt2 direct_intens reflected_intens
save nav75_251300115_5458.mat intens115 taxi yearday tt1 tt2 direct_intens reflected_intens
save nav75_25130012_5458.mat intens12 taxi yearday tt1 tt2 direct_intens reflected_intens
%%

% choose one point to show the line plot
figure(4)
plot(taxi,intens11(:,1))
xlim([6 7])


figure(5)
good1=find(tt1(1,:)>0); % 11 kHz
good2=find(tt1(2,:)>0); % 11.5 kHz
good3=find(tt1(3,:)>0); % 12 kHz

plot(yearday(good1),tt1(1,good1),'b')
hold on
plot(yearday(good2),tt1(2,good2),'r')
hold on
plot(yearday(good3),tt1(3,good3),'g')
%xlim([min(yearday) max(yearday)])
xlim([280 290])
grid off
ylim([6 6.5])
legend('11','11.5','12')
%%
load('mooring.mat')
% Distance between interrogator and transponder
L0=L0(3,:); %HM75

% Distance along the direct path
L1=L1(3,:); %HM75

% Distance along the reflected path
L2=L2(3,:); %HM75

% depth
z=rcvz(3,:); %HM75
%%
fc=11;
% clean time with no interference
if fc == 11
tc1=414;tc2=416;
ii1=find( yearday > tc1 & yearday < tc2);

tc1=493;tc2=496;
ii2=find( yearday > tc1 & yearday < tc2);

tc1=542.5;tc2=544.5;
ii3=find( yearday > tc1 & yearday < tc2);
end

if fc == 11.5
% clean time with no interference
tc1=268.5;tc2=270.5;
ii1=find( yearday > tc1 & yearday < tc2);

tc1=284.5;tc2=286.75;
ii2=find( yearday > tc1 & yearday < tc2);

tc1=416;tc2=418;
ii3=find( yearday > tc1 & yearday < tc2);

tc1=598;tc2=600;
ii4=find( yearday > tc1 & yearday < tc2);
end

if fc == 12
tc1=269;tc2=272;
ii1=find( yearday > tc1 & yearday < tc2);

tc1=286;tc2=289;
ii2=find( yearday > tc1 & yearday < tc2);

tc1=419;tc2=422;
ii3=find( yearday > tc1 & yearday < tc2);

tc1=600;tc2=602;
ii4=find( yearday > tc1 & yearday < tc2);
end


