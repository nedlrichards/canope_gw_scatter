function [rcvx, tax, jdn ]=HMreadnvx(fname, varargin)
%function [rcvx, tax, jdn ]=HMreadnvx(fname, varargin)
%
% rcvx vector of acoustic data in uPa
% tax  time axis relative to jdn
% jdn  matlab date number 
%

  bits2volts=2.5/2^23;        % A/D full scale is 5 Vpeak-peak or 2.5 Vpeak 
  ampgain=10^(12/20);
  hysens=10^(-168/20);

  if nargin==1
    getdata=1;
  else
    getdata=varargin{1};
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  nccgbl=netcdf.getConstant('NC_GLOBAL');
  ncid=netcdf.open( fname, 'NC_NOWRITE');
  dtype = netcdf.getAtt(ncid,nccgbl,'ChannelType');

  year = netcdf.getAtt(ncid,nccgbl,'year', 'double');
  sernum = netcdf.getAtt(ncid,nccgbl,'serial_number');

  dimid = netcdf.inqDimID(ncid,'nsamples');
  [dimname, nsamples] = netcdf.inqDim(ncid,dimid);

  if getdata
    varid = netcdf.inqVarID(ncid,'samples');
    rcv = netcdf.getVar(ncid,varid,'double');
    Nm=length(rcv);
  else
    rcv=[];
    Nm=nsamples;
  end
  varid = netcdf.inqVarID(ncid,'start_time');
  start_time = netcdf.getVar(ncid,varid,'double');
  varid = netcdf.inqVarID(ncid,'sample_rate');
  sample_rate = netcdf.getVar(ncid,varid,'double');
  varid = netcdf.inqVarID(ncid,'jamset_time');
  jamsettime = netcdf.getVar(ncid,varid,'double');
  varid = netcdf.inqVarID(ncid,'star_jamset_latch_time');
  starjamsetlatch = netcdf.getVar(ncid,varid,'double');

  netcdf.close(ncid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if any(1.05*abs(rcv(:))>2^23)
    %fprintf('Saturation warning\n');
  end

  rcvx=bits2volts*rcv(:).';    %Convert bits to volts
  rcvx=rcvx/ampgain/hysens;   % Convert volts to input pressure (uPa)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fs=10000000/256;   %sampling frequency
  fs=sample_rate;   %sampling frequency
  %fprintf( 'sample_rate: %f\n', sample_rate);
  dt=1/fs;
  fn=fs/2;

  yd0=datenum([year 0 0 0 0 0 ]); %time at beginnning of year
    
  immdly=0.00934;  %imm flat delay from STAR to HM,
                   % measured with 100 ohm loop resistance
  immdly=0.00954;  %imm flat delay from STAR to HM,
                   % measured with 1k ohm loop resistance
  drdly=129/fs;    %ADS1274 data ready delay, from datasheet, fig73
  grpdly=38/fs;    %ADS1274 group  delay, from datasheet, p4

  %fprintf('jamsettime: %9d start_time: %9d ', jamsettime , start_time);
  %fprintf('starjamsetlatch: %14.6f ', starjamsetlatch);
  %fprintf('\n');

  tax=dt*[0:Nm-1];
  tofs=(start_time+drdly-grpdly);
  tofsx=datenum(yd0+tofs/24/3600);
  %datestr(tofsx, 'dd-mmm-yyyy HH.MM.SS.FFF')
  taxhm=tax+tofs;

  jmxtime=datenum(yd0 + jamsettime/24/3600);
  jmdly=starjamsetlatch+immdly-jamsettime;
  taxstar=taxhm+jmdly;

  %{
  load ../HMclock/hmclock
  ix=find([hmclock.sn]==sernum);
  hmc=hmclock(ix);

  xdly=0;
  for iv=1:length(hmc.ev)
    jmtime=hmc.ev(iv).jmtime;
    hmtime=hmc.ev(iv).hmtime;
    jmdly=hmc.ev(iv).jmdly;
    hmdly=hmc.ev(iv).hmdly;

    hmtime=[jmtime; hmtime(:)];
    hmdly=[jmdly; hmdly(:)];

    %datestr(hmtime(1), 'dd-mmm-yyyy HH.MM.SS.FFF')
    %datestr(hmtime(end), 'dd-mmm-yyyy HH.MM.SS.FFF')
    %clf
    %%plot( hmtime-jmtime, hmdly-jmdly, 'o')
    %plot( hmtime-yd0, hmdly, 'o')
    %grid on
    %xlabel( 'yearday');
    %ylabel( 'hmclock error (s)');
    %title([sprintf('SN %d  ', sernum) datestr(jmtime)] )
    %drawnow
    %pause
 
    if jmtime==jmxtime 
      xdly=interp1( hmtime, hmdly, tofsx, 'linear', 'extrap');
      break;
    end
  end

  if xdly
    taxstar=taxhm+xdly;
  end 
  %}

  tofs=ceil(taxstar(1));
  jdn=datenum(yd0+tofs/24/3600);
  tax=taxstar-tofs;

  return

