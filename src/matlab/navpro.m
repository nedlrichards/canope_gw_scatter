function [mf]=navpro( rcv, fc)
%
% AFAIK: this is the processing that DH uses with STAR navigation data. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %signal replica
  fc=1000*fc;

  fs=10000000/256;   %sampling frequency
  dt=1/fs;

  T=0.009;     %pulse length
  N=floor(T/dt);
  taxis=dt*[0:N];
  sg=sin(2*pi*fc*taxis);  %signal replica
  sgb=ones(1,N);  %signal replica
  %sgb=kaiser(N,2.5).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%design demodulation filter
  Wn=2000.0;    %demodulation filter cutoff
  uwarp = 2*fs*tan(pi*Wn/fs);
  [zero,pole,gain]   = besselap(3);
  [ssa,ssb,ssc,ssd]  = zp2ss(zero,pole,gain);
  [ssa,ssb,ssc,ssd]  = lp2lp(ssa,ssb,ssc,ssd,uwarp);
  [ssa,ssb,ssc,ssd]  = bilinear(ssa,ssb,ssc,ssd,fs);
  fira = poly(ssa);
  firb = poly(ssa-ssb*ssc)+(ssd-1)*fira;

  %clf
  %subplot( 211)
  %[H,fax]=freqz( firb, fira, 1000, fs);
  %semilogx( fax, 20*log10(abs(H)));
  %grid on
  %set( gca, 'xlim', [0 1000]);
  %set( gca, 'ylim', [-90 10]);

  %subplot( 212)
  %[irp,tax]=impz( firb, fira, 1000, fs);
  %plot( tax, irp);
  %grid on
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  rcv=rcv(:).';
  Nm=length(rcv);
  tax=dt*[0:Nm-1];

  %mf=filter( sg, 1, rcv);
  %mfr=mf.*cos(-2*pi*fc*tax);
  %mfi=mf.*sin(-2*pi*fc*tax);
  %mfr=filter(firb, fira, mfr);
  %mfi=filter(firb, fira, mfi);
  %mf=(mfr.*mfr + mfi.*mfi);
 
  
  
 %%design demodulation filter
%   Wn=2000.0;    %demodulation filter cutoff
%   uwarp = 2*fs*tan(pi*Wn/fs);
%   [zero,pole,gain]   = besselap(3);
%   [ssa,ssb,ssc,ssd]  = zp2ss(zero,pole,gain);
%   [ssa,ssb,ssc,ssd]  = lp2lp(ssa,ssb,ssc,ssd,uwarp);
%   [ssa,ssb,ssc,ssd]  = bilinear(ssa,ssb,ssc,ssd,fs);
%   fira = poly(ssa); , %returns the coefficients of the polynomial whose roots are the elements of r
%   firb = poly(ssa-ssb*ssc)+(ssd-1)*fira;
%   
  dmx=rcv.*exp(-j*2*pi*fc*tax);
  dm=filter(firb, fira, dmx); %filters the input data x using a rational transfer function
  %mf = dm;
  mf=filter( sgb, 1, dm);
  %mf=abs(mf).^2;

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

