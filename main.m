clear;close all; clc

  
  g=9.81;
  damping=0.05;
  period=0.85;
  
  %%%% load ground motion - here an old one from a .dat
  temp = load('Kobe95.dat');
  time=temp(:,1);
  signal=temp(:,2);

  
spectrum=makeSpectrum();
spectrum.damping=damping ;
spectrum.time=time;
spectrum.signal=signal;

spectrumSpa=target.spectrum;
spectrumPeriod=target.period; 


%%% look for the wanted pseudo acceleartion for your calculation
  
index=find(spectrumPeriod>=period,1);
Spa=spectrumSpa(index);
disp([Spa at ' num2str(period) 's is ' num2str(Spa/g) 'g'])


%%% plot
plot(spectrumPeriod,spectrumSpa/g*,'r--')
plot(period,Spa/g*,'.')

axis([0 4 0 inf])
legend({'elastic spectrum','selected period'})
xlabel('Period [s]')
ylabel('Sa [g]')
