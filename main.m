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

spectrumSpa=spectrum.spectrum;
spectrumPeriod=spectrum.period; 


%%% look for the wanted pseudo acceleartion for your calculation
  
index=find(spectrumPeriod>=period,1);
Spa=spectrumSpa(index);
disp(['Spa at ' num2str(period) 's is ' num2str(Spa/g) 'g'])


%%% plot
figure(1); hold on
plot(spectrumPeriod,spectrumSpa/g,'r','LineWidth',1)
plot(period,Spa/g,'.k','MarkerSize',10)
plot([0 period],[Spa/g Spa/g],'--k')
plot([period period],[0 Spa/g],'--k')

axis([0 4 0 inf])
legend({'elastic spectrum','selected period'})
xlabel('Period [s]')
ylabel('Sa [g]')
