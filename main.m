clear;close all; clc

  %%%% load ground motion - here an old one from a txt
  damping=0.05;

  
target=makeSpectrum();
target.damping=damping ;
target.time=time;
target.signal=acceleration;

targetSpectrum=target.spectrum;
periodSpectrum=target.period; 




%%% look for the wanted pseudo acceleartion for your calculation
  
index=find(periodSpectrum>=period,1);
Spa=targetSpectrum(index);
disp(['questo è il valore di Spa ' num2str(Spa)])
