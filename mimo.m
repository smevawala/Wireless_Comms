%% Wireless Comms Mini Matlab MIMO
%Neema Aggarwal
%Shivam Mevawala
%Nicobitch

close all;

SNR = -4:1:8; %list of SNR values to run algorithm
%intialize vecs
BERc=zeros(length(SNR));
tblen =16; %will handle delay for convolution coder

n=3072; %msg length, must be mult of 6,4, and 128
m=2; %QPSK is 2-QAM whick is BPSK

%use the SNR to calculate EbNo
EbNo_c = SNR -10*log10(log2(m));

hMIMO = comm.MIMOChannel(...
  'PathDelays',                1 * 1e-5,...
  'AveragePathGains',          -4,...
  'SpatialCorrelation',        false,...
  'NumTransmitAntennas',       3,...
  'NumReceiveAntennas',        3,...
  'PathGainsOutputPort',       true);

%use the SNR to calculate EbNo for the normal sytem and the convolutional
%coder
EbNo = SNR -10*log10(log2(m));
% eq = dfe(5, 5, rls(.99)); %Construct a decision feedback equalizer object
% eq.SigConst=qammod(0:1,2); %Set the constellation to 4-qam
bers=zeros(1,10); %initialize bers

%loop over SNR values
for k=1:length(SNR)
X_bin=randi([0 m-1],n,3);

% trellis = poly2trellis(7,[171 133]);
%encode (based on which coderate we want)
% code = convenc(X_bin,trellis);

%modulate
Xm=qammod(X_bin,m);


[Fc, pathGains] = step(hMIMO, Xm);


%add noise
Yc=awgn(Fc, SNR(k),'measured');

pginv=zeros(size(pathGains,1),size(pathGains,3), size(pathGains,4));
PGinv=zeros(3,3);
Ye=zeros(size(Yc));
for kkk=1:n
   PGinv(:,:)=pathGains(kkk,1,:,:);
   Ye=Yc*pinv(PGinv);
end



%demodulate
Zc=qamdemod(Ye,m);
%decode based on coderate
% d = vitdec(Zc,trellis,tblen,'trunc','hard');
%calculate bit error rate
ber=biterr(Zc,X_bin)/(3*n);
BERc(k)=ber;
end

% plot BER
figure

semilogy(EbNo_c, BERc,'kx');

xlabel('EbNo (dB)')
ylabel('BER')
title('Waterfall Plot (Convolutional Coder)')
% legend('theoretical', 'actual')


