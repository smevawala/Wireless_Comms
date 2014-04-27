%% Wireless Comms Mini Matlab 1
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
nsubc=64; % number of subcarriers

% coderate = 1/2;

%use the SNR to calculate EbNo
EbNo_c = SNR -10*log10(log2(m));


delayVector = [0 1 2 3 4] * 1e-5; % Discrete delays of four-path channel (s)
gainVector = [0 -4 -6 -9 -14];
rchan_sel=rayleighchan(1e-5,1, delayVector, gainVector);%Set up the
%channel fading object with delay and gain vecs
rchan_sel.StoreHistory = 1;
rchan_sel.StorePathGains = 1;

%use the SNR to calculate EbNo for the normal sytem and the convolutional
%coder
EbNo = SNR -10*log10(log2(m));
eq = dfe(5, 5, rls(.99)); %Construct a decision feedback equalizer object
eq.SigConst=qammod(0:1,2); %Set the constellation to 4-qam
bers=zeros(1,10); %initialize bers

%loop over SNR values
for k=1:length(SNR)
X_bin=randi([0 m-1],1,n);

% trellis = poly2trellis(7,[171 133]);
%encode (based on which coderate we want)
% code = convenc(X_bin,trellis);

%modulate
Yc=qammod(X_bin,m);

% parallelize
Pc=zeros(length(Yc)/nsubc,nsubc); 
for kk=1:length(Yc)/nsubc
Pc(kk,:)=Yc(1,((kk-1)*nsubc+1):((kk)*nsubc));
end
% IFFT
ifft_sig=ifft(Pc',nsubc)';

% Adding Cyclic Extension to enable circular conv

cext_data=zeros(length(Yc)/nsubc,nsubc+16);
cext_data(:,(1:16))=ifft_sig(:,(nsubc-15:nsubc));
for i=1:nsubc
    
    cext_data(:,(i+16))=ifft_sig(:,i);
    
end

%filter through chan
Fc = zeros(size(cext_data));
PG = zeros(size(cext_data,2), 5, size(cext_data,1));
for kkk=1:size(cext_data, 1)
    Fc(kkk, :)=filter(rchan_sel,cext_data(kkk,:));
    PG(:,:,kkk)=rchan_sel.PathGains;
end

%add noise
Ac=awgn(Fc, SNR(k),'measured');



%Removing Cyclic Extension
rxed_sig=zeros(length(Yc)/nsubc,nsubc);
for i=1:nsubc
    rxed_sig(:,i)=Ac(:,i+16);    
end


% FFT
ff_sig=fft(rxed_sig',nsubc)';

%serialize
Sc=zeros(1,length(Yc));
for kkk=1:length(Yc)/nsubc
    
    Sc((kkk-1)*nsubc+1:(kkk)*nsubc)=ff_sig(kkk,:);
end

%demodulate
Zc=qamdemod(Sc,m);
%decode based on coderate
% d = vitdec(Zc,trellis,tblen,'trunc','hard');
%calculate bit error rate
ber=biterr(Zc,X_bin)/(n);
BERc(k)=ber;
end

% plot BER
figure

semilogy(EbNo_c, BERc,'kx');

xlabel('EbNo (dB)')
ylabel('BER')
title('Waterfall Plot (Convolutional Coder)')
% legend('theoretical', 'actual')


