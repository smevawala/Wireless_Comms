%% Wireless Comms mini Matlab 1
%Neema Aggarwal
%Shivam Mevawala
%Nicobitch

close all;
rate=0;

SNR = -4:1:8; %list of SNR values to run algorithm
%intialize vecs
% BER=zeros(length(SNR));
BERc=zeros(length(SNR));
tblen =32; %will handle delay for convolution coder

n=3000; %number of samples
m=4; %QPSK is 4-QAM
punc= [1 1 1 0 0 1];
if rate
    coderate = 3/4;
else
    coderate = 1/2;
end
%use the SNR to calculate EbNo for the normal sytem and the convolutional
%coder
EbNo = SNR -10*log10(log2(m));
EbNo_c = SNR -10*log10(log2(m)*coderate);

%loop over SNR values
for k=1:length(SNR)
%generate a random vector of 4 symbols
X=randi([0 m-1],1,n);
%modulate
% Y=qammod(X,m);
%add noise
% A=awgn(Y, SNR(k),'measured');
%demodulate
% Z=qamdemod(A,m);
%calculkate bit error rate
% ber=biterr(Z,X)/(2*n);
% BER(k)=ber;

% Convolutional Coder
%convert symbols to binary
X_bin = reshape((de2bi(X, 2,'left-msb')).',1,n*2); 
%define a trellis (default chosen) with coderate .5
trellis = poly2trellis(7,[171 133]);
%encode
if rate
    code = convenc(X_bin,trellis, punc);
else
    code = convenc(X_bin,trellis);
end
%modulate
Yc=qammod(bin2dec([num2str(code(1:2:end-1)') num2str(code(2:2:end)')])',m);
%add noise
Ac=awgn(Yc, SNR(k),'measured');
%demod
Zc=reshape(de2bi(qamdemod(Ac,m),2,'left-msb').',1,length(Ac)*2);
%decode
if rate
    d = vitdec(Zc,trellis,tblen,'trunc','hard', punc);
else
    d = vitdec(Zc,trellis,tblen,'trunc','hard');
end
%calculate bit error rate
ber=biterr(d,X_bin)/(2*n);
BERc(k)=ber;
end

%plots
% figure
% semilogy(EbNo,berawgn(EbNo, 'qam', 4),'b-');
% hold on;
% semilogy(EbNo,BER, 'rx');
% 
% xlabel('EbNo (dB)')
% ylabel('BER')
% title('Waterfall Plots')
% legend('theoretical', 'actual')


figure
SPECT = distspec(trellis,7);
semilogy(EbNo_c,bercoding(EbNo_c, 'conv', 'hard', 1/2, SPECT),'m-');
hold on;
semilogy(EbNo_c, BERc,'kx');

xlabel('EbNo (dB)')
ylabel('BER')
title('Waterfall Plots (Convolutional Coder)')
legend('theoretical', 'actual')


