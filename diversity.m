%% Wireless Comms mini Matlab 2
%Neema Aggarwal
%Shivam Mevawala
%Nicolas Castro

close all;
clear all;
SNR = -4:2:20; %list of SNR values to run algorithm
%intialize vecs
BER_reg=zeros(length(SNR));
BER_mrrc2=zeros(length(SNR));
BER_mrrc4=zeros(length(SNR));


n=307200; %number of samples
m=2; %BPSK is 2-QAM

% delayVector = [0 1 2 3 4] * 1e-5; % Discrete delays of four-path channel (s)
% gainVector = [0 -4 -6 -9 -14];

rchan_flat=rayleighchan(1e-5,1e4); %Set up the channel fading object
rchan_flat.StoreHistory = 1; %If this value is 1, channel state information
%needed by the channel visualization tool is stored. 
rchan_flat.StorePathGains = 1;%If set to 1, the complex path gain vector is
%stored.
rchan_flat2=rayleighchan(1e-5,1e4); %Set up the channel fading object
rchan_flat2.StoreHistory = 1; %If this value is 1, channel state information
%needed by the channel visualization tool is stored. 
rchan_flat2.StorePathGains = 1;%If set to 1, the complex path gain vector is
%stored.

rchan_flat3=rayleighchan(1e-5,1e4); %Set up the channel fading object
rchan_flat3.StoreHistory = 1; %If this value is 1, channel state information
%needed by the channel visualization tool is stored. 
rchan_flat3.StorePathGains = 1;%If set to 1, the complex path gain vector is
%stored.
rchan_flat4=rayleighchan(1e-5,1e4); %Set up the channel fading object
rchan_flat4.StoreHistory = 1; %If this value is 1, channel state information
%needed by the channel visualization tool is stored. 
rchan_flat4.StorePathGains = 1;%If set to 1, the complex path gain vector is
%stored.

rchan_flat5=rayleighchan(1e-5,1e4); %Set up the channel fading object
rchan_flat5.StoreHistory = 1; %If this value is 1, channel state information
%needed by the channel visualization tool is stored. 
rchan_flat5.StorePathGains = 1;%If set to 1, the complex path gain vector is
%stored.

rchan_flat6=rayleighchan(1e-5,1e4); %Set up the channel fading object
rchan_flat6.StoreHistory = 1; %If this value is 1, channel state information
%needed by the channel visualization tool is stored. 
rchan_flat6.StorePathGains = 1;%If set to 1, the complex path gain vector is
%stored.

rchan_flat7=rayleighchan(1e-5,1e4); %Set up the channel fading object
rchan_flat7.StoreHistory = 1; %If this value is 1, channel state information
%needed by the channel visualization tool is stored. 
rchan_flat7.StorePathGains = 1;%If set to 1, the complex path gain vector is
%stored.

rchan_flat8=rayleighchan(1e-5,1e4); %Set up the channel fading object
rchan_flat8.StoreHistory = 1; %If this value is 1, channel state information
%needed by the channel visualization tool is stored. 
rchan_flat8.StorePathGains = 1;%If set to 1, the complex path gain vector is
%stored.

rchan_flat9=rayleighchan(1e-5,1e4); %Set up the channel fading object
rchan_flat9.StoreHistory = 1; %If this value is 1, channel state information
%needed by the channel visualization tool is stored. 
rchan_flat9.StorePathGains = 1;%If set to 1, the complex path gain vector is
%stored.

rchan_flat10=rayleighchan(1e-5,1e4); %Set up the channel fading object
rchan_flat10.StoreHistory = 1; %If this value is 1, channel state information
%needed by the channel visualization tool is stored. 
rchan_flat10.StorePathGains = 1;%If set to 1, the complex path gain vector is
%stored.


%use the SNR to calculate EbNo for the normal sytem and the convolutional
%coder
EbNo = SNR -10*log10(log2(m));


%loop over SNR values
for k=1:length(SNR)
    %generate a random vector of 4 symbols
    X=randi([0 m-1],1,n);
    %modulate
    Y=qammod(X,m);
    Y0=zeros(1,length(Y)/2);
    Y1=zeros(1,length(Y)/2);
    
    for kk=1:length(Y)
        if mod(kk,2) == 0 % number is even
            Y1(kk-1)=Y(kk);
            Y0(kk)=-conj(Y(kk));
        else
            %number is odd
            Y0(kk)=Y(kk);
            Y1(kk+1)=conj(Y(kk));
        end
        
    end
% 
%     for kk=1:length(Y)
%         if mod(kk,2) == 0 % number is even
%             Y1(kk/2)=Y(kk);
%         else
%             %number is odd
%             Y0(kk/2+.5)=Y(kk);
%         end
%         
%     end
    %add noise and rayleigh fading
    A=filter(rchan_flat,Y);
    A2=filter(rchan_flat2,Y);
    A3=filter(rchan_flat3,Y);
    A4=filter(rchan_flat4,Y);    
    
    filter(rchan_flat5,Y0); %instantiate the Path Gains
    filter(rchan_flat6,Y1);   
%     Anewc0=filter(rchan_flat6,conj(Y0));
%     Anewc1=filter(rchan_flat5,-conj(Y1));
    
    Ag = awgn(A, SNR(k),'measured');
    Ag2 = awgn(A2, SNR(k),'measured'); 
    Ag3 = awgn(A3, SNR(k),'measured'); 
    Ag4 = awgn(A4, SNR(k),'measured'); 
    
    h0 = rchan_flat5.PathGains.'; 
    h1 = rchan_flat6.PathGains.';
    for kkk=1:(length(Y)/2)
        h0(2*kkk) = h0(2*kkk-1);
        h1(2*kkk) = h1(2*kkk-1);
    end

    Anew0=h0.*Y0;
    Anew1=h1.*Y1; 
    
    Rnew=awgn(Anew0+Anew1,SNR(k),'measured');
    DNew=zeros(1,length(Y));
    for kkk=1:(length(Y)/2)
        r0=Rnew(2*kkk-1);
        r1=Rnew(2*kkk);

        Dnew(2*kkk-1)= conj(h0(2*kkk-1))*r0 + h1(2*kkk)*conj(r1);
        Dnew(2*kkk)= conj(h1(2*kkk-1))*r0-h0(2*kkk)*conj(r1);
    end
    
    Ae=Ag.*conj(rchan_flat.PathGains.');
    Ae2=Ag2.*conj(rchan_flat2.PathGains.');
    Ae3=Ag3.*conj(rchan_flat3.PathGains.');
    Ae4=Ag4.*conj(rchan_flat4.PathGains.');
   
    %demodulate
    
    Z=qamdemod(Ae,m);
    Z_mrrc2=qamdemod(Ae+Ae2,m);
    Z_mrrc4=qamdemod(Ae+Ae2+Ae3+Ae4,m);
    Z_new=qamdemod(Dnew,m);
    
    %calculate bit error rate
    BER_reg(k)=biterr(Z,X)/(n);
    BER_mrrc2(k)=biterr(Z_mrrc2,X)/(n);
    BER_mrrc4(k)=biterr(Z_mrrc4,X)/(n);
    BER_new(k)=biterr(Z_new,X)/(n);
end

%plots
figure 
semilogy(EbNo, BER_reg,'kx');
hold on;
semilogy(EbNo, BER_mrrc2,'go');
semilogy(EbNo, BER_mrrc4,'rx');
semilogy(EbNo, BER_new,'bo');
xlabel('EbNo (dB)')
ylabel('BER')
title('Waterfall Plots- Flat Fading')
legend('no diversity', 'MRRC 2 RX', 'MRRC 4 RX','new 2tr')

