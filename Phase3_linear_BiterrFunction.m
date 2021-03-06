 % Phase 3: Basic error control coding to improve the performance

clear all; close all; clc; 

% Var Declarations
fc = 10000;                                     % carrier freq
fs = 16*fc;                                     % sample freq
bbDataRate = 1000;                              % baseband DataRate 1kbps
nBits = 1024;                                   % data bits
SNRdb = -15:1:15;                               % SNR in db
SNR = (10.^(SNRdb/10));                         % convert SNR from db to dec
samplePeriod = fs/bbDataRate;
encoded_nBits = nBits/4*7;                      % Data bits for linear Code (7,4) 

sigLen = fs*nBits/bbDataRate +1;                % signal length
sigLen_encoded = fs*encoded_nBits/bbDataRate +1; % Encoded signal length
amp = 2;                                        % amplitude
t = 0: 1/fs : nBits/bbDataRate;                 % time intervals --not encoded
t_encoded = 0: 1/fs : encoded_nBits/bbDataRate;  % time intervals --encoded

carrier = amp .* cos(2*pi*fc*t);
carrierLow_encoded = amp .* cos(2*pi*fc*t_encoded);             % for BFSK 1st lower freq carrier
carrierHigh_encoded = amp .* cos(2*pi*5*fc*t_encoded);          % for BFSK 2nd higher freq carrier
carrier_encoded = amp .* cos(2*pi*fc*t_encoded);      %For Encoded Carrier 

[b,a] = butter(6,0.2);                          % 6th Order low pass filter coefficients


% create arrays to store Bit Error Rates BER
OOKerrorArr = zeros(1,length(SNR));             
encoded_OOKerrorArr = zeros(1,length(SNR));      
encoded_BPSKerrorArr = zeros(1,length(SNR));     
encoded_BFSKerrorArr = zeros(1,length(SNR));       

runCycles = 20;

plotSNRdb = 5;                                          % change in multiple of 5 to match; for plotting

for i = 1 : length(SNR)                                 % loop for diff SNR values
        aveOOKerror = 0;
        encoded_aveOOKerror = 0;
        encoded_aveBFSKerror = 0;
        encoded_aveBPSKerror = 0;   

    for j = 1 : runCycles                               % loop to cal ave error rate over runCycles
        data = round(rand(1,nBits));                    % data Gen
        pol = cyclpoly(7,4);
        parmat = cyclgen(7,pol);
        genmat = gen2par(parmat);                       
        linear_sig = encode(data,7,4,'linear/binary', genmat); %linear data
        % signal Gen
        sig = zeros(1,sigLen);
        sig_encoded = zeros(1,sigLen_encoded);
        for k = 1 : sigLen - 1
            sig(k) = data(ceil(k*bbDataRate/fs));
        end
        for k = 1 : sigLen_encoded - 1
            sig_encoded(k) = linear_sig(ceil(k*bbDataRate/fs));
        end
        sig(sigLen) = sig(sigLen-1);
        sig_encoded(sigLen_encoded) = sig_encoded(sigLen_encoded - 1);
        
        %***** On-Off Keying (OOK) *****
           % modulation
        encoded_OOKsig = sig_encoded .* carrier_encoded;                        % if 1 = carrier; if 0 = 0;
        OOKsig = sig .* carrier;
           
            % noise Gen
        OOKsigPwr = bandpower(OOKsig);                  % sum of absolute square/length
        encoded_OOKsigPwr = bandpower(encoded_OOKsig);

        OOKnoisePwr = OOKsigPwr ./ SNR(i);
        OOKnoise = sqrt(OOKnoisePwr/2) .* randn(1,sigLen); %%%% align with Phase1
        encoded_OOKnoisePwr = encoded_OOKsigPwr ./ SNR(i);
        encoded_OOKnoise = sqrt(encoded_OOKnoisePwr/2) .* randn(1,sigLen_encoded); 

            % transmission
        OOKreceived = OOKsig + OOKnoise;
        encoded_OOKreceived = encoded_OOKsig + encoded_OOKnoise;

            % demodulation
        OOKdemod = OOKreceived .* (2.* carrier);        %coherent
        OOKdemod = filtfilt(b, a, OOKdemod);
        encoded_OOKdemod = encoded_OOKreceived .* (2.* carrier_encoded);        %coherent
        encoded_OOKdemod = filtfilt(b, a, encoded_OOKdemod);

            % sampling and decision logic
        OOK = samplingANDdecision(OOKdemod, samplePeriod, nBits, amp/2);
        encoded_OOK = samplingANDdecision(encoded_OOKdemod, samplePeriod, encoded_nBits, amp/2);
        
            %decoding linear Code
        OOK_decoded = decode(encoded_OOK,7,4,'linear/binary', genmat);



        %***** Binary Phase Shift Keying (BPSK) *****
            % modulation
        encoded_BPSKsig = sig_encoded .* 2 - 1;                         % if 1 -> 1*2-1 = 1; if 0 -> 0*2-1 = -1
        encoded_BPSKsig = encoded_BPSKsig .* carrier_encoded;

            % noise Gen
        encoded_BPSKsigPwr = bandpower(encoded_BPSKsig);
        encoded_BPSKnoisePwr = encoded_BPSKsigPwr ./ SNR(i);
        encoded_BPSKnoise = sqrt(encoded_BPSKnoisePwr/2) .* randn(1,sigLen_encoded);

            % transmission
        encoded_BPSKreceived = encoded_BPSKsig + encoded_BPSKnoise;

            % demodulation
        encoded_BPSKdemod = encoded_BPSKreceived .* (2.* carrier_encoded);
        encoded_BPSKdemod = filtfilt(b, a, encoded_BPSKdemod);

            % sampling and decision logic
        encoded_BPSK = samplingANDdecision(encoded_BPSKdemod, samplePeriod, encoded_nBits, 0);

            %decoding linear Code
        BPSK_decoded = decode(encoded_BPSK,7,4,'linear/binary', genmat);
        
        %***** Binary Frequency Shift Keying *****
            % modulation
        encoded_BFSKsigLow = carrierLow_encoded .* (sig_encoded == 0);
        encoded_BFSKsigHigh = carrierHigh_encoded .* (sig_encoded == 1);
        encoded_BFSKsig = encoded_BFSKsigLow + encoded_BFSKsigHigh;

            % noise Gen
        encoded_BFSKsigPwr = bandpower(encoded_BFSKsig);
        encoded_BFSKnoisePwr = encoded_BFSKsigPwr ./ SNR(i);
        encoded_BFSKnoise = sqrt(encoded_BFSKnoisePwr/2) .* randn(1,sigLen_encoded);

            % transmission
        encoded_BFSKreceived = encoded_BFSKsig + encoded_BFSKnoise;

            % demodulation
        encoded_BFSKdemodLow = encoded_BFSKreceived .* (2.* carrierLow_encoded);
        encoded_BFSKdemodLow = filtfilt(b, a, encoded_BFSKdemodLow);
        encoded_BFSKdemodHigh = encoded_BFSKreceived .* (2.* carrierHigh_encoded);
        encoded_BFSKdemodHigh = filtfilt(b, a, encoded_BFSKdemodHigh);
        encoded_BFSKdemod = encoded_BFSKdemodHigh - encoded_BFSKdemodLow;

            % sampling and decision logic
        encoded_BFSK = samplingANDdecision(encoded_BFSKdemod, samplePeriod, encoded_nBits, 0);

             %decoding linear Code
        BFSK_decoded = decode(encoded_BFSK,7,4,'linear/binary', genmat);
        
        %Bit error
        encoded_OOKerror =  biterr(OOK_decoded, data) ./encoded_nBits;
        OOKerror = biterr(OOK, data) ./nBits;
        encoded_BPSKerror = biterr(BPSK_decoded, data) ./encoded_nBits;
        encoded_BFSKerror = biterr(BFSK_decoded, data) ./encoded_nBits;
    
        encoded_aveOOKerror = encoded_OOKerror + encoded_aveOOKerror;
        aveOOKerror = OOKerror + aveOOKerror ;
        encoded_aveBPSKerror = encoded_BPSKerror + encoded_aveBPSKerror;
        encoded_aveBFSKerror = encoded_BFSKerror + encoded_aveBFSKerror;

    end     % end for runCycles loop
 
    OOKerrorArr(i) = aveOOKerror/runCycles;
    encoded_OOKerrorArr(i) = encoded_aveOOKerror/runCycles;
    encoded_BPSKerrorArr(i) = encoded_aveBPSKerror/runCycles;
    encoded_BFSKerrorArr(i) = encoded_aveBFSKerror/runCycles;
    
    % store variable for specified SNR value for plotting
    if (SNRdb(i) == plotSNRdb)
        dataPlot = data;
        linearData = linear_sig;
        sigPlot = sig;
        linearsigPlot = sig_encoded;
        modOOKplot = encoded_OOKsig;
        recOOKplot = encoded_OOKreceived;
        demodOOKplot = encoded_OOK;
        decodeOOKplot = OOK_decoded;

        modBPSKplot = encoded_BPSKsig;
        recBPSKplot = encoded_BPSKreceived;
        demodBPSKplot = encoded_BPSK;
        decodeBPSKplot = BPSK_decoded;

        modBFSKplot = encoded_BFSKsig;
        recBFSKplot = encoded_BFSKreceived;
        demodBFSKplot = encoded_BFSK;
        decodeBFSKplot = BFSK_decoded;
    end
    
end     % end for SNR loop


% error semilogy plot
figure(1);
p1 = semilogy (SNRdb,encoded_OOKerrorArr,'r-*'); hold on;
p2 = semilogy (SNRdb,OOKerrorArr,'k-*'); hold;hold;
p3 = semilogy (SNRdb,encoded_BPSKerrorArr,'g-*'); hold;hold;
p4 = semilogy (SNRdb,encoded_BFSKerrorArr,'b-*'); hold;
title('Bit Error Rate for different SNR');
legend([p1(1) p2(1) p3(1) p4(1)], {'Linear/OOK', 'Unencoded/OOK', 'Linear/BPSK', 'Linear/BFSK'}); ylabel('BER'); xlabel('SNR(dB)');
xlim([0 50]);

% data plot
figure (2);
subplot(211);
plot(dataPlot);
title('Data waveform');
ylim([-0.1 1.1]); xlim([0 length(data)]);

subplot(212);
plot(linearData);
title('Linear Encoded Data waveform');
ylim([-0.1 1.1]); xlim([0 length(linear_sig)]);

% modulated signals plot
figure (3);
subplot(411); plot(linearsigPlot); title('Linear Encoded Data signal'); ylim([-0.09 1.1]); xlim([0 1800]);
subplot(412); plot(modOOKplot); title('OOK modulated signal'); ylim([-2.1 2.1]); xlim([0 1800]); 
subplot(413); plot(modBPSKplot); title('BPSK modulated signal'); ylim([-2.1 2.1]); xlim([0 1800]); 
subplot(414); plot(modBFSKplot); title('BFSK modulated signal');  ylim([-2.1 2.1]); xlim([0 1800]);

% received signals plot
figure (4);
subplot(311); plot(recOOKplot); title(['OOK received signal at ' num2str(plotSNRdb) ' dB SNR']); xlim([0 1800]);
subplot(312); plot(recBPSKplot); title(['BPSK received signal at ' num2str(plotSNRdb) ' dB SNR']); xlim([0 1800]);
subplot(313); plot(recBFSKplot); title(['BFSK received signal at ' num2str(plotSNRdb) ' dB SNR']); xlim([0 1800]);

% demodulated signals plot
figure (5);
subplot (311); plot(demodOOKplot); title('OOK demodulated signal'); xlim([0 1800]); ylim([-0.09 1.1]);
subplot (312); plot(demodBPSKplot); title('BPSK demodulated signal'); xlim([0 1800]); ylim([-0.09 1.1]);
subplot (313); plot(demodBFSKplot); title('BFSK demodulated signal'); xlim([0 1800]); ylim([-0.09 1.1]);

% decoded signals plot
figure (6);
subplot (411); plot(dataPlot); title('Data waveform'); ylim([-0.09 1.1]); xlim([0 length(data)]);
subplot (412); plot(decodeOOKplot);  title ('OOK decoded signal'); ylim([-0.1 1.1]); xlim([0 length(data)]); ylim([-0.09 1.1]);
subplot (413); plot(decodeBPSKplot);  title ('BPSK decoded signal'); ylim([-0.1 1.1]); xlim([0 length(data)]); ylim([-0.09 1.1]);
subplot (414); plot(decodeBFSKplot);  title ('BFSK decoded signal'); ylim([-0.1 1.1]); xlim([0 length(data)]); ylim([-0.09 1.1]);


function [binary] = samplingANDdecision(x, samplePeriod, nBits, threshold)
    sampled = zeros(1, nBits);
    binary = zeros(1, nBits);
    for n = 1 : nBits
        sampled(n) = x((2 * n - 1) * samplePeriod / 2);
        if (sampled(n) > threshold)
            binary(n) = 1;
        else
            binary(n) = 0;
        end
    end
end