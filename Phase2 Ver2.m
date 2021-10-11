% Phase 2: Modulation for communication

clear all; close all; clc; 

% Var Declarations
fc = 10000;                                     % carrier freq
fs = 16*fc;                                     % sample freq
bbDataRate = 1000;                              % baseband DataRate 1kbps
nBits = 1024;                                   % data bits
SNRdb = -15:1:15;                               % SNR in db
SNR = (10.^(SNRdb/10));                         % convert SNR from db to dec
samplePeriod = fs/bbDataRate;

sigLen = fs*nBits/bbDataRate +1;                % signal length
amp = 2;                                        % amplitude
t = 0: 1/fs : nBits/bbDataRate;                 % time intervals
carrier = amp .* cos(2*pi*fc*t);
carrierHigh = amp .* cos(2*pi*5*fc*t);          % for BFSK 2nd higher freq carrier
[b,a] = butter(6,0.2);                          % low pass filter coefficients

% create arrays to store Bit Error Rates BER
OOKerrorArr = zeros(1,length(SNR)); aveOOKerror = 0;
BPSKerrorArr = zeros(1,length(SNR)); aveBPSKerror = 0;
BFSKerrorArr = zeros(1,length(SNR)); aveBFSKerror = 0;

runCycles = 5;

plotSNRdb = 5;                                        % for plotting at the specified SNR value 

for i = 1 : length(SNR)                                 % loop for diff SNR values
    totalOOKerror = 0; totalBPSKerror = 0; totalBFSKerror = 0;
    
    for j = 1 : runCycles                               % loop to cal ave error rate over number of runCycles
        % data Gen
        data = round(rand(1,nBits));
        % signal Gen
        sig = zeros(1,sigLen);
        for k = 1 : sigLen - 1
            sig(k) = data(ceil(k*bbDataRate/fs));
        end
        sig(sigLen) = sig(sigLen-1);
        
        %***** On-Off Keying (OOK) *****
            % modulation
        OOKsig = sig .* carrier;                        % if 1 = carrier; if 0 = 0;
            % noise Gen
        OOKsigPwr = bandpower(OOKsig);                  % sum of absolute square/length
        OOKnoisePwr = OOKsigPwr ./ SNR(i);
        OOKnoise = sqrt(OOKnoisePwr/2) .* randn(1,sigLen);
            % transmission
        OOKreceived = OOKsig + OOKnoise;
            % demodulation
        OOKdemod = OOKreceived .* (2.* carrier);        %coherent
        OOKdemod = filtfilt(b, a, OOKdemod);
            % sampling and decision logic
        OOK = samplingANDdecision(OOKdemod, samplePeriod, nBits, amp/2);
        
        %***** Binary Phase Shift Keying (BPSK) *****
            % modulation
        BPSKsig = sig .* 2 - 1;                         % if 1 -> 1*2-1 = 1; if 0 -> 0*2-1 = -1
        BPSKsig = BPSKsig .* carrier;
            % noise Gen
        BPSKsigPwr = bandpower(BPSKsig);
        BPSKnoisePwr = BPSKsigPwr ./ SNR(i);
        BPSKnoise = sqrt(BPSKnoisePwr/2) .* randn(1,sigLen);
            % transmission
        BPSKreceived = BPSKsig + BPSKnoise;
            % demodulation
        BPSKdemod = BPSKreceived .* (2.* carrier);
        BPSKdemod = filtfilt(b, a, BPSKdemod);
            % sampling and decision logic
        BPSK = samplingANDdecision(BPSKdemod, samplePeriod, nBits, 0);
        
        %***** Binary Frequency Shift Keying *****
            % modulation
        BFSKsigLow = carrier .* (sig == 0);
        BFSKsigHigh = carrierHigh .* (sig == 1);
        BFSKsig = BFSKsigLow + BFSKsigHigh;
            % noise Gen
        BFSKsigPwr = bandpower(BFSKsig);
        BFSKnoisePwr = BFSKsigPwr ./ SNR(i);
        BFSKnoise = sqrt(BFSKnoisePwr/2) .* randn(1,sigLen);
            % transmission
        BFSKreceived = BFSKsig + BFSKnoise;
            % demodulation
        BFSKdemodLow = BFSKreceived .* (2.* carrier);
        BFSKdemodLow = filtfilt(b, a, BFSKdemodLow);
        BFSKdemodHigh = BFSKreceived .* (2.* carrierHigh);
        BFSKdemodHigh = filtfilt(b, a, BFSKdemodHigh);
        BFSKdemod = BFSKdemodHigh - BFSKdemodLow;
            % sampling and decision logic
        BFSK = samplingANDdecision(BFSKdemod, samplePeriod, nBits, 0);
        
        % bit error
        OOKerror = 0; BPSKerror = 0; BFSKerror = 0;
        for l = 1 : nBits -1
            if (OOK(l) ~= data(l))
                OOKerror = OOKerror + 1;
            end
            if (BPSK(l) ~= data(l))
                BPSKerror = BPSKerror + 1;
            end
            if (BFSK(l) ~= data(l))
                BFSKerror = BFSKerror + 1;
            end
        end
        totalOOKerror = totalOOKerror + OOKerror;
        totalBPSKerror = totalBPSKerror + BPSKerror;
        totalBFSKerror = totalBFSKerror + BFSKerror;
 
    end     % end for runCycles loop
    
    aveOOKerror = totalOOKerror/runCycles;
    aveBPSKerror = totalBPSKerror/runCycles;
    aveBFSKerror = totalBFSKerror/runCycles;
    OOKerrorArr(i) = aveOOKerror/nBits;
    BPSKerrorArr(i) = aveBPSKerror/nBits;
    BFSKerrorArr(i) = aveBFSKerror/nBits;
    
    % store variable for specified SNR value for plotting
    if (SNRdb(i) == plotSNRdb)
        dataPlot = data;
        sigPlot = sig;
        modOOKplot = OOKsig;
        recOOKplot = OOKreceived;
        demodOOKplot = OOKdemod;
        decodeOOKplot = OOK;
        modBPSKplot = BPSKsig;
        recBPSKplot = BPSKreceived;
        demodBPSKplot = BPSKdemod;
        decodeBPSKplot = BPSK;
        modBFSKplot = BFSKsig;
        recBFSKplot = BFSKreceived;
        demodBFSKplot = BFSKdemod;
        decodeBFSKplot = BFSK;
    end
    
end     % end for SNR loop


% error semilogy plot
figure(1);
semilogy (SNRdb,OOKerrorArr,'-x'); hold on;
semilogy (SNRdb,BPSKerrorArr,'-*'); hold; hold;
semilogy (SNRdb,BFSKerrorArr,'-+'); hold;
title('Bit Error Rate for different SNR');
legend('OOK', 'BPSK', 'BFSK'); ylabel('BER'); xlabel('SNR(dB)');

% data plot
figure (2);
plot(dataPlot);
title('Data waveform');
ylim([-0.1 1.1]); xlim([0 length(data)]);

% modulated signals plot
figure (3);
subplot(411); plot(sigPlot); title('Data signal'); ylim([-0.1 1.1]); xlim([0 3000]);
subplot(412); plot(modOOKplot); title('OOK modulated signal'); xlim([0 3000]);
subplot(413); plot(modBPSKplot); title('BPSK modulated signal'); xlim([0 3000]);
subplot(414); plot(modBFSKplot); title('BFSK modulated signal'); xlim([0 3000]);

% received signals plot
figure (4);
subplot(311); plot(recOOKplot); title(['OOK received signal at ' num2str(plotSNRdb) ' dB SNR']); xlim([0 3000]);
subplot(312); plot(recBPSKplot); title(['BPSK received signal at ' num2str(plotSNRdb) ' dB SNR']); xlim([0 3000]);
subplot(313); plot(recBFSKplot); title(['BFSK received signal at ' num2str(plotSNRdb) ' dB SNR']); xlim([0 3000]);

% demodulated signals plot
figure (5);
subplot (311); plot(demodOOKplot); title('OOK demodulated signal'); xlim([0 3000]);
subplot (312); plot(demodBPSKplot); title('BPSK demodulated signal'); xlim([0 3000]);
subplot (313); plot(demodBFSKplot); title('BFSK demodulated signal'); xlim([0 3000]);

% decoded signals plot
figure (6);
subplot (411); plot(dataPlot); title('Data waveform'); ylim([-0.1 1.1]); xlim([0 length(data)]);
subplot (412); plot(decodeOOKplot);  title ('OOK decoded signal'); ylim([-0.1 1.1]); xlim([0 length(data)]);
subplot (413); plot(decodeBPSKplot);  title ('BPSK decoded signal'); ylim([-0.1 1.1]); xlim([0 length(data)]);
subplot (414); plot(decodeBFSKplot);  title ('BFSK decoded signal'); ylim([-0.1 1.1]); xlim([0 length(data)]);


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