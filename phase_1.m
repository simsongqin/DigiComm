clear all; close all; clc;

%declare variables
data_length = 1024;
signal_power = 1;
type = 'power';

SNR_dB = 0:1:50;

SNR = convert_dB_to_dec(SNR_dB, type);

noise_powers = signal_power ./ SNR;

%generate our 1024 data bits
data = generate_data(data_length);

threshold = 0;
test_samples = 20;
bit_errors = [];

for i = 1 : length(SNR)
    bit_errors(i) = 0;
    for j = 1 : test_samples
        noise = generate_noise(data_length, noise_powers(i));
        
        received_signal = data + noise;
        
        received_signal = 2*(received_signal >= threshold)-1;
        
        error_signal = (received_signal~=data);
        
        bit_errors(i) = bit_errors(i) + mean(error_signal);
    end
    bit_errors(i) = bit_errors(i)/test_samples;
end

%print noise mean, SD and variance
disp(mean(noise))
disp(std(noise))
disp(var(noise))

% calculate theoretical BER
theory_rate = (1 / 2) * erfc(sqrt(SNR));
    
%plot the result           
figure(1)
semilogy (SNR_dB, theory_rate,'r', 'linewidth', 1.5);
ylabel('BER');
xlabel('SNR (dB)')
title('BER vs SNR (dB) - Step Size: 1');
hold on
semilogy (SNR_dB, bit_errors,'bx', 'linewidth', 2);
axis([0 50 1/(test_samples*data_length) 1]);
legend('Theoretical BER','Real BER');
hold off

%plot the result           
figure(2)
semilogy (SNR_dB(1:5:50), theory_rate(1:5:50),'r', 'linewidth', 1.5);
ylabel('BER');
xlabel('SNR (dB)')
title('BER vs SNR (dB) - Step Size: 5');
hold on
semilogy (SNR_dB(1:5:50), bit_errors(1:5:50),'bx', 'linewidth', 2);
legend('Theoretical BER','Real BER');
axis([0 50 0 1]);
hold off

%data generation
figure(3) 
subplot(311)
plot(data);
xlim([0 1024]);
title('Generated Data (1024 Bits)')
%noise generation
subplot(312) 
plot(noise);
xlim([0 1024]);
title('Generated Noise (1024 Bits)')
%received data generation
subplot(313) 
plot(received_signal);
xlim([0 1024]);
title('Received Data (1024 Bits)')