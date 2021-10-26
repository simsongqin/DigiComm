function [generated_noise] = generate_noise(data_length, noise_power)
% this function generates an array of noise
% data_length: number of noise samples to be generated
% noise_power: noise variance of the signal
generated_noise = sqrt(noise_power/2) .* randn(1, data_length);
end

