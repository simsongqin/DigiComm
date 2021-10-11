function [generated_noise] = generate_noise(data_length, noise_power)
%This function generates an array of noise
%   Params:
%   data_length: number of noise samples
%   noise_power: noise variance of the signal
generated_noise = sqrt(noise_power/2) .* randn(1, data_length);
end

