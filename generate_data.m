function [generated_data] = generate_data(data_length)
% this function generates an array of data
% data_length - number of data to be generated
generated_data = round(rand(1, data_length)) .* 2 - 1;
end
