function [ output_args ] = TX_shuffle_sentence(sound_name,local_bin,global_bin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


[matrix fs] = wavread(sound_name);

wav_length = length(matrix);

local_bin = local_bin / 1000 * fs;

if exist('global_bin')
    global_bin = global_bin / 1000 * fs;
else
    global_bin = wav_length;
end


ceil( wav_length / global_bin )


global_order = randperm( ceil( wav_length / global_bin ));

local_order = zeros(length(global_order,global_bin/local_bin));

for i = 1:length(global_order)
    local_order(i,:) = randperm(global_bin/local_bin);
end



for i = 1:floor(wav_length / bin)
    start_point = (i - 1) * bin + 1;
    end_point = i * bin;
    matrix(end_point:-1:start_point) = matrix(start_point:end_point);
end

start_point = i * bin + 1;
matrix(end:-1:start_point) = matrix(start_point:end);

wavwrite(matrix,fs,['sentence_' num2str(time_bin(r)) ])


end

