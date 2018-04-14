function [] = TX_reverse_sentence(sound_name,time_bin,savename )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here





for r = 1:length(time_bin)
    
    [matrix,fs] = audioread(sound_name);
    
    wav_length = length(matrix);
    
    
    
    bin =floor( time_bin(r) / 1000 * fs);
    
    
    for i = 1:floor(wav_length / bin)
        start_point = (i - 1) * bin + 1;
        end_point = i * bin;
        matrix(end_point:-1:start_point) = matrix(start_point:end_point);
    end
    
    

    
    start_point = i * bin + 1;
    if length(matrix(start_point:end)) < bin
        matrix(end:-1:start_point) = matrix(start_point:end);
    end
    
   
    
    audiowrite(savename,matrix*0.5,fs)
end



end

