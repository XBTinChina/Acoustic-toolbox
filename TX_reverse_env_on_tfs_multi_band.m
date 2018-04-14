function [reverse_fine] = TX_reverse_env_on_tfs_multi_band(x1,env_binsize, env_Fco,fine_Fco Fs)
% MULTI_BAND_CHIMERA - Synthesize pair of multi-band "auditory chimeras"
% by dividing each signal into frequency bands, then interchanging
% envelope and fine structure in each band using Hilbert transforms.
%
% Usage:  [e1_fs2, e2_fs1] = multi_band_chimera(x1, x2, Fco, Fs, refilter)
%	x1, x2	original signals
%	Fco	    band cutoff frequencies (Nbands+1), or filter bank
%	Fs	    sampling rate in Hz (default 1)
%   refilter  set to 1 to filter again after exchange operation (default 0)
% 	e1_fs2  chimera with envelope of x1, fine structure of x2 in each band
% 	e2_fs1  chimera with envelope of x2, fine structure of x1 in each band
%
%	Copyright Bertrand Delgutte, 1999-2000
%
if nargin < 3, error('Specify original signals and cutoff frequencies'); end
if nargin < 4, Fs = 1; end
if nargin < 5, refilter = 0; end

if min(size(x1)) == 1, x1 = x1(:); end	% make column vectors


nchan = size(x1, 2);


env_bin = floor( env_binsize / 1000 * Fs);

% Because the Hilbert transform and the filter bank are both linear operations,
% they commute and associate.  We create a bank of complex FIR filters whose
% real and imaginary parts are in quadrature (cos and sin).  This complex filter
% is directly applied the original signals. The Hilbert envelope in each band
% is the absolute value of the complex filter output.
% This approach avoids computation of large FFTs as in Matlab's 'hilbert'.

if min(size(fine_Fco)) > 1 | isstr(Fs),
    fine_b = fine_Fco;	% kluge to specify filters
else fine_b = quad_filt_bank_without_hilbert(fine_Fco, Fs);
end


if min(size(env_Fco)) > 1 | isstr(Fs),
    env_b = env_Fco;	% kluge to specify filters
else env_b = quad_filt_bank(env_Fco, Fs);
end


reverse_fine = zeros(size(x1));
reverse_env = zeros(size(x1));
fine_structure = zeros(size(x1));
% loop over filters in bank



for k = 1:size(fine_b,2)
    
    zfilt1 = fftfilt(fine_b(:,k), x1);
    
    temp_fs = cos(angle(zfilt1));
    
    fine_structure = fine_structure + temp_fs;
    
end




for k = 1:size(env_b,2),
    
    zfilt1 = fftfilt(env_b(:,k), x1);
    
    
    temp_fs = cos(angle(zfilt1));
    
    fine_structure = fine_structure + temp_fs;
    
    for i = 1:floor(length(x1) / env_bin)
        start_point = (i - 1) * env_bin + 1;
        end_point = i * env_bin;
        temp_env(end_point:-1:start_point) = temp_env(start_point:end_point);
    end
    
    
    if   end_point ~= length(x1)
        
        
        start_point = i * env_bin + 1;
        if length(x1(start_point:end)) < env_bin
            
            
            temp_env(end:-1:start_point) = temp_env(start_point:end);
            
            
        end
        
        
    end
    
    
    zfilt2 = fftfilt(fine_b(:,k), x1);
    
    
    
    % interchange envelope and fine structure
    temp_reverse_fs_env = temp_env .* temp_fs;
    % refilter backward to avoid delay accumulation
    
    
    % accumulate over frequency bands
    reverse_fine = reverse_fine + temp_reverse_fs_env;
    
    
    
end
