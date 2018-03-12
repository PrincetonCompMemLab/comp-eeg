function [filtered_signal] = stable_butterworth_bandpass(input_signal, cutoff_freqs, sample_rate, order)

low_freq = cutoff_freqs(1);
high_freq = cutoff_freqs(2);

h  = fdesign.bandpass('N,F3dB1,F3dB2', order, low_freq, high_freq, sample_rate);
Hd_constructed = design(h, 'butter');

% here we filter using a zero-phase filter function
filtered_signal = filtfilthd(Hd_constructed, input_signal);

end
