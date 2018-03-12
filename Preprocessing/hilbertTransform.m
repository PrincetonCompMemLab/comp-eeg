function hilpow = hilbertTransform(EEG, bpfind)
chans = EEG.nbchan;
hilpow = nan(chans, length(bpfind)-1,EEG.pnts);
% hilpha = nan(chans, length(bpfind)-1,EEG.pnts);

for chan = 1:chans;
    for bpfreq = 1:length(bpfind)-1
        bpfs = stable_butterworth_bandpass(double(EEG.data(chan, :)), ...
            [bpfind(bpfreq) bpfind(bpfreq+1)], EEG.srate, 4);
        
        h = hilbert(bpfs);
        hilpow(chan, bpfreq, :) = abs(h);
        %         hilpha(chan, bpfreq, :) = angle(h);
    end
end
end