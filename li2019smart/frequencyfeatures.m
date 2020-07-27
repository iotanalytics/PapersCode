function [spec,centroid,spks,slogs] = frequencyfeatures(signalNor)
        fs = 500;
        amountpksF = 5;
        
        N = length(signalNor);
        fSignal = fft(signalNor);
        rSignal = real(fSignal);
        rSignal = rSignal(1:N/2+1);
        aSignal = abs(fSignal);
        aSignal = aSignal(1:N/2+1);
        
%       plot(real(fSignal));
%       figure;
%       plot(rSignal);
%       figure;
%       plot(aSignal);
        %peaks Real Signal
        [pks,locs] = findpeaks(rSignal);
        
        %getting the first 5 longest
        [spks, spks_order] = sort(pks,'descend');
        slogs = locs(spks_order,:);

%        amountpksF = min(amountpksF,length(spks));
        spks  = spks(1:amountpksF);
        slogs = slogs(1:amountpksF);
        slogs = sort(slogs);
        length(slogs);
        %Centroid
        T=abs(fSignal); % Take real parts of transform
        re = sum(T);
        N = size(T);
        N = N(1);
        acum = 0;
        for n = 1:N 
            f = n * fs / N;
            acum = acum + (f*T(n)); 
        end
        centroid = acum / re;
        
        %spectrum
        N = length(signalNor);
        psdx = (1/(fs*N)) * aSignal.^2;
        spec = sum(psdx);
end