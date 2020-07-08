function [stds,entropyS,pksre,maxPeak,locMaxPeak,before,after] = timefeatures(signalNor)
      %max amount of peaks to get for the signal
      amountpks = 5;
      amoutdatafrommaxpeak = 10;

      %standart deviation
      stds = std(signalNor);
      
      %Signal Normalization to get the entropy
      x_normalized = signalNor/max(abs(signalNor));
      entropyS = entropy(x_normalized);
      
      %finding and getting peaks
      pks = findpeaks(signalNor);
      
      %getting the maximun peaks
      maxPeak = max(pks);
      pks = sort( pks ,'descend');
      pksres = zeros(5,1);
      numpks = length(pks);
      amountp = min(numpks,amountpks);
      pksre = pks(1:amountp);
      pksre = [pksre ; zeros(amountpks - length(pksre), 1)];
      
      %max peaks location
      locMaxPeak = find(signalNor==maxPeak);
      signalNor(locMaxPeak);
      
      %partial data before and after max peak
      from = max(1,locMaxPeak-(amoutdatafrommaxpeak));                    %este
      before = signalNor(from:locMaxPeak-1);
      if( length(before) < amoutdatafrommaxpeak)
          before = [ zeros(amoutdatafrommaxpeak - length(before), 1) ; before];
      end
      mincopy = min(locMaxPeak+(amoutdatafrommaxpeak),length(signalNor)); %este
      after  = signalNor(locMaxPeak+1:mincopy);
      if( length(after) < amoutdatafrommaxpeak)
          after = [ after ; zeros(amoutdatafrommaxpeak - length(after), 1)];
      end
      
end