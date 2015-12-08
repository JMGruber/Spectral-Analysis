function [alldbl4,alldbl5,allvibfit] = findsmalldbl(mat,alldbl,allvibfit,singlespecfit,allsingle,meanvibs,meansinglespecs,peak1,fwhm1,skewred,intfarred,farredminwav,allSNRfraction)

% This function identifies remaining double-band spectra after subtracting
% an estimated vibrational band from previously identified single-band spectra
% alldbl5 = 1D matrix of indices, corresponding to a typically far-red bands
% alldbl4 = 1D matrix of indices, corresponding to double bands similar to
%           previously identified ones (contained in alldbl)

alldbl5=[]; alldbl4=[];
for spec=allsingle
    if singlespecfit(3,spec)<peak1+fwhm1/2
        data=mat(:,spec+1);
        if (~isempty(meanvibs))
            index=find(mat(:,1)>=meansinglespecs(3),1,'first');
            amplblue=sum(data(index-1:index+1))/3;
            vibpeakampl=amplblue*meanvibs(1)/meansinglespecs(1);
            vibband=gaussian(mat(:,1),vibpeakampl,meanvibs(2),meanvibs(3)); 
            data=data-vibband;
            %after subtracting average vib wing, check to see if another (small) peak is present
            bordermin=find(mat(:,1)>peak1+2*fwhm1,1,'first');
            bordermin2=find(mat(:,1)>farredminwav,1,'first');
            intred=sum(data(bordermin:end));
            intred2=sum(data(bordermin2:end));  %only for case of negative baseline
%                                if (SNR>SNRfarred)&&(intred>intfarred)
            if (intred>intfarred)||(intred2>intfarred)%(SNR>SNRfarred) %large prob that this denotes an additional small red band
                alldbl5=[alldbl5 spec];
                allvibfit(:,spec)=[vibpeakampl,meanvibs(2),meanvibs(3)];
            end
        elseif ~isempty(alldbl)
            if (singlespecfit(4,spec)>skewred/2)&&(allSNRfraction(spec)>1)   %perhaps some small ones will also be fitted
                alldbl4=[alldbl4 spec];
            end
        end
    end
end