function [alldbl Qdbl] = addQdbl(alldbl,SNRfraction)

% This function identifies doubles having too small intensities to be resolved
% and add them to alldbl

% USAGE:
% SNRfraction = allSNR/SNRdbl, where allSNR is a 1D array of SNR values for
%               the spectral sequence under consideration
% alldbl = 1D array, containing the indices of double-band spectra
% Qdbl = 1D array, containing the indices of the identified small doubles

difalldbl=diff(alldbl);
gaps=find(difalldbl>1);
Qdbl=[];
for j=gaps
    for k=alldbl(j)+1:alldbl(j+1)-1
        if SNRfraction(k)<1
            Qdbl=[Qdbl k];
        end
    end
end
if ~isempty(Qdbl)
    alldbl=sort([alldbl Qdbl]);
end