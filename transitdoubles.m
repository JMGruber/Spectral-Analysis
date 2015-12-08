function [trans,transi] = transitdoubles(alldbl,singlespecfit,allspecfit,peak1,fwhm1)
%{ 
 This function removes from alldbl all entries that correspond to
 (i) transitory doubles, i.e., an intermediate between a single- and double-band profile
 (ii) transitory singles, i.e., double-band spectra that arise from 
      a single-band spectrum shifting within the current time bin

USAGE:  alldbl = 1D array, containing the indices of double-band spectra
        singlespecfit = 4 x m matrix, containing single-band spectral fit parameters
        allspecfit = 8 x m matrix, containing double-band spectral fit parameters
        peak1 = estimated blue/single peak position
        fwhm1 = estimated fwhm of blue/single band

        trans = 1D array, containing entries of alldbl that correspond to transitions
        trans = 1D array, containing indices of alldbl that correspond to transitions

 ((i) and (ii) could be combined under one for loop to shorten execution time!)  
%} 
    
trans=[];
transi=[]; 
ti=0; % index in alldbl

% (i)
ti=0;
for spec=alldbl(1:end)
    ti=ti+1;
    if (spec>1)&&(spec<size(allspecfit,2))&&(singlespecfit(3,spec)<peak1+fwhm1)      % neglects strongly shifted bands
        if (((allspecfit(1,spec-1)>0)&&(allspecfit(1,spec+1)==-3)... % single to dbl
                 &&(singlespecfit(2,spec)<(singlespecfit(2,spec-1)+singlespecfit(2,spec+1))/2))...
             || ((allspecfit(1,spec+1)>0)&&(allspecfit(1,spec-1)==-3)... %dbl to single
                 &&(singlespecfit(2,spec)<(singlespecfit(2,spec+1)+singlespecfit(2,spec-1))/2)))

            trans=[trans spec];
            transi=[transi ti];
        end
    end
end

% (ii)
ti=0;
for spec=alldbl
    ti=ti+1;
    if (spec>1)&&(spec<size(allspecfit,2))
        if (allspecfit(5,spec-1)==-2)&&(allspecfit(5,spec+1)==-2)&&(abs(singlespecfit(3,spec+1)-singlespecfit(3,spec-1))>fwhm1/2) 
            trans=[trans spec];
            transi=[transi ti];
        end
    end
end