function alldbl = removetransitdoubles(alldbl,singlespecfit,allspecfit,peak1,fwhm1)

% This function removes from alldbl all entries that correspond to
% (i) quneched spectra (too small intensities) - remove this option, otherwise chance to miss far reds!!
% (ii) transitory doubles, i.e., an intermediate between a single- and
% double-band profile
% (iii) transitory singles, i.e., double-band spectra that arise from 
%  a single-band spectrum shifting within the current time bin


trans=[];
transi=[]; 
ti=0; % index in alldbl

% (i)
% for spec=alldbl
%     ti=ti+1;
%     if singlespecfit(1,spec)==0
%         trans=[trans spec];
%         transi=[transi ti];
%     end
% end

if length(alldbl)>2
    ti=1;
    for spec=alldbl(2:end-1)
        ti=ti+1;
        if ((allspecfit(5,spec-1)==-2)&&(allspecfit(5,spec+1)==-2)&&(abs(singlespecfit(3,spec+1)-singlespecfit(3,spec-1))>fwhm1/2))...  % (iii) %check this!!
            || ((singlespecfit(3,spec)<peak1+fwhm1)... % neglects strongly shifted bands      (ii)
               &&   (((allspecfit(1,spec-1)>0)&&(allspecfit(1,spec+1)==-3)... % single to dbl
                     &&(singlespecfit(2,spec)<(singlespecfit(2,spec-1)+singlespecfit(2,spec+1))/2))...
                  ||((allspecfit(1,spec+1)>0)&&(allspecfit(1,spec-1)==-3)... %dbl to single
                     &&(singlespecfit(2,spec)<(singlespecfit(2,spec+1)+singlespecfit(2,spec-1))/2))))
                    
           trans=[trans spec];
           transi=[transi ti];
        end
    end
end

for spec=1:length(trans)
    if trans(spec)<alldbl(end)
        alldbl(transi(spec)-spec+1:end-1)=alldbl(transi(spec)-spec+2:end);
    end
    alldbl=alldbl(1:end-1);
end