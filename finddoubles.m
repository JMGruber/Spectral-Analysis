function [alldbl1 alldbl3] = finddoubles(alldbl,mat,singlespecfit,allSNR,noise,fwhmjump,firstgood,lastgood,SNRdbl,Qthr,fwhmdbl,broadthr,peak1,fwhm1,skewred)

% This function identifies double-band spectra in 3 ways:
% Part A - Relative threshold: Jumps in fitted widths > fwhmjump (indices contained in 1D matrix alldbl1)
% Part B - Absolute threshold for fitted skewness and width of single-band spectra (indices contained in 1D matrix alldbl2)
% Part C - Red, single-band spectral fits are checked for (small) blue bands (indices contained in 1D matrix alldbl3)
% All indices are combined into alldbl1

% Parts A & B operate only for spectra with SNR > SNRdbl, 
% while Part C looks for generally smaller peaks, with SNR > Qthr

% Warning: Make sure vib band is removed before this function is employed, since it often distorts a single-band fit!!


alldbl1=alldbl; alldbl2=[]; alldbl3=[]; %doubles identified by different algorithms are separated to facilitate debugging

%Part A: Relative threshold: jumps in widths from single-band fits
if lastgood>1
    allwidths=singlespecfit(2,:);
    difwidths=diff(allwidths(allwidths>0&(allSNR>SNRdbl)'));
    dif2widths=difwidths(1:end-1)+difwidths(2:end); %for transitory shifts
    fwhmjumpnew=fwhmjump; %because value can be changed
    if (~isempty(difwidths))&&((max([difwidths dif2widths])>fwhmjumpnew)||(abs(min([difwidths dif2widths]))>fwhmjumpnew))    %double bands are present
        dblsp=false;
        firstnormal=firstgood;
        while allwidths(firstnormal)<=0
            firstnormal=firstnormal+1;
        end
        prevwidth=allwidths(firstnormal);
        prev2width=prevwidth; %2 back, use for transitory shifts
        firstdbl=firstnormal;
        singles=[];
        meanswidth=0;
        for specnr=firstnormal+1:lastgood
            if (allwidths(specnr)>0)&&(allwidths(specnr)<broadthr)
                if (allwidths(specnr)-min(prevwidth,prev2width)>fwhmjumpnew)&&(~dblsp)&&(isempty(singles))&&(singlespecfit(3,specnr)<peak1+fwhm1)...
                     &&(allSNR(specnr)>SNRdbl)&&(singlespecfit(2,specnr)>fwhmdbl(specnr)) 
                    dblsp=true;
                    singles=firstdbl:specnr-1;
                    meanswidth=mean(allwidths(singles));
                    if length(singles)>1
                        stdswidth=std(allwidths(singles));
                        fwhmjumpnew=min(fwhmjumpnew,stdswidth*4);
                    end
                    firstdbl=specnr;
                elseif (allwidths(specnr)-max(prevwidth,prev2width)<-fwhmjumpnew)&&(dblsp)  %jump back to single
                    if ((meanswidth>0)&&(allwidths(specnr)-meanswidth<fwhmjumpnew))||...
                            (allwidths(specnr)-fwhm1<fwhmjumpnew) %if no previous single bands exist
                        alldbl1=[alldbl1 firstdbl:specnr-1];
                        dblsp=false;
                    end
                end
                prevwidth=allwidths(specnr);
                if allwidths(specnr-1)>0
                    prev2width=allwidths(specnr-1);
                end
            end
        end
        if dblsp
            alldbl1=[alldbl1 firstdbl:specnr];
        end
    end
end

for specnr=firstgood:lastgood
    if specnr==10
    end
    %Part B: Absolute threshold for skewness and width from single-band fits
    if (allSNR(specnr)>SNRdbl)&&(singlespecfit(1,specnr)>-.1)&&(singlespecfit(2,specnr)<broadthr)  
        if    (((singlespecfit(3,specnr)<peak1+fwhm1/4)&&(singlespecfit(4,specnr)>skewred-.1)&&(allSNR(specnr)>1.5*SNRdbl)) ||(singlespecfit(2,specnr)>=fwhmdbl(specnr)-.01))
            alldbl2=[alldbl2 specnr];
        end
    end
    
    %Part C: Checks if a red-shifted band contains an additional (small) blue band
    if (allSNR(specnr)>Qthr)&&(singlespecfit(1,specnr)>0)&&(singlespecfit(2,specnr)<broadthr) 
        if singlespecfit(3,specnr)>peak1+fwhm1/2
            bordermin=find(mat(:,1)>peak1-fwhm1,1,'first');
            bordermax=find(mat(:,1)>peak1+fwhm1,1,'first');
            bluepart=mat(bordermin:bordermax,specnr+1);
            [~,I]=max(bluepart);
            if (I>2)&&(I<length(bluepart)-2)
                SNR=(sum(bluepart(I-2:I+2))/5)/noise;
            else
                SNR=Qthr;
            end
            intblue=sum(bluepart);      %check if this really works!!
            if (SNR>Qthr)&&(intblue>200)     %large prob that this denotes the second, blue band
                alldbl3=[alldbl3 specnr];
            end
        end
    end
    
    %Part D: Identify doubles by relatively intense vib band - supposed to
    %be a subset of alldbl5!
%     if (allvibfit(1,specnr)>0)&&(allSNR(specnr)>SNRdbl)
%         if singlespecfit(1,specnr)/allvibfit(1,specnr)<mainvsvibfactor
%             alldbl4=[alldbl4 specnr];
%         end
%     end
end

if ~isempty(alldbl2)
    alldbl1=[alldbl1 alldbl2];
end
if ~isempty(alldbl3)
    alldbl1=[alldbl1 alldbl3];
end

% remove multiple entries
if ~isempty(alldbl1)
    alldbl1=sort(alldbl1);
    difalldbl=diff(alldbl1);
    alldbl1=[alldbl1(difalldbl>0) alldbl1(end)];
end