function newmatrix = bgsubtr(mat,subtrbl,bgspec,excwav)

%{
This function subtracts different backgrounds from data, particularly useful for a 
negative baseline or spectra with erroneous backgrounds

mat = the whole matrix of data, with the first column the wavelengths and the rest the different intensities
subtrbl = boolean, indicating whether a baseline should be subtracted
bgspec = 1D array of background data
        OR peak position of spectral fit (real) + width
        OR width
excwav = wavelength of excitation peak. Leave out when leakage of excitation is absent


Steps if no background spectrum is given
(1) A linear baseline is subtracted, defined by the average of the first 10 and last 10 data points

(2) if bgspec consists of two values (peak + width),
    (a) the whole matrix is subtracted by the spectrum with the 
        largest negative intensity (this applies only if a quenched state was accessed) 
    (b) if (a) does not apply, the average values of all the data in the range 
        peak position +- width are subtracted - this is the most obvious position where a negative peak will be

(3) if bgspec consists of one value (width)
    (a) the whole matrix is subtracted by the spectrum with the 
        largest negative intensity (this applies only if a quenched state was accessed) 
    (b) if (a) does not apply, the minimum intensity (minint) is found.
        If the average of all values within minint +- width is negative, these values are subtracted for all spectra

(4) If the intensity of any spectrum is still less than zero, the whole matrix is 
increased by a constant until int>0

TPJK, 2009-2012
%}


if (exist('subtrbl','var')==0)
    disp('erroneous baseline input');
elseif exist('bgspec','var')==0 
    disp('erroneous input');
elseif (length(bgspec)<=2)
    minint=-200; specbase=[];
    wav=mat(:,1);
    for spnr=2:size(mat,2)
        data=mat(:,spnr);
        if exist('excwav','var')~=0
            excnr=find(wav>excwav,1);
            deltawav=wav(2)-wav(1);
            if (sum(data(excnr-10:excnr+10))>20)
                [~,i]=max(data(excnr-10:excnr+10));
                data(excnr-10+i-round(6/deltawav):excnr-10+i+round(6/deltawav))=0;
            elseif (sum(data(excnr-10:excnr+10))<-20)
                [~,i]=min(data(excnr-10:excnr+10));
                data(excnr-10+i-round(5/deltawav):excnr-10+i+round(5/deltawav))=0;                
            end
        end
        if subtrbl
            baselineslope=(mean(data(length(data)-9:length(data)))-mean(data(1:10)))/length(data);
            linbaseline=mean(data(1:10))+baselineslope.*(1:length(data));
            data=data-linbaseline';
            mat(:,spnr)=data;
            int=sum(data);
            if (int<minint)
                  minint=int;
                  specbase=data;
                  minnr=spnr;
            end
        end
    end
    if isempty(specbase)
        if length(bgspec)==2
          bmin=find(mat(:,1)>bgspec(1)-bgspec(2),1);
          bmax=find(mat(:,1)>bgspec(1)+bgspec(2),1);
          band=mean(mat(bmin:bmax,2:size(mat,1)),2); 
          if sum(band)<0;
             specbase=[zeros(1,bmin-1) band zeros(1,size(mat,1)-bmax)]';
             minnr=0;
          end
        elseif length(bgspec)==1
            [~,spnr]=min(min(mat(:,2:size(mat,2))));
            [~,I]=min(mat(:,spnr+1));
            if (I>bgspec)&&(I<size(mat,1)-bgspec)
                band=mat(I-bgspec:I+bgspec,spnr+1);
                if sum(band)<-50 %should be modifyable???
                   specbase=[zeros(1,I-bgspec-1) band' zeros(1,size(mat,1)-(I+bgspec))]';
                   minnr=0;
               end
            end
        end
    end
    if ~isempty(specbase)
        for spnr=2:size(mat,2)
            mat(:,spnr)=mat(:,spnr)-specbase;
        end
        if (minnr>1)&&(minnr<size(mat,2))
            mat(:,minnr)=mat(:,minnr+1)+mat(:,minnr-1)./2;
        end;
    end
else
    if length(bgspec)~=length(mat)
        disp('erroneous background!');
    else
        for spnr=2:size(mat,2)
            mat(:,spnr)=mat(:,spnr)-bgspec;
        end
    end
    if subtrbl
        for spnr=2:size(mat,2)
            data=mat(:,spnr);
            baselineslope=(mean(data(length(data)-9:length(data)))-mean(data(1:10)))/length(data);
            linbaseline=mean(data(1:10))+baselineslope*(1:length(data));
            data=data-linbaseline';
            mat(:,spnr)=data;
        end
    end
end

minint=min(sum(mat(:,2:size(mat,2))));
% while minint<0
%     mat(:,2:size(mat,2))=mat(:,2:size(mat,2))+0.1;
%     minint=min(sum(mat(:,2:size(mat,2))));
% end
newmatrix=mat;