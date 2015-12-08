function [singlespecfit,dataminvib] = guessvibband(mat,allvibfit,singlespecfit,dataminvib,badsinglesnr,fwhmdbl,meanvibs,meansinglespecs,newlbs,newubs,noise,mainvsvibfactor,fwhmvib,peakvib)

% This function estimates the vibrational band of misfitted ones.
% If there are no previous good fits, the user-given fit parameters are used

for specnr=badsinglesnr
    data=mat(:,specnr+1);
    newubs(2)=fwhmdbl(specnr);
    if (~isempty(meanvibs))&&(~isempty(meansinglespecs))
        index=find(mat(:,1)>=meansinglespecs(3),1,'first');
        amplblue=sum(data(index-1:index+1))/3;
        vibpeakampl=amplblue*meanvibs(1)/meansinglespecs(1);
        allvibfit(:,specnr)=[vibpeakampl,meanvibs(2),meanvibs(3)]; %approximation
        data=data-gaussian(mat(:,1),vibpeakampl,meanvibs(2),meanvibs(3)); 
        specfit=skewgaussfit4([mat(:,1),data],[],3,newlbs,newubs); %fittype=3
    else % estimate vib fit parameters, based on user-given values
        maxdata=max(data);
        amplblue=maxdata-1.5*(sqrt(maxdata)+noise);
        vibpeakampl=amplblue/(mainvsvibfactor*1.5);
        allvibfit(:,specnr)=[vibpeakampl,fwhmvib,peakvib];
        data=data-gaussian(mat(:,1),vibpeakampl,fwhmvib,peakvib); 
        specfit=skewgaussfit4([mat(:,1),data],[],3,newlbs,newubs); %fittype=3
    end
    if (specfit(1)<newlbs(1)+.1)||(specfit(2)>newubs(2)-.1) %misfit
        specfit=[-1 -1 -1 -1];
    end
    singlespecfit(:,specnr)=specfit;
    dataminvib(:,specnr+1)=data;
end