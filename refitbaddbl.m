function [allspecfit alldbl] = refitbaddbl(mat,allspecfit,alldbl,lbds,ubds,fitpar,fwhmdbl,constrainblue,SQthr,mainvsvibfactor,scalewidth_wav,scalef_wav)

% if constrainblue = true: use also single-band data to fit double-band spectra

peak2=fitpar(2); fwhm2=fitpar(3); skewmax=fitpar(4); %fwhmdbl=fitpar(5);

alldblfit=allspecfit(:,alldbl);
goodbool=alldblfit(7,:)>0;
gooddblspecs=alldbl(goodbool);
Qdbl=alldblfit(7,:)==0;
Qdblspecs=alldbl(Qdbl);
badbool=alldblfit(7,:)<0; %to exclude small int
baddblspecs=alldbl(badbool);
[~, allsingle]=find(allspecfit(5,:)==-2);
if ~isempty(allsingle)
    allsfit=allspecfit(1:4,allsingle);
    goodsingles=allsingle.*(allsfit(3,:)>0 & allsfit(3,:)<fitpar(1)+fitpar(3)/4);   %use single-band data
else
    goodsingles=[];
end
goodsingles=goodsingles(goodsingles>0);

if sum(goodbool)==0 %no good initial fits
    if ~isempty(goodsingles)
        peak1avg=mean(allspecfit(3,goodsingles),2);
        peak1std=std(allspecfit(3,goodsingles),0,2)/2; %make small, because normal std already tried in main script
        fwhm1avg=mean(allspecfit(2,goodsingles),2);
        fwhm1std=std(allspecfit(2,goodsingles),0,2)/2;
        skew1avg=mean(allspecfit(4,goodsingles),2)/3;  %restrict skewness of blue band!
        skew1std=std(allspecfit(4,goodsingles),0,2)/2;
    else
        peak1avg=lbds(3)+2;   %guess starting values
        fwhm1avg=lbds(2)+2;
        skew1avg=lbds(4)/3;
        peak1std=2;           %guess starting values
        fwhm1std=2;
        skew1std=.5;

    end
    if length(goodsingles)==1
        peak1std=2;           %guess starting values
        fwhm1std=2;
        skew1std=.5;
    end
    peak2avg=peak2;
    fwhm2avg=fwhm2;
    skew2avg=(lbds(8)+ubds(8))/4;
else
    if (constrainblue)&&(~isempty(goodsingles))
        goodspecs=[goodsingles gooddblspecs]; 
    else
        goodspecs=gooddblspecs;
    end
    peak1avg=mean(allspecfit(3,goodspecs),2);
    fwhm1avg=mean(allspecfit(2,goodspecs),2);
    skew1avg=mean(allspecfit(4,goodspecs),2)/3; %restrict skewness of blue band!
    if length(goodspecs)>1
        peak1std=std(allspecfit(3,goodspecs),0,2); %make small, because normal std already tried in main script
        fwhm1std=std(allspecfit(2,goodspecs),0,2);  
        skew1std=std(allspecfit(4,goodspecs),0,2);
    else
        peak1std=2;           %guess starting values
        fwhm1std=2;
        skew1std=.5;
    end
    peak2avg=mean(allspecfit(7,gooddblspecs),2);
    fwhm2avg=mean(allspecfit(6,gooddblspecs),2);
    skew2avg=mean(allspecfit(8,gooddblspecs),2);
end

lb=[SQthr fwhm1avg-fwhm1std peak1avg-peak1std*2 skew1avg-skew1std lbds(5:8)];
ub=[Inf fwhm1avg+fwhm1std peak1avg+peak1std*2 min(skew1avg+skew1std,skewmax) ubds(5) max(fwhm2avg+fwhm1std*5,ubds(6)) ubds(7:8)]; 
stillmisfits=[]; notmisfits=[];
countdbl=1;
for spec=baddblspecs
    if spec==9
    end
    data=mat(:,spec+1);
    dblspecfit=dblskewgaussfit2([mat(:,1) data],[fwhm1avg,peak1avg,skew1avg,fwhm2avg,peak2avg,skew2avg],lb,ub);
    nlb=lb;
    nub=ub;
    if (sum(dblspecfit([2,3,5,6,7])<lb([2,3,5,6,7])+.1)>0)||...
          (sum(dblspecfit([2,3,6,7])>ub([2,3,6,7])-.1)>0)||...
          (dblspecfit(4)<lb(4)+.01)||(dblspecfit(8)>ub(8)-.01)   %misfit - restrict red band too
       if sum(goodbool)>1
            peak2std=std(allspecfit(7,gooddblspecs),0,2);
            fwhm2std=std(allspecfit(6,gooddblspecs),0,2);
            skew2std=std(allspecfit(8,gooddblspecs),0,2);
            nlb(5:8)=[0 fwhm2avg-2*fwhm2std peak2avg-2*peak2std skew2avg-skew2std/2];           %red bands are restricted - give option/do later??
            nub(5:8)=[Inf fwhm2avg+2*fwhm2std peak2avg+2*peak2std min(skew2avg+skew2std/2,skewmax)]; %without /2??
            dblspecfit=dblskewgaussfit2([mat(:,1) data],[fwhm1avg,peak1avg,skew1avg,fwhm2avg,peak2avg,skew2avg],nlb,nub);
       end
       inc1=(nub(7)-nlb(3))/10;
       inc2=(nub(6)-nlb(6))/10;
       while (nlb(6)<nub(6)-3)&&(nlb(7)<nub(7)-1)&&...
            (((dblspecfit(2)<nlb(2)+.1)&&(dblspecfit(6)>nub(6)-.1))||...
                ((dblspecfit(7)<nlb(7)+.1)&&(dblspecfit(6)>nub(6)-.1))||...
                ((dblspecfit(2)<nlb(2)+.1)&&(dblspecfit(7)<nlb(7)+.1)))  %most likely strong overlap
            nlb(7)=nlb(7)+inc1;
            nlb(6)=nlb(6)+inc2;                
            dblspecfit=dblskewgaussfit2([mat(:,1) data],[fwhm1avg,peak1avg,skew1avg,fwhm2avg,peak2avg,skew2avg],nlb,nub);
        end
%         if (nlb(6)>=nub(6)-3)
%             nlb(7)=lb(7)+2;
%             dblspecfit=dblskewgaussfit2([mat(:,1) data],[fwhm1avg,peak1avg,skew1avg,fwhm2avg,peak2avg,skew2avg],nlb,nub);
%         elseif (nlb(7)>nub(7)-1)
%             nlb(6)=lb(6)+2;
%             dblspecfit=dblskewgaussfit2([mat(:,1) data],[fwhm1avg,peak1avg,skew1avg,fwhm2avg,peak2avg,skew2avg],nlb,nub);
%         end
        newfwhm1std=fwhm1std;
        newpeak1std=peak1std;
        newskew1std=skew1std;
        while (newfwhm1std>0.4)&&...
            ((sum(dblspecfit([2,3,6,7])<nlb([2,3,6,7])+.1)>0)||...
          (sum(dblspecfit([2,3,6,7])>nub([2,3,6,7])-.1)>0)||...
          (dblspecfit(4)<nlb(4)+.01)||(dblspecfit(8)>nub(8)-.01))   %misfit
            
            newfwhm1std=newfwhm1std/2;
            newpeak1std=newpeak1std/2;
            newskew1std=newskew1std/2;
            nlb(1:4)=[0 fwhm1avg-newfwhm1std peak1avg-newpeak1std skew1avg-newskew1std];
            nub(1:4)=[Inf fwhm1avg+newfwhm1std peak1avg+newpeak1std skew1avg+newskew1std]; 
            dblspecfit=dblskewgaussfit2([mat(:,1) data],[fwhm1avg,peak1avg,skew1avg,fwhm2avg,peak2avg,skew1avg],nlb,nub);
            
%            %fixed values: FWHMq1, Offsetq1, Skewq1 - doesn't work
%            dblspecfit=dblskewgaussfit([mat(:,1) data],dblspecfit(1),fwhm1avg,peak1avg,skew1avg,dblspecfit(2),fwhm2avg,peak2avg,skew2avg);
%            dblspecfit=[dblspecfit(1) fwhm1avg peak1avg skew1avg dblspecfit(2:5)];
        end
        % If an additional small red peak is fitted in the tail, fit with single peak only
        if (allspecfit(2,spec)<fwhmdbl(spec))...
            &&((dblspecfit(1)/dblspecfit(5)>5)||((dblspecfit(1)/dblspecfit(5)>mainvsvibfactor)&&(dblspecfit(7)-dblspecfit(3)<fwhm2))) % check if this is really necessary!!
            dblspecfit=[allspecfit(1:4,spec)' -2 -2 -2 -2];
            baddblspecs=trimarray(baddblspecs,countdbl);
            countdbl=countdbl-1;
            disp('Warning: increase fwhmdoubles or scalewidth_int!!');
        % As a last resort, use some estimation    
        elseif (sum(dblspecfit([2,3,6,7])<nlb([2,3,6,7])+.1)>0)...
                  ||(sum(dblspecfit([2,3,6,7])>nub([2,3,6,7])-.1)>0)...
                  ||(dblspecfit(4)<nlb(4)+.01)||(dblspecfit(8)>nub(8)-.01)   %misfit
            specfit1=skewgaussfit4([mat(:,1),data],[fwhm1avg,peak1avg,skew1avg],3,[],[]);
            ampl=specfit1(1);
            singlepeak1=skewgaussian3([max(dblspecfit(1),ampl*.75),fwhm1avg,peak1avg,skew1avg],mat(:,1));
            singlepeak2=data-singlepeak1;
            border=find(mat(:,1)>specfit1(3),1,'first');
            single=false;
            specfit2=skewgaussfit4([mat(border:end,1),singlepeak2(border:end)],[],3,nlb(5:8),nub(5:8));
            if (specfit2(3)>specfit1(3)+specfit1(2)+specfit2(2))&&(scalewidth_wav) %only for well separated peaks
                nub(6)=nub(6)+(specfit2(3)-specfit1(3))/scalef_wav; %scale fwhmdbl according to wavelength
                border=ceil((border+find(mat(:,1)>specfit2(3)-specfit2(2),1,'first'))/2);
                if specfit2(1)<SQthr+.1
                    specfit2=[gaussfit4([mat(border:end,1),singlepeak2(border:end)],specfit2(2),specfit2(3),[],[]) 0];
                    if specfit2(2)<nlb(6)+.1
                        single=true;
                    end
                else
                    specfit2=skewgaussfit4([mat(border:end,1),singlepeak2(border:end)],[],3,nlb(5:8),nub(5:8));
                end
            end
            if (~single)||(specfit2(2)<nlb(6)-.1)||(specfit2(1)<SQthr+.1)  
                specfit2=[gaussfit4([mat(border:end,1),singlepeak2(border:end)],specfit2(2),specfit2(3),nlb(5:7),nub(5:7)) 0]; % misfit for small signals - use constraints!
                if ((specfit2(2)<nlb(6)+.1)||(specfit2(3)<nlb(7)+.1)||(specfit2(3)>nub(7)-.1))...  % then this is likely a single-band spectrum
                   &&((allspecfit(2,spec)<fwhmdbl(spec))&&(allspecfit(4,spec)<ubds(4)/2))
                    dblspecfit=[allspecfit(1:4,spec)' -2 -2 -2 -2];
                    baddblspecs=trimarray(baddblspecs,countdbl);
                    countdbl=countdbl-1;
                    single=true;
                    disp('Warning: increase fwhmdoubles or scalewidth_int!!');
                end
            end
            if ~single
                n=.75;
                while (n<.9)&&((specfit2(2)>nub(6)-.1)||(specfit2(3)<nlb(7)+.1))   %most likely strong overlap
                    n=n+.04;
                    singlepeak1=skewgaussian3([max(dblspecfit(1),ampl*n),fwhm1avg,peak1avg,skew1avg],mat(:,1));
                    singlepeak2=data-singlepeak1;
                    specfit2=skewgaussfit4([mat(:,1),singlepeak2],[],3,nlb(5:8),nub(5:8));
                end
                singlepeak1=data-skewgaussian3([specfit2(1),specfit2(2),specfit2(3),specfit2(4)],mat(:,1));
                specfit1=skewgaussfit4([mat(:,1),singlepeak1],[],3,lb(1:4),ub(1:4));
                dblspecfit=[specfit1 specfit2];
                if (specfit2(1)<0)||(specfit2(2)>nub(6)-.1)||(specfit2(3)<nlb(7)+.1)
                     stillmisfits=[stillmisfits spec];
                else
                    notmisfits=[notmisfits spec];
                end
            end
        else 
            notmisfits=[notmisfits spec]; 
        end
    end
    allspecfit(:,spec)=dblspecfit;
    countdbl=countdbl+1;
end
if ~isempty(stillmisfits)
    notmisfits=[notmisfits gooddblspecs];
    if ~isempty(notmisfits)
        peak2avg=mean(allspecfit(7,notmisfits),2);
        fwhm2avg=mean(allspecfit(6,notmisfits),2);
        skew2avg=mean(allspecfit(8,notmisfits),2);
        if size(notmisfits,2)>1
            peak2std=std(allspecfit(7,notmisfits),0,2);
            fwhm2std=std(allspecfit(6,notmisfits),0,2);
            skew2std=std(allspecfit(8,notmisfits),0,2);
        else
            peak2std=2; %guess starting values
            fwhm2std=2;
            skew2std=.2;
        end
        lb(5:8)=[0 fwhm2avg-fwhm2std peak2avg-peak2std min(skew2avg-skew2std/2,0.01)];
        ub(5:8)=[Inf fwhm2avg+fwhm2std*2 peak2avg+peak2std min(max(skew2avg+skew2std,0.01),skewmax)]; 
        for spec=stillmisfits
           data=mat(:,spec+1);
           dblspecfit=dblskewgaussfit2([mat(:,1) data],[fwhm1avg,peak1avg,skew1avg,fwhm2avg,peak2avg,skew2avg],lb,ub);
           allspecfit(:,spec)=dblspecfit;
        end
    end
end
if isempty(baddblspecs)
    alldbl=gooddblspecs;
else
    alldbl=sort([gooddblspecs baddblspecs Qdblspecs]);
end