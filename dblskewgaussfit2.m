function fitdata = dblskewgaussfit2(data,startq,lb,ub)

% Use function together with dblgaussfit.m, dblskewgaussian2.m, skewgaussian3.m, 
%                           gaussfit3.m, and dblgaussian2.m
% data = xy matrix
% startq = [FWHMq1,Offsetq1,Skewq1,FWHMq2,Offsetq2,Skewq2],
%           corresponding to starting values of FWHM, Offset (FLP), and Skewness
%           of blue peak (q1s) and red peak (q2s)
% These values are improved b.m.o. the output from a double-gaussian fit
%
% lb = lower bound of fitting parameters = 1x8 array [Ampl_lb1 FWHM_lb1 Offset_lb1 Skew_lb1 Ampl_lb2 FWHM_lb2 Offset_lb2 Skew_lb2]
% ub = upper bound of fitting parameters = 1x8 array [Ampl_ub1 FWHM_ub1 Offset_ub1 Skew_ub1 Ampl_ub2 FWHM_ub2 Offset_ub2 Skew_ub2]
%
% startq can be estimated from lb and ub, and vice versa. Just leave the
% corresponding arrays empty (but not all of them simultaneously!)
%
% fitdata = [Amplitude1, FWHM1, Offset1, Skewness1, Amplitude2, FWHM2, Offset2, Skewness2];
%
% Written by TPJK
warning off;

if isempty(startq)
    if isempty([lb ub])
        disp('Warning: Lacking starting values or boundary conditions!!');
        fitdata=zeros(1,8);
        return 
    else
        startq=(ub+lb)/2;
        startq=[startq(2:4) startq(6:8)];
        if min(ub-lb)==0
            ub=ub+.15; lb=lb-.1; lb(1:3)=abs(lb(1:3));
        end
    end
else
    if startq(3)==0
        startq(3)=0.01;
    end
    if startq(6)==0
        startq(6)=0.01;
    end
    if isempty([lb ub])
        lb=[0 startq(1)-sqrt(startq(1)) startq(2)-startq(1) startq(3)-.1 0 startq(4)-sqrt(startq(4)) startq(5)-startq(4) startq(6)-.1];
        ub=[Inf startq(1)+sqrt(startq(1)) startq(2)+startq(1) startq(3)+.1 Inf startq(4)+sqrt(startq(4)) startq(5)+startq(4) startq(6)+.1];
    elseif isempty(lb)
        lb=ub-abs(startq);
    elseif isempty(ub)
        ub=lb+abs(startq);
    else
        if min(ub-lb)==0
            ub=ub+.15; lb=lb-.1; lb(1:3)=abs(lb(1:3));
        end
        if (min(startq-[lb(2:4) lb(6:8)])<0)||(min([ub(2:4) ub(6:8)]-startq)<0)
            startq=(ub+lb)/2;
            startq=[startq(2:4) startq(6:8)];
        end
    end
end
    
lbd=[lb(1:3) lb(5:7)];
ubd=[ub(1:3) ub(5:7)];

fitdataq=dblgaussfit(data,[startq(1:2) startq(4:5)],lbd,ubd);
% A1=fitdataq(1);
% FWHM1=fitdataq(2);
% Offset1=fitdataq(3);
% A2=fitdataq(4);
% FWHM2=fitdataq(5);
% Offset2=fitdataq(6);

%start=[A1,FWHM1,Offset1,Skewq1,A2,FWHM2,Offset2,Skewq2];
start=[fitdataq(1:3) startq(3) fitdataq(4:6) startq(6)];
start(6)=max(start(6),startq(4)); %FWHM2

fitdata =lsqcurvefit(@dblskewgaussian2,start,data(:,1),data(:,2),lb,ub);
fitdata=real(fitdata);

% Step=mean(diff(data(:,1)))/30; 
% GaussX=min(data(:,1)):Step:max(data(:,1));  
% plot(data(:,1),data(:,2),'-ro');  
% hold on;  
% plot(GaussX,dblskewgaussian2([fitdata(1),fitdata(2),fitdata(3),fitdata(4),fitdata(5),fitdata(6),fitdata(7),fitdata(8)],GaussX'),'b');  
% plot(GaussX,skewgaussian3([fitdata(1),fitdata(2),fitdata(3),fitdata(4)],GaussX'),':k');
% plot(GaussX,skewgaussian3([fitdata(5),fitdata(6),fitdata(7),fitdata(8)],GaussX'),':k');
% title(['DOUBLEDKEWGAUSSFIT:  Width1: ', num2str(fitdata(2)), '    Center1: ', num2str(fitdata(3)),'  Skew1: ', num2str(fitdata(4)),'  Width2: ', num2str(fitdata(6)), '    Center2: ', num2str(fitdata(7)),'  Skew2: ', num2str(fitdata(8))])  
% %    title(['DOUBLEDKEWGAUSSFIT:   Center1: ', num2str(fitdata(2)),'  Skew1: ', num2str(fitdata(3)),'  Width2: ', num2str(fitdata(5)), '    Center2: ', num2str(fitdata(6)),'  Skew1: ', num2str(fitdata(7))])  
% xlabel('Wavelength');  
% ylabel('Intensity');  
% grid; 
% hold off 
% axis tight; 

return 