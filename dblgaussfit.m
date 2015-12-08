function fitdata = dblgaussfit(data,startq,lb,ub)

% Use together with dblgaussian2.m, and gaussian.m
%
% data = xy matrix
% startq = [FWHMq1,Offsetq1,FWHMq2,Offsetq2], corresponding to starting
%           values of the FWHM and the Offset (FLP) for the blue peak
%           (first 2) and the red peak (last 2).
%
% lb = lower bound of fitting parameters = 1x6 array [Ampl_lb1 FWHM_lb1 Offset_lb1 Ampl_lb2 FWHM_lb2 Offset_lb2]
% ub = upper bound of fitting parameters = 1x6 array [Ampl_ub1 FWHM_ub1 Offset_ub1 Ampl_lb2 FWHM_lb2 Offset_lb2]
%
% fitdata = [Amplitude1, FWHM1, Offset1, Amplitude2, FWHM2, Offset2];
%
% Written by TPJK

warning off;

wav=data(:,1);
x=data(:,2);
midwavx=(startq(2)+startq(4))/2;
midwav=find(wav>midwavx,1,'first');
[A1,Index1]=max(x(1:midwav));
[A2,Index2]=max(x(midwav:end));
if Index1==midwav
    A2=A2*0.7;
end
if Index2==midwav
    A1=A1*0.7;
end

start=[A1,startq(1:2),A2,startq(3:4)];

%fitdata = nlinfit(data(:,1),data(:,2),@dblgaussian2,[A1,FWHMq1,Offsetq1,A2,FWHMq2,Offsetq2]);
fitdata =lsqcurvefit(@dblgaussian2,start,data(:,1),data(:,2),lb,ub);
fitdata=real(fitdata);

%     Step=mean(diff(data(:,1)))/30; 
%     GaussX=min(data(:,1)):Step:max(data(:,1));  
%     plot(data(:,1),data(:,2),'-ro');  
%     hold on;  
%     plot(GaussX,dblgaussian(GaussX,fitdata(1),fitdata(2),fitdata(3),fitdata(4),fitdata(5),fitdata(6)),'b');  
%     title(['DOUBLEGAUSSFIT:  Width1: ', num2str(fitdata(2)), '    Center1: ', num2str(fitdata(3)),'  Width2: ', num2str(fitdata(5)), '    Center2: ', num2str(fitdata(6))])  
%     xlabel('Wavelength');  
%     ylabel('Intensity');  
%     grid; 
%     hold off 
%     axis tight; 

return