	function fitdata = gaussfit4(data,FWHMq,Offsetq,lb,ub)  

% This function requires the use of gaussian2.m
%
% Usage:    gaussfit4(data,FWHMq,Offsetq,lb,ub);  
% Where: data is a xy-matrix, FWHMq is the esimated FWHM  
% 	   Offsetq is the esimated offset (FLP) of the peak (X-axis)  
%  
% lb = lower bound of fitting parameters = 1x3 array [Ampl_lb FWHM_lb Offset_lb]
% ub = upper bound of fitting parameters = 1x3 array [Ampl_ub FWHM_ub Offset_ub]
%
% Output:   DataVector(1,3)  
% Where: DataVector(1) is the fitted Amplitude, A  
% 	     DataVector(2) is the fitted FHWM   
%	     DataVector(3) is the fitted Offset  
% 
%	Data: 10/4/99 
%	Written by DSL, modified by TPJK
%
%	The input parameters will be guestimated from data if
%	desired....  just don't put them in...
%

warning off;
ShowDataVar=1;

[Aq,Index]=max(data(:,2));  
Option=[0;5e-8;5e-8]; 
Option(14)=20000; 

if exist('Offsetq','var')==0 
	Offsetq=data(Index,1);
end

if exist('FWHMq','var')==0 
	TEMPDATA=data;
%	TEMPDATA(:,2)=normal(data(:,2));
	TEMPDATA(:,2)=(data(:,2))./max(data(:,2));
    
	%assuming the data is somewhat resolved....
	A=find(( TEMPDATA(:,2) > 0.3) & (TEMPDATA(:,2) < 0.7));
	A1=A(A<Index);
	A2=A(A>Index);

	if (isempty(A1) || isempty(A2))       
		FWHMq=20;
	else
		FWHMq=TEMPDATA(round(mean(A2)),1)-TEMPDATA(round(mean(A1)),1);
	end
end

start=[Aq,FWHMq,Offsetq];

if (~isempty(lb))&&(~isempty(ub))
    if min(ub-lb)==0
        ub=ub+.1; lb=lb-.1; lb(1:3)=abs(lb(1:3));
    end
end

fitdata = lsqcurvefit(@gaussian2,start,data(:,1),data(:,2),lb,ub);

% if ShowDataVar==1
%     Step=mean(diff(data(:,1)))/30; 
%     GaussX=min(data(:,1)):Step:max(data(:,1));  
%     plot(data(:,1),data(:,2),'r-');  
%     hold on;  
%     plot(GaussX,gaussian(GaussX,fitdata(1),fitdata(2),fitdata(3)),'b');  
%     title(['GAUSSFIT:  Width: ', num2str(fitdata(2)), '    Center: ', num2str(fitdata(3))])  
%     xlabel('X-axis');  
%     ylabel('Intensity');  
%     grid; 
%     hold off 
%     axis tight; 
% end

return 