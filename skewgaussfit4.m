function fitdata = skewgaussfit4(data,startq,fntype,lb,ub)  

% Use function together with skewgaussian1.m/2.m/3.m (depending on fntype),
% and gaussfit3.m

% data = xy matrix
% startq = 1x3 array with starting values of FWHM, Offset (FLP), and Skewness
%          Leave empty if the values should be determined b.m.o. a Gaussian fit.

% fntype = type of skewed Gaussian used
% fntype=3: FS69 function (Frazer and Suzuki 1969)
% fntype=2: Mikas's function (more flexible)
% fntype=1: Composite function

% lb = lower bound of fitting parameters = 1x4 array [Ampl_lb FWHM_lb Offset_lb Skew_lb]
% ub = upper bound of fitting parameters = 1x4 array [Ampl_ub FWHM_ub Offset_ub Skew_ub]
% Leave lb and/or ub empty if you don't want to constrain the parameters

% fitdata = [Amplitude, FWHM, Offset, Skewness];

% Written by TPJK



warning off;

if isempty(startq)
    paraest=gaussfit3(data);
    start=[paraest(1) paraest(2) paraest(3) 0.01];
else
    start=[max(data(:,2)),startq];
    if start(4)==0
        start(4)=0.01;
    end
end

if (~isempty(lb))&&(~isempty(ub))
    if min(ub-lb)==0
        ub=ub+.1; lb=lb-.1; lb(1:3)=abs(lb(1:3));
    end
end

switch fntype
    case 1
        fitdata = lsqcurvefit(@skewgaussian1,start,data(:,1),data(:,2),lb,ub);
    case 2
        fitdata = lsqcurvefit(@skewgaussian2,start,data(:,1),data(:,2),lb,ub);
    case 3
        fitdata = lsqcurvefit(@skewgaussian3,start,data(:,1),data(:,2),lb,ub);
end

% Step=mean(diff(data(:,1)))/30; 
% GaussX=min(data(:,1)):Step:max(data(:,1));  
% plot(data(:,1),data(:,2),'-r');  
% hold on;  
% plot(GaussX,skewgaussian1([fitdata(1),fitdata(2),fitdata(3),fitdata(4)],GaussX),'b');  
% title(['SkewGAUSSFIT:  Width: ', num2str(fitdata(2)), '    Center: ',...
%     num2str(fitdata(3)),'   Skewness:  ',num2str(fitdata(4)) ]);  
% xlabel('X-axis');  
% ylabel('Intensity');  
% grid; 
% hold off 
% axis tight; 