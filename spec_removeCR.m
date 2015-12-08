function newmat=spec_removeCR(mat,threshold,numavg)

% Remove outliers, usually due to Cosmic Rays, b.m.o. linear interpolation 
% Compares the data with a set threshold
% USAGE:
%       mat = full matrix, i.e., [wav spec1 spec2 ...]
%       threshold = intensity/pixel for an outlier
%       numavg = number of data points on each side of outlier to interpolate
%          it is set large (3-4) if SNR is low and/or spectral resolution is high

if ~isempty(find(mat(:,2:end)>threshold,1)) %slightly faster than sum(sum(mat(:,2:end)>threshold))>0 
    newdata=zeros(size(mat));
    for j=2:size(mat,2)
        data=mat(:,j);
        if sum(data>threshold)>0 %or: max(data)>threshold (which is faster?)
            px=find(data>threshold);
            %if cosmic ray hit too close to edge of window, make values zero
            if (px(1)<numavg)||(px(end)>length(data)-numavg)
                newval=zeros(length(px),1);
            else  %interpolate by using mean of numavg values before and after
                if data(px(1)-1)-data(px(1)-2)>threshold/2 %smearing to adjacent pixels
                    data(px(1)-1)=data(px(1)-2);
                elseif data(px(end)+1)-data(px(end)+2)>threshold/2 %smearing to adjacent pixels
                    data(px(1)+1)=data(px(1)+2); 
                end
                startslope=mean(data(px(1)-numavg:px(1)-1)); 
                endslope=mean(data(px(end)+1:px(end)+numavg));
                slope=(endslope-startslope)/length(px);
                newval=startslope+slope*(1:length(px));
            end
            data(px)=newval;
        end
        newdata(:,j)=data;
    end
    newmat=[mat(:,1) newdata(:,2:end)];
else
    newmat=mat;
end