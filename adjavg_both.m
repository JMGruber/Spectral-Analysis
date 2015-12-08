function avgdata=adjavg_both(data,n)

% This function averages adjacent data points to smoothen the data
% n = number of points to average ON BOTH SIDES of every data point
newdata=data;
if n>=1
    for i=1:length(data)
        j=n;
        while (i<=j)||(i>length(data)-j)
            j=j-1;
        end
        if j==0
            newdata(i)=data(i);
        else
            newdata(i)=mean(data(i-j:i+j));
        end;
    end
end
avgdata=newdata;