function array = remvalarray(array,values)

% This function determines whether any values exist in "array" and removes them
% array, values = 1D array
% Important condition: both input arrays should be sorted!!

indices=[];
for i=1:length(array)
    j=1;
    valexists=false;
    while (j<=length(values))&&(values(j)<=array(i))&&(~valexists)
        if array(i)==values(j)
            valexists=true;
        end
        j=j+1;
    end
    if valexists
        indices=[indices i];
    end
end
array=trimarray(array,indices);