function array = trimarray(array,indices)

% removes elements corresponding to input indices from an array
% indices = array

for i=1:length(indices)
    if (indices(i)>1)&&(indices(i)<length(array)) %to disregard faulty indices
        array=[array(1:indices(i)-1) array(indices(i)+1:end)];
        indices=indices-1;
    elseif indices(i)==1
        if length(array)==1
            array=[];
        else
            array=array(2:end);
            indices=indices-1;
        end
    elseif indices(i)==length(array)
        array=array(1:end-1);
    end
end