function newmat = spectimebin(mat,timebin)

% This function bins consecutive spectra
% mat = [wavelength data]
% data = n x m matrix
% timebin = number of consecutive spectra to be binned

n=floor((size(mat,2)-1)/timebin);
mat2=zeros(length(mat),n);
for j=1:n
    k=mat(:,(j-1)*timebin+2:j*timebin+1); %1st column is wavelength
    mat2(:,j)=sum(k,2)./timebin;
end
newmat=[mat(:,1) mat2];