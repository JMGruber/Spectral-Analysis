complex=1;
readdir='D:\Dropbox\Michal\OCP\';
bgspec=1;   %true if a file containing a bg spectrum exists.
bgfilename1='Back led on_1.txt';    %is used when bgspec=true. File should be in readdir.
bgfilename2='Back spec_1.txt';    %is used when bgspec=true. File should be in readdir.
xmin=610; xmax=700; % maximum wavelength window considered 
subtrbl=1;
timebin=1;

   
mat=dlmread(fullfile(readdir,['spec' int2str(complex)]));

if timebin>1         
mat=spectimebin(mat,timebin); %bin along time
end

bgfile1=dlmread(fullfile(readdir,bgfilename1));
bgfile2=dlmread(fullfile(readdir,bgfilename2));
bgfile1=bgfile1(16:end,2);
bgfile1=bgfile1-sum(bgfile1)/length(bgfile1);
bgfile2=bgfile2(16:end,2);
bgfile2=bgfile2-sum(bgfile2)/length(bgfile2);
mat=mat(1:length(bgfile(:,1)),:);

    if subtrbl
        for spnr=2:length(mat(1,2:end))
        data=mat(:,spnr);
        baselineslope=(mean(data(length(data)-29:length(data)))-mean(data(1:30)))/length(data);
        linbaseline=mean(data(1:30))+baselineslope.*(1:length(data));
        data=data-linbaseline';
        mat(:,spnr)=data;
        int=sum(data);
        end
    end
    
nmin=find(mat(:,1)>=xmin,1,'first');
nmax=find(mat(:,1)<=xmax,1,'last');
mat=mat(nmin:nmax,:);
wav=mat(:,1);
turn=22;
mat(:,[1,2:turn])=bgsubtr(mat(:,[1,2:turn]),0,bgfile1(nmin:nmax));
mat(:,[1,turn+1:end])=bgsubtr(mat(:,[1,turn+1:end]),0,bgfile2(nmin:nmax));
%mat=bgsubtr(mat,0,bgfile(nmin:nmax));
%negative=find(mat<0);
% while ~isempty(negative)
     mat(:,2:end)=mat(:,2:end);
%     negative=find(mat<0);
% end
%mat(negative)=0;
size=length(mat(1,2:end));
sum1=[];
sum2=[];
centroid=[];

for i=2:size+1
centroid(i-1)=(sum(wav.*mat(:,i))/sum(mat(:,i)));
sum1(i-1)=sum(wav.*mat(:,i));
sum2(i-1)=sum(mat(:,i));
end
centroid;
plot(centroid)
ylim([630 710]) 