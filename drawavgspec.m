readdir='D:\Sheela\Mutant 3 5C POC flushed e-5';
writedir='D:\Sheela\Mutant 3 5C POC flushed e-5';
%bg=dlmread(fullfile(readdir,'bgall.txt'));
mat_avg=zeros(249,1);
includefiles=[20 30 32 47 48 61 68 79 80 87 105 115 125 143 161 184 1059 1061 1072 1073]
% WT includefiles=[79 80 126 142 160 164 169 168 200 231]

Meanval=zeros(249,1);
for specnumber=5:941;
j=1;
    while (j<=length(includefiles))&&(includefiles(j)<=specnumber)
        if specnumber==includefiles(j)
            mat=dlmread(fullfile(readdir,['spec' int2str(specnumber)]));
            Meanval=(Meanval+mean(mat(:,2:size(mat,2)),2))/2;
        end
        j=j+1;
    end
    


end
j



for specnumber=5:941;
    mat=dlmread(fullfile(readdir,['spec' int2str(specnumber)]));
    
    mat_2=zeros(249,1);
    for j=2:31;
        if sum(mat(:,j)-Meanval)>300;
        if sum(mat(1:53,j)-Meanval(1:53))<300;
        mat_2=mat_2+mat(:,j)-Meanval;
        end
        %mat_2=smooth(mat_2,3);
        end
    end
    plot(mat(:,1),mat_2)
    mat_avg=mat_avg+mat_2;
    
    
   
 
   % h=image(mat(:,1),1:size(mat,2),mat(:,2:end)');
   % axis([600 800 0.5 size(mat,2)+.5])
   % xlabel('Wavelength (nm)');
   % ylabel('Illumination time (s)');
   % colormap(jet(min(round(max(max(mat(:,2:end)))*5/6),256)));
   % saveas(h,fullfile(writedir,['contourplot' int2str(specnumber) '.jpg']));
end
%mat_avg=smooth(mat_avg,15);
plot(mat(:,1),mat_avg)
                            
   