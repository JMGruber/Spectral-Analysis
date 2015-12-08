readdir='D:\data\Werk\VU\data\FCP\Copy of FCPb Tjaart\spectra\pH5.5';
writedir='D:\data\Werk\VU\data\FCP\Copy of FCPb Tjaart\spectra\pH5.5\Analysis2';
for specnumber=1:108;
    mat=dlmread(fullfile(readdir,['spec' int2str(specnumber)]));
    h=image(mat(:,1),1:size(mat,2),mat(:,2:end)');
    axis([600 800 0.5 size(mat,2)+.5])
    xlabel('Wavelength (nm)');
    ylabel('Illumination time (s)');
    colormap(jet(min(round(max(max(mat(:,2:end)))*5/6),256)));
    saveas(h,fullfile(writedir,['contourplot' int2str(specnumber) '.jpg']));
end
                            
close all;