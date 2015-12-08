    goodspecs1=allcint1>0 & allcpeaks1>0 & (allcpeaks2<=0 | allcpeaks2>allcpeaks1+fwhm1/2);
    goodspecs2=allcint2>0 & allcpeaks2>allcpeaks1+fwhm1/2;

    % plot histograms
    % (1) Distribution of avg blue and red peak positions per spectral sequence
    % (2) Std of (1)
    % (3) Distribution of all peak positions
    % (4) Spectral jumps
    histpeakmin=floor(min([avgpeaks1(avgpeaks1>0) avgpeaks2(avgpeaks2>0)])/binFLP)*binFLP;
    histpeakmax=ceil(max([avgpeaks1(avgpeaks1>0) avgpeaks2(avgpeaks2>0)])/binFLP)*binFLP;
    markavgpeak=histpeakmin:binFLP:histpeakmax;
    histpeaks1=hist(avgpeaks1(avgpeaks1>0),markavgpeak);
    histpeaks2=hist(avgpeaks2(avgpeaks2>0),markavgpeak);
    h1=figure;
    subplot(2,2,1); bar(markavgpeak,histpeaks1,'b');
    xlabel('avg FLP (nm)'); ylabel('# spectra');
    hold on; bar(markavgpeak,histpeaks2,'r');
    axis tight;

    markstdpeak=0:binFLP:ceil(max([stdpeaks1(stdpeaks1>0) stdpeaks2(stdpeaks2>0)]));
    histstdpeaks1=hist(stdpeaks1(stdpeaks1>0),markstdpeak);
    histstdpeaks2=hist(stdpeaks2(stdpeaks2>0),markstdpeak);
    if ~isempty(histstdpeaks1)
        subplot(2,2,2); bar(markstdpeak, histstdpeaks1,'b');
        xlabel('avg SDFP (nm)'); ylabel('# spectra');
    end
    if ~isempty(histstdpeaks2)
        hold on; bar(markstdpeak, histstdpeaks2,'r');
    end
    axis tight;
        
    %markers=xmin:binFLP:xmax;
    histpeakmin=floor(min([allcpeaks1(goodspecs1>0);allcpeaks2(goodspecs2>0)])/binFLP)*binFLP;
    histpeakmax=ceil(max([allcpeaks1(goodspecs1>0);allcpeaks2(goodspecs2>0)])/binFLP)*binFLP;
    markallpeaks=histpeakmin:binFLP:histpeakmax;
    histallpeaks1=hist(allcpeaks1(goodspecs1>0),markallpeaks);
    histallpeaks2=hist(allcpeaks2(goodspecs2>0),markallpeaks);
    subplot(2,2,3); bar(markallpeaks,histallpeaks1,'b'); 
    xlabel('all FLP (nm)'); ylabel('# spectra');
    hold on; bar(markallpeaks, histallpeaks2,'r'); 
     axis tight;

    %Three types of jumps:
    %1: blue band 
    %2: red band
    %3: difference between blue and red bands - green
    minjump=floor(min(min([allcjumps1 allcjumps12],[],2))/binJumps)*binJumps;
    maxjump=ceil(max(max([allcjumps1 allcjumps2 allcjumps12],[],2))/binJumps)*binJumps;  
    markjumps=minjump:binJumps:maxjump;
    jumps1=hist(allcjumps1(abs(allcjumps1)>jumpthr),markjumps);
    if ~isempty(jumps1)
        subplot(2,2,4); bar(markjumps,jumps1,'b');
        xlabel('Jump size (nm)'); ylabel('Occurrence');
        jumps2=hist(allcjumps2(abs(allcjumps2)>jumpthr),markjumps);
        if ~isempty(jumps2)
            hold on; bar(markjumps,jumps2,'r');
        end
        jumps12=hist(allcjumps12(abs(allcjumps12)>jumpthr),markjumps);
        if ~isempty(jumps12)
            hold on; bar(markjumps,jumps12,'g');
        end
        axis tight;
    end
    
    % Plot blue yield vs red wavelength - if no relationship, complex jumped back and forth between a single-band blue and red state!!
    h2=figure; plot(allcpeaks2(goodspecs2>0 & allcfi1>0),allcfi1(goodspecs2>0 & allcfi1>0),'o');
    xlabel('Red FLP (nm)'); ylabel('Blue yield (%)');

    % Plot fit parameters
    % colour code:
    % red: redder peak of double-band spectra
    % blue: bluer peak of double-band spectra
    % green: single bands of complexes exhibiting double bands
    % black: single-band spectra
    h3=figure; 
    %goodspecs1=allpeaks1>peakmin+.1 & allfwhm1<fwhmmax2;
    %goodspecs2=allpeaks2>peakmin+fwhmmin & allpeaks2<peakmax & allfwhm2>fwhmmin & allskew2<skewmax2-.01 & allskew2>skewmin1+.01 & allint2>SNRdbl/2;
    subplot(2,2,1); 
    %plot(allpeaks1(allpeaks2==-2&goodspecs1),allint1(allpeaks2==-2&goodspecs1),'k.');
    %hold on; plot(allpeaks1(goodspecs2),allint1(goodspecs2),'b.'); 
    %plot(allpeaks2(goodspecs2),allint2(goodspecs2),'r.');
    alldoubles=doubles(doubles>0); 
    doubles_bool=zeros(maxspecs,matrixsize);
    doubles_bool(:,alldoubles)=1; % %can use doubles instead of alldoubles before zeros added at end of script
    doubles_bool=doubles_bool>0;
    doubles_singles=allcpeaks2==-2 & goodspecs1 & doubles_bool;
    singles_singles=allcpeaks2==-2 & goodspecs1 & ~doubles_bool;
    %allreds=(avgpeaks1>0).*doubles_bool;
    plot(allcpeaks1(singles_singles),allcint1(singles_singles),'k.');
    hold on; plot(allcpeaks1(doubles_singles),allcint1(doubles_singles),'g.');
    plot(allcpeaks1(goodspecs2),allcint1(goodspecs2),'b.'); 
    plot(allcpeaks2(goodspecs2),allcint2(goodspecs2),'r.');
    xlabel('FLP (nm)'); ylabel('intensity (cps)');
    subplot(2,2,2); 
    plot(allcpeaks1(allcpeaks2==-2&goodspecs1),allcfwhm1(allcpeaks2==-2&goodspecs1),'k.');
    hold on; plot(allcpeaks2(goodspecs2),allcfwhm2(goodspecs2),'r.');
    plot(allcpeaks1(goodspecs2),allcfwhm1(goodspecs2),'b.'); 
    %plot(allpeaks2(allpeaks1==-2&allpeaks2>0),allfwhm2(allpeaks1==-2&allpeaks2>0),'k.');
    xlabel('FLP (nm)'); ylabel('fwhm (nm)');
    subplot(2,2,3); 
    plot(allcpeaks1(allcpeaks2==-2&goodspecs1),allcskew1(allcpeaks2==-2&goodspecs1),'k.');
    hold on; plot(allcpeaks2(goodspecs2),allcskew2(goodspecs2),'r.');
    plot(allcpeaks1(goodspecs2),allcskew1(goodspecs2),'b.'); 
    %plot(allpeaks2(allpeaks1==-2&allpeaks2>0),allskew2(allpeaks1==-2&allpeaks2>0),'k.');
    xlabel('FLP (nm)'); ylabel('skewness');
    subplot(2,2,4); 
    plot(allcfwhm1(allcpeaks2==-2&goodspecs1),allcint1(allcpeaks2==-2&goodspecs1),'k.');
    hold on; plot(allcfwhm2(goodspecs2),allcint2(goodspecs2),'r.');
    plot(allcfwhm1(goodspecs2),allcint1(goodspecs2),'b.'); 
    %plot(allfwhm2(allpeaks1==-2&allpeaks2>0),allint2(allpeaks1==-2&allpeaks2>0),'k.');
    xlabel('fwhm (nm)'); ylabel('intensity (cps)');

    % Plot fit parameters of vib band
    if fitvib
        h4=figure;
        goodvibs=(allcvibampl>0 & goodspecs1); 
        subplot(2,2,1); plot(allcpeaks1(goodvibs),allcvibampl(goodvibs),'.');
        xlabel('FLP blue band (nm)'); ylabel('Ampl vib band (nm)');
        subplot(2,2,2); plot(allcpeaks1(goodvibs),allcvibfwhm(goodvibs),'.');
        xlabel('FLP blue band (nm)'); ylabel('FWHM vib band (nm)');
        subplot(2,2,3); plot(allcpeaks1(goodvibs),allcvibpeaks(goodvibs),'.');
        xlabel('FLP blue band (nm)'); ylabel('FLP vib band (nm)');
        subplot(2,2,4); plot(allcpeakjumps4vib(allcvibjumps>5),allcvibjumps(allcvibjumps>5),'.'); % arbitrary number to avoid small jumps that result from fitting inaccuracies
        xlabel('FLP jumps blue band (nm)'); ylabel('FLP jumps vib band (nm)');
    end

    % Plot intensities of singles and different bands of doubles
    h5=figure; subplot(1,2,1);
    plot(allcpeaks1(singles_singles),allcint1(singles_singles),'k.');
    hold on; plot(allcpeaks1(doubles_singles),allcint1(doubles_singles),'g.');
    plot(allcpeaks2(goodspecs2),allcint1(goodspecs2)+allcint2(goodspecs2),'m.');
    xlabel('FLP (nm)'); ylabel('total intensity (cps)');
    avgintss=mean(allcint1(singles_singles));
    avgintds=mean(allcint1(doubles_singles));
    avgintdd=mean(allcint1(goodspecs2)+allcint2(goodspecs2));
    avgintb=mean(allcint1(goodspecs2));
    avgintr=mean(allcint2(goodspecs2));
    subplot(1,2,2); 
    bar([avgintss,0,0,0,0],'k');
    hold on; bar([0,avgintds,0,0,0],'g');
    bar([0,0,avgintdd,0,0],'m'); 
    bar([0,0,0,avgintb,0],'b'); 
    bar([0,0,0,0,avgintr],'r');
    xlabel('ss    ds    dd    db    dr'); ylabel('avg. intensity (cps)');