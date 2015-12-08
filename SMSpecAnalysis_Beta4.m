% The parameters can be set as to fit maximally 2 spectral bands + 1 vibrational wing,
% the latter of which corresponds to the "blue" spectral band
% Each spectral sequence is categorised in one of 4 or 5 groups, i.e.,
% deads, doubles, broads, goods (and optionally blues) - see output file "sortedpeaks"
%
% The user is given the option to let the script automatically
% - identify spectra exhibiting spectral bluing
% - subtract a background spectrum and/or a baseline
% - identify double-band spectra
% - resolve a vibrational band
%
% For a single spectral sequence, the indices of the spectra within a spectral
% sequence that correspond to double-band spectra can be given manually.
%
% USAGE: Just run, but don't forget to change the necessary parameters first
%        (currently tuned for LHCII trimers)
%
%        Use together with:
%        - skewgaussfit4.m, skewgaussian1.m, skewgaussian2.m, skewgaussian3.m
%        - dblgaussfit.m, dblskewgaussfit2.m, dblskewgaussian2.m, dblgaussian2.m
%        - gaussfit3.m, gaussfit4.m, mingauss2.m, gaussian.m, gaussian2.m
%        - bgsubtr.m, spec_removeCR.m
%        - spectimebin.m, adjavg_both.m
%        - guessvibband.m
%        - finddoubles.m, transitdoubles.m, addQdbl.m, findsmalldbl.m, refitbaddbl.m
%        - trimarray.m, remvalarray.m
%
% OUTPUT 
% - Fit parameters (amplitudes, peak positions, fwhm, skewness),
% intensities, SNR, and RMS of single- and double-band spectra and vibrational bands
% - Averages of spectral fit parameters and spectral jumps
% - Many secondary parameters
% 
% MAIN CHALLENGES
% 1. Accurately identify and resolve double-band spectra, especially for
%  - strongly overlapping bands
%  - small red bands
% 2. Differentiate between a vibrational band and an additional red band
% 3. Differentiate between double-band and broadened, quenched single-band states
%
% ALGORITHM
% For every complex:
% 1. The average spectrum of all the unquenched states is fitted. 
%    Based on the spectral profile, excessively blue-shifted and broad spectra are identified. 
% 2. Each spectrum is fitted with a single skewed Gaussian. Based on the profile,
%    blues, broads, double, and red-shifted spectra are identified and for
%    the remaining profiles, the corresponding vibrational bands are resolved.
% 3. Good fits to single-band spectra and vibrational bands are identified. 
%    Means and fitting boundaries are calculated.
% 4. A vibrational band for all misfitted ones is estimated, based on the good 
%    fits or, alternatively, user-given values
% 5. Double-band spectra are identified by means of three different algorithms,
%    apart from those already found
% 6. Misfitted singles and vibrational bands are refitted, based on fit parameters of good fits,
%    except those identified as doubles.
% 7. All single-band spectra are scanned for possible double bands after the
%    well estimated vibrational bands are subtracted from the single-band profiles.
% 8. Double-band spectra are fitted.
% 9. Misfitted doubles are refitted, based on fitting parameters of good fits.
% 10. Average fitting parameters and jumps are calculated
%
% Methods to identify double-band states
% 1. Significantly broadened and/or positively skewed spectra after
%    subtraction of vibrational band or its whole spectral region (ALGORITHM Step 2).
% 2. Relative threshold in fitted widths of single-band spectra
% 3. Absolute threshold for fitted skewness and width of single-band spectra
% 4. Presence of (small) blue bands in addition to identified single, red bands
%    (2-4: ALGORITHM Step 5, invoking external function "finddoubles.m")
% 5. Misfitted single-band spectra after refitting these and the vibrational bands
%    (ALGORITHM Step 6).
% 6. Presence of (small) red bands in addition to identified single, blue bands,
%    after subtracting the vibrational band (ALGORITHM Step 7, "findsmalldbl.m").
% 7. Significantly positively skewed single-band profiles after subtracting
%    vibrational band (ALGORITHM Step 7, "findsmalldbl.m").
%
% Some additional remarks
% - The two spectral bands are fitted with skewed Gaussians, while that of the
%  (small) vibrational band with a normal Gaussian.
% - A double-band specturm that is fitted with a single skewed Gaussian is
%   often the result of a too much contrained width of the blue/red band
% - Negative skewness doesn't make much sense for far reds, but it's often
%   incorporated to improve the effectiveness/accuracy of the fit.
% - To save execution time, it's (almost) impossible to find all doubles
%  (alldbl) by means of one or two external functions.
% - The accuracy of the regressions of course decreases with the noise
%
% Proposed improvements
% - Exponential Gaussian in some cases
% - Routine for resolving double bands, especially by differentiating
%   between well separable (easily resolvable) and strongly overlapping bands.
% - Possible identification of distinct (special), user-defined peaks
% - Average of Q spectra fitted. (For that, Qint will be increased.)
% - Correlation of intensity traces with spectral sequences.
%
% Output of fit results:
% -3: double-band spectrum (supposed to have been resolved)
% -2: single peak
% -1: misfit
% 0: intensity too small
%
% Tjaart Krüger
% Vrije Universiteit Amsterdam
% 2009 - 2012

tic; 
clear all; close all;

%% Initialisation

% Input files
allfiles=50; %A:B or matrix of individual spec files
skipfiles=[];    %numbers don't need to be in any particular order; matrix can also be empty
includefiles=[];    %good complexes that shouldn't be skipped, but are identified as blue or broad, e.g. due to an impurity at some point
dbl2norm=[];   %force identified double-banded complex to be single-banded
readdir='D:\Dropbox\Michal\test\';
writedir='D:\Dropbox\Michal\test\analysis';

% Background subtraction
bgspec=0;           %true if a file containing a bg spectrum exists. Useful if background intensity < 0. Otherwise background subtraction is done automatically.
bgfilename='bg50.txt';    %is used when bgspec=true. File should be in readdir.
bgi=false;              %true if there's a bg file (bg1, bg2, ...) corresponding to each spectral file (spec1, spec2, ...). False if not used or only one bg file should be used.
subtrbl=0;          %true if a baseline should be subtracted from all spectra
subtrpeak=false;         %true if the background sometimes contains a real spectrum

% Intensity and SNR thresholds
CRthr=450;    %threshold of individual intensity values. Larger values are definitely due to cosmic rays.
Qthr=2;       %max SNR of Q spectra. If SNR is below this threshold, chances are good for a misfit of single-band spectra. Can be as small as 1.5 - 2, since Q states are not ocnsidered yet.
SNRdead=2.5;    %max SNR when complex is considered "dead" (if large SNR, set ~3.5, otherwise ~2.5)
SNRbulk=50;     %min SNR of clearly >1 complex

testbluing=true;    %true if spectral bluing (excessively blue-shifted spectra) should be identified
constrainblue=true; %if true, then the mean of good fits to single bands or blue parts is used to improve the fits of single-band spectra, vib wings, and double-band spectra for strong band overlaps
                    %fit parameters of red-shifted states are not used.

fittype=3; %type of skewed Gaussian used for fitting
% fittype=3: FS69 function (Frazer and Suzuki 1969)
% fittype=2: Mikas's function (more flexible)
% fittype=1: Composite function

% Double-band spectra
alldbl=[];          %manual input for indices of double-band spectra for a SINGLE spectral sequence. Leave empty if it should be determined automatically.
finddbls=1;  %checks for double bands. If false, then only necessary parameters need to be defined.
SNRdbl=5;           %SNR above which the script considers a double band
fwhmjump=8;         %jump in fwhm to consider double bands (8 is approx max for Lhcbs)
fwhmdouble=25;      %max fwhm of a single peak before the script considers double bands (~25 for Lhcbs, 26 for minors).
rmsthr=10;
findfarreds=0;   %typically small bands
farredminwav=730;   %minimum wavelength of far-red bands
intfarred=50;      %minimum intensity of a far-red band.
%skewblue=-0.1; %max skewness of a single peak before the script considers double bands to the blue - not incorporated yet
skewred=1.2; %same as skewblue, but for red tails

scalewidth_wav=true; %scale fwhmdouble according to peak wavelength
scalef_wav=4;        %if scalewidth_wav = true, fwhmdouble = fwhmdouble+(fittedpeakmax-peak1)/scalef_wav
scalewidth_int=true; %scale fwhmdouble according to intensity
scalef_int=16;       %if scalewidth_int = true, fwhmdouble = fwhmdouble+(scalef_int/SNR)^2 (minors: 15)

%vibrational band
fitvib=true;            %true if an additional small red spectrum should be fitted
peakvib=730; fwhmvib=50; %typical peak position and fwhm of vib band
peakvibmin=720; peakvibmax=750; 
fwhmvibmax=75;
minvib=715;             %lower boundary of vibrational band (nm)
mainvsvibfactor=4;      %minimum ratio of main (blue) peak vs. vibrational band, before the latter is likely a real red band

% first estimations for fitting. If finddbls=false, values for fwhm2 and peak2 will be ignored
fwhm1=20; fwhm2=30;
skewness=0.03;
peak1=680; peak2=700; %typical peak positions for double-band spectra;
%specialpeakmin=679; specialpeakmax=681; %special category of peaks - not incorporated yet

% Boundary conditions for fitting
xmin=620; xmax=825; % maximum wavelength window considered for all calculations. Set e.g. to 0 and 1000, resp., if full range is required.
peakmin=672; peakmax=820; %note if < peakmin, complex prob underwent spectral bluing
peak2min=0; % lower boundary of red peak. Set to 0 if it should be determined automatically
fwhmmin=12; fwhmmax1=25; fwhmmax2=80; % 1 and 2 refer to blue and red peaks, respectively. Better to choose smaller max values and scale with wavelength
skewmin1=-.3; skewmin2=-.3; skewmax1=skewred; skewmax2=skewred;

% Binning: input data
wavbin=0;   % wavelength binning: number of adjacent datapoints in both directions
wavbinsmall=1; % binning for small peaks (see SNRsmall)
SNRsmall=2*Qthr; % when SNR < SNRsmall, larger wavelength binning is employed
timebin=0;  % number of consecutive spectra to be binned

% Binning: output data
binFLP=1; %bin size of FL peak position (in nm)
binJumps=1; %bin size of jump sizes (in nm)
jumpthr=2; %jump size threshold (in nm) - two consecutive spectral peaks that differ by more than this number are considered a jump

maxspecs = 45;  % max number of spectra per complex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           User may change parameter values until here                 %%%                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preallocate matrix sizes

nfiles=size(allfiles);
allfiles=sort(allfiles);
matrixsize=allfiles(end);
allcpeaks1=zeros(maxspecs,matrixsize);
allcpeaks2=zeros(maxspecs,matrixsize);
allcfwhm1=zeros(maxspecs,matrixsize);
allcfwhm2=zeros(maxspecs,matrixsize);
allcskew1=zeros(maxspecs,matrixsize);
allcskew2=zeros(maxspecs,matrixsize);
allcfi1=zeros(maxspecs,matrixsize);
allcint1=zeros(maxspecs,matrixsize);
allcint2=zeros(maxspecs,matrixsize);
allcSNR=zeros(maxspecs,matrixsize);
allcjumps1=zeros(maxspecs,matrixsize);
allcjumps2=zeros(maxspecs,matrixsize);
allcjumps12=zeros(maxspecs,matrixsize); %difference between peaks1 and peaks2
allcvibampl=zeros(maxspecs,matrixsize);
allcvibpeaks=zeros(maxspecs,matrixsize);
allcvibfwhm=zeros(maxspecs,matrixsize);
allcvibjumps=zeros(maxspecs,matrixsize);
allcpeakjumps4vib=zeros(maxspecs,matrixsize); %jumps in peak1 to correlate with jumps in vib band
allcrms=zeros(maxspecs,matrixsize);

avgint1=zeros(1,matrixsize); avgfwhm1=zeros(1,matrixsize); avgskew1=zeros(1,matrixsize); avgpeaks1=zeros(1,matrixsize); 
avgint2=zeros(1,matrixsize); avgfwhm2=zeros(1,matrixsize); avgskew2=zeros(1,matrixsize); avgpeaks2=zeros(1,matrixsize);
stdint1=zeros(1,matrixsize); stdfwhm1=zeros(1,matrixsize); stdskew1=zeros(1,matrixsize); stdpeaks1=zeros(1,matrixsize); 
stdint2=zeros(1,matrixsize); stdfwhm2=zeros(1,matrixsize); stdskew2=zeros(1,matrixsize); stdpeaks2=zeros(1,matrixsize);

deads=zeros(1,matrixsize); blues=zeros(1,matrixsize); doubles=zeros(1,matrixsize); broads=zeros(1,matrixsize); goods=zeros(1,matrixsize); % 5 spectral categories

%% Initial settings

skipfiles=sort(skipfiles);
includefiles=sort(includefiles);
dbl2norm=sort(dbl2norm);

if ~isempty(alldbl)
    if (nfiles-numel(skipfiles))>1
        disp('Warning: manual input for doubles applies to one spectral sequence only!');
        break
    else
        manualdbl=true;
    end
else
    manualdbl=false;
end

if (fitvib)&&(Qthr>SNRdbl)  %for now, Qthr <= SNRdbl, to distinguish double-band fits
    Qthr=SNRdbl;
end

for i=allfiles
    complex=i;
    j=1;
    goodfile=true; %i.e. not to be skipped
    while (j<=length(skipfiles))&&(skipfiles(j)<=complex)&&(goodfile)
        if complex==skipfiles(j)
            goodfile=false;
        end
        j=j+1;
    end
    if goodfile
        j=1;
        incfile=false;
        while (j<=length(includefiles))&&(includefiles(j)<=complex)&&(~incfile)
            if complex==includefiles(j)
                incfile=true;
            end
            j=j+1;
        end
        
        if (isempty(dbl2norm))||(~finddbls)
            finddbltemp=finddbls;
        else
            j=1;
            finddbltemp=true;
            while (j<=length(dbl2norm))&&(dbl2norm(j)<=complex)&&(finddbltemp)
                if complex==dbl2norm(j)
                    finddbltemp=false;
                end
                j=j+1;
            end
        end
        
        mat=dlmread(fullfile(readdir,['spec' int2str(complex)]));
%        mat=dlmread(fullfile(readdir,['spec (' int2str(complex) ')']));
        if timebin>1         
            mat=spectimebin(mat,timebin); %bin along time
        end
        mat=spec_removeCR(mat,CRthr,3); %correct for cosmic-ray impacts
        nmin=find(mat(:,1)>=xmin,1,'first');
        nmax=find(mat(:,1)<=xmax,1,'last');
        mat=mat(nmin:nmax,:);
        wav=mat(:,1);
        peakmx=min(peakmax,wav(end-5));
        if bgspec   %background subtraction
            if bgi
                bgfile=dlmread(fullfile(readdir,['bg' int2str(complex)]));
            else
                bgfile=dlmread(fullfile(readdir,bgfilename));
            end
            mat=bgsubtr(mat,subtrbl,bgfile(nmin:nmax));    
        else
            if subtrpeak
                mat=bgsubtr(mat,subtrbl,fwhm1);
            else
                mat=bgsubtr(mat,subtrbl); %change that baseline is subtracted!!
            end
        end
        data=mat(:,2:end);
        
        if ~manualdbl
            alldbl=[];
        end
        
        bgregionmax=max(peak1-(fwhmmax1*1.5),xmin+20);
        bgregionmaxx=find(wav>=bgregionmax,1,'first');
        noise=std(data(1:bgregionmaxx));        
        Sdbl=SNRbulk*noise;
        SQthr=Qthr*noise;
        
        dead=false; blue=false; broad=false;
        
        % find first and last good spectra in terms of SNR
        maxdata=max(data);
%        lastgood=find(maxdata-(1.5*sqrt(maxdata)+noise)>SNRdead*noise,1,'last');   % factor 1.5 compensates for outliers. Not the best approx, e.g. small, single reds are missed!
        lastgood=[];
        for spec=size(data,2):-1:1
            [~,I]=max(data(:,spec));
            if (I>2)&&(I<length(data)-2)
                peakampl=sum(data(I-2:I+2,spec))/5;
                if peakampl/noise > SNRdead
                    lastgood=spec;
                    break;
                end
            end
        end
  
        if (isempty(lastgood)) 
            dead=true;
            deads(i)=complex;
        else
            firstgood=find(maxdata-(1.5*sqrt(maxdata)+noise)>Sdbl,1,'last')+1;
            if isempty(firstgood)
                firstgood=1;
            end
            dataminvib=mat(:,1:lastgood+1);             % data minus vib band
            [~,I]=max(mat(:,2:lastgood+1));
            I(I<3)=3; I(I>length(mat)-3)=length(mat)-3; % For indices of next command. This is anyway in the noise zone.
            peakampl=zeros(1,lastgood);
            for spec=firstgood:lastgood
                peakampl(spec)=sum(mat(I(spec)-2:I(spec)+2,spec+1))/5;
            end
            unQspecsnr=[false peakampl>SQthr];  % first index for ignoring wavelengths
            unQspecs=mat(:,unQspecsnr);         % All unQ specs
            %Alternative, less effective way
            %minvibn=find(mat(:,1)<=minvib,1,'last'); %otherwise vib wing can give misfits and confuse normal spectra with blue ones
            %unQspecsnr=[false,maxdata-(1.5*sqrt(maxdata)+noise)>SQthr];  
            %unQspecs=mat(1:minvibn,unQspecsnr); 
            if (isempty(unQspecs)) || (firstgood>lastgood)
                dead=true;
                deads(i)=complex;
            else
                % fitting boundaries for single peaks                
                lbs = [SQthr fwhmmin peakmin skewmin1];
                ubs = [Sdbl fwhmdouble peakmx skewmax1]; %width can be scaled with wavelength
                
                % fitting boudaries for vibration band
                if fitvib
                    lbvib = [noise 2/3*fwhmvib peakvibmin];
                    ubvib = [Sdbl fwhmvibmax peakvibmax];
                end 

                if manualdbl
                    toosmall=find(alldbl>lastgood);
                    alldbl=trimarray(alldbl,toosmall);
                end

          %%% 1. Identify blue or broad spectra by fitting average of all good (unQ) spectra %%%
                if size(unQspecs,2)>1
                    avgdata=sum(unQspecs,2)./size(unQspecs,2); %first index being the wavelength
                else
                    avgdata=unQspecs;
                end
                specfit=skewgaussfit4([wav,avgdata],[],fittype,[SQthr 0 lbs(3:4)],[Sdbl Inf ubs(3:4)]);
                if (~incfile)&&((specfit(2)>peakmx-peak1)||(specfit(2)>2*fwhmmax1+fwhmmax2))%&&(specfit(4)<skewmin1+.001))
                    broad=true;
                    broads(i)=complex;
                elseif (testbluing)&&(~incfile)  % vib band may cause ambiguity!!
                    if (specfit(3)<peakmin+.1)||...
                            ((specfit(3)<peak1-fwhm1/2.5)&&(specfit(2)>fwhmdouble-.1))||((specfit(4)<skewmin1+.01)&&(specfit(3)<peak1+fwhm1/2)) %normal+blue not included here
                       blue=true;
                       blues(i)=complex; 
                    end
                end

          %%% 2. Each spectrum is fitted with a single band. Blues and broads are again identified  %%%
          %%%    Vib bands are resolved except for clearly double and red-shifted spectra           %%%

                allint1=zeros(maxspecs,1);
                allint2=zeros(maxspecs,1);
                allSNR=zeros(maxspecs,1);
                spec=firstgood;
                singlespecfit=-1*ones(4,lastgood);  %fit parameters of single-band fits
                allvibfit=-1*ones(3,lastgood);      %fit parameters of vib band
                fwhmdbl=zeros(1,lastgood);          %scales for every spectrum according to wavelength and intensity

                while (spec<=lastgood)&&(~blue)&&(~broad)
                    if spec==7
                    end
                    data=mat(:,spec+1);
                    int=sum(data);
                    SNR=peakampl(spec)/noise;                    
                    single=true;
                    if manualdbl
                        j=1;
                        while (j<=length(alldbl))&&(alldbl(j)<=spec)&&(single)
                            if spec==alldbl(j)
                                single=false;
                            end
                            j=j+1;
                        end
                    end
                    if single
                        vibbandfit=[0 0 0];
                        if SNR<Qthr     %SNR too small to resolve spectrum
                            specfit=[0 0 0 0];
                        else
                            if SNR>SNRsmall
                                data=adjavg_both(data,wavbin);
                            else
                                data=adjavg_both(data,wavbinsmall);
                            end
                            if scalewidth_int
                                ubs(2)=max(fwhmdbl(spec),fwhmdouble+(scalef_int/SNR)^2); %scale fwhmdbl according to intensity (Q states often become broader)
                            end    
                            fwhmdbl(spec)=ubs(2);
                            if SNR>SNRsmall
                                specfit=skewgaussfit4([wav,data],[],fittype,lbs,[Sdbl Inf ubs(3:4)]); % fitting without constraints is faster, but less effective
                            else
                                specfit=skewgaussfit4([wav,data],[fwhm1 peak1 skewness],fittype,lbs,[Sdbl Inf ubs(3:4)]);
                            end
%                            if (specfit(2)<lbs(2)+.1)      %misfit
%                                specfit=skewgaussfit4([wav,data],[],fittype,[],[]);
%                            end                            
                            if (specfit(1)<SQthr)||(specfit(3)<lbs(3)+.1)||(specfit(3)>ubs(3))  %SNR too small to resolve spectrum
                               specfit=[0 0 0 0];
                            elseif ((specfit(3)>peak1+fwhm1/2)&&(specfit(2)<2*fwhmmax1+fwhmmax2)) %if clearly double, don't fit vib band to save time
                                if (specfit(4)<lbs(4)+.1) || (specfit(4)>ubs(4)-.1) || (specfit(3)>ubs(3)-.1) || (specfit(2)>fwhmmax1+fwhmmax2)                           
                                    border=find(wav>specfit(3)-fwhm2,1,'first');
                                    specfit=skewgaussfit4([wav(border:end),data(border:end)],[fwhm2 specfit(3) skewness],fittype,[],[]);
                                    if (specfit(4)<lbs(4)+.1) || (specfit(4)>ubs(4)-.1) || (specfit(3)<max(wav(border),peak1+fwhm1/2)) || (specfit(3)>ubs(3)-.1)
                                         specfit=gaussfit3([wav(border:end),data(border:end)],fwhm2,peak2);
                                         specfit(4)=0;
                                         if (~incfile)&&(SNR>SNRdbl*2)&&(specfit(2)>max(peakmx-peak1,2*fwhmmax1+fwhmmax2))&&(specfit(3)<ubs(3))&&(specfit(3)>lbs(3)) %SNR bygevoeg
                                             broad=true;
                                             broads(i)=complex;
                                             break;
%                                         elseif (specfit(2)<lbs(2)+.1)||(specfit(2)>ubs(2)-.1)||(specfit(3)<max(wav(border),lbs(3)+.1))||(specfit(3)>ubs(3)-.1) %problematic for finding doubles!!
                                         elseif (specfit(2)<lbs(2)+.1)||(specfit(3)<max(wav(border),lbs(3)+.1))
                                             specfit=[-1 -1 -1 -1];
                                         end
                                    end
                                end
                            elseif (~incfile)&&(SNR>SNRdbl*2)&&(specfit(2)>2*fwhmmax1+fwhmmax2)%(specfit(2)>max(peakmx-peak1,2*fwhmmax1+fwhmmax2))
                                 broad=true;
                                 broads(i)=complex;
                                 break;
                            elseif (testbluing)&&(~incfile)&&(SNR>1.5*SNRdbl)&&(specfit(1)>SQthr)&&... % vib band may cause ambiguity!!
                                 ((specfit(3)<peakmin+.1)||((specfit(3)<peak1-fwhm1/2.5)&&(specfit(2)>fwhmdouble-.1)&&(specfit(4)<skewmax1-.01)))
                                blue=true;
                                blues(i)=complex;
                                break
                            else
                                if (fitvib)&&(SNR>2*Qthr) %fit vib band
                                    singlemax=find(wav>=minvib,1,'first');
                                    mainspecfit=skewgaussfit4([wav(1:singlemax) data(1:singlemax)],[fwhm1 peak1 skewness],fittype,lbs,ubs);
                                    if (SNR>2*SNRdbl)&&(testbluing)&&(~incfile)...
                                       &&((mainspecfit(4)<lbs(4)/2)||((mainspecfit(4)<lbs(4)/3)&&(mainspecfit(2)>ubs(2)-.1))) % potential bluing
                                        diffspec=data(1:singlemax)-skewgaussian3(mainspecfit,wav(1:singlemax));
                                        [A,I]=max(diffspec);
                                        peaki=find(wav>=peak1-fwhm1/2,1,'first');
                                        if (A>lbs(1))&&(I<peaki)
                                            blue=true;
                                            blues(i)=complex;
                                            break
                                        end
                                    end
                                    if (finddbltemp)&&(SNR>2*SNRdbl)&&((mainspecfit(2)>ubs(2)-.1)||(mainspecfit(4)>ubs(4)-.01))
                                        alldbl=[alldbl spec];
                                    else
                                        % dblskewgaussfit delivers too broad fits for vib band!
                                        vibmin=find(wav>=(minvib+peakvibmin)/2,1,'first'); %???
                                        vibband=data-gaussian(wav,mainspecfit(1),mainspecfit(2),mainspecfit(3)); 
                                        vibbandfit=gaussfit4([wav(vibmin:end) vibband(vibmin:end)],fwhmvib,peakvib,lbvib,ubvib);
                                        if vibbandfit(1)<lbvib(1)+.1
                                            vibbandfit=[0 0 0];
                                        %elseif (vibbandfit(3)<lbvib(3)+.1) - potential red band but not always
                                        elseif (vibbandfit(2)>ubvib(2)-.1)||(vibbandfit(3)<lbvib(3)+.1)||(vibbandfit(3)>ubvib(3)-.1) %misfit
                                            vibmin=find(wav>=minvib,1,'first'); 
                                            vibbandfit=gaussfit4([wav(vibmin:end) vibband(vibmin:end)],fwhmvib,peakvib,lbvib,ubvib);
                                            if (vibbandfit(2)>ubvib(2)-.1)||(vibbandfit(3)<lbvib(3)+.1)||(vibbandfit(3)>ubvib(3)-.1) %misfit
                                                vibmin=find(wav>=peakvibmin,1,'first');
                                                vibbandfit=gaussfit4([wav(vibmin:end) vibband(vibmin:end)],fwhmvib,peakvib,lbvib,ubvib);
                                                if (vibbandfit(2)>ubvib(2)-.1)||(vibbandfit(3)<lbvib(3)+.1)||(vibbandfit(3)>ubvib(3)-.1) %misfit
                                                    vibbandfit=[0 0 0];
                                                end
                                            end
                                        end
                                        if vibbandfit(1)==0
                                            specfit=mainspecfit;    %to avoid unnecessary fitting routine
                                        else
                                            if (mainspecfit(1)/vibbandfit(1)>=mainvsvibfactor) || (allSNR(spec)<=SNRdbl)
                                                data=data-gaussian(wav,vibbandfit(1),vibbandfit(2),vibbandfit(3));
                                                %specfit=skewgaussfit4([wav,data],[],fittype,lbs,[Sdbl Inf peakmax Inf]); %borders necessary for testing blues. But I don't want to restrict the width now
                                                specfit=skewgaussfit4([wav,data],[fwhm1 peak1 skewness],fittype,lbs,ubs);
                                            else
                                                vibbandfit=[0 0 0];
                                                alldbl=[alldbl spec]; %saves some execution time
                                            end
                                        end
                                    end
                                end
                                dataminvib(:,spec+1)=data;                            
    %                            if (testbluing)&&(intminvib>intdbl)&&(specfit(2)>fwhmdbl)&&(specfit(3)<(peak1+peak2)/2)&&(specfit(1)>10) %normal+blue
    %                               if (specfit(3)-specfit(2)<peak1-2*fwhm1)&&(specfit(4)<skewmax1/2)
%                                 if (testbluing)&&(~incfile)&&(SNR>SNRdbl*2)&&(specfit(1)>SQthr)
%                                    if (specfit(3)<peakmin+.1)||((specfit(3)<(peak1+peakmin)/2)&&(specfit(2)>fwhmdouble-.1)&&(specfit(4)<skewmax1))...  %doesn't work always!
%                                            ||((specfit(4)<skewmin1+.01)&&(specfit(3)<peak1+fwhm1/2))...
%                                            ||((specfit(3)-specfit(2)/2<peak1-fwhm1)&&(specfit(4)<skewmax1/2)) %fwhmmin ipv fwhm1??%normal+blue. What about 680+small 675???
%                                        blue=true;
%                                        blues(i)=complex;   
%                                    end
%                                 end
                            end
                        end
                        if ~blue
                            allint1(spec)=int;
                            singlespecfit(:,spec)=specfit';
                            allvibfit(:,spec)=vibbandfit';
                        end
                    end
                    allSNR(spec)=SNR;
                    spec=spec+1;
                end
                if sum(singlespecfit(1,:))==0
                    dead=true;
                    deads(i)=complex;
                end

                if (~blue)&&(~broad)&&(~dead)

            %%% 3. Good fits to single-band spectra and vib band are identified. Means and fitting boundaries are calculated. %%%

                    goodvibs=[]; badsinglesbool=[]; goodsinglesbool=[];
                    meansinglespecs=[]; stdsinglespecs=[]; meanvibs=[];
                    if constrainblue
    %                    goodsingles=singlespecfit(:,singlespecfit(1,:)>0 & singlespecfit(2,:)<fwhmmax1 & singlespecfit(3,:)<peak1+fwhm1/2);
                        if scalewidth_int
                            goodsinglesbool=(singlespecfit(1,:)>0 & singlespecfit(2,:)<fwhmdbl & singlespecfit(3,:)<peak1+fwhm1/2);
                        else
                            goodsinglesbool=(singlespecfit(1,:)>0 & singlespecfit(2,:)<ubs(2) & singlespecfit(3,:)<peak1+fwhm1/2);
                        end
                        goodsinglesbool(alldbl)=0;
                        if fitvib
                            %badvibsbool=allvibfit(2,:)>fwhmvibmax-.1 | (allvibfit(3,:)<peakvibmin+.1 & allvibfit(3,:)>0) | allvibfit(3,:)>peakvibmax-.1 | singlespecfit(3,:)>peak1+fwhm1/2;
                            badvibsbool=allvibfit(1,:)<=0 & singlespecfit(3,:)<peak1+fwhm1/2; 
                            goodvibsbool=~badvibsbool;
    %                        goodvibsbool(alldbl)=0; %contained in next command
                            goodvibsbool=goodsinglesbool.*goodvibsbool;
                            goodsinglesbool=goodvibsbool; %misfit in one <=> misfit in other
                            goodvibs=allvibfit(:,goodvibsbool.*(1:lastgood)>0); %all spec parameters before firstgood are set to 0
                            goodvibs=goodvibs(:,goodvibs(1,:)>0);
                            if size(goodvibs,2)>1
                                meanvibs=mean(goodvibs,2);
                                stdvibs=std(goodvibs,0,2);
                            elseif size(goodvibs,2)==1
                                meanvibs=goodvibs;
                                stdvibs=[sqrt(goodvibs(1)),abs(meanvibs(2)-fwhmvib),abs(goodvibs(3)-peakvib)];
                            end
                        end
                        goodsingles=singlespecfit(:,goodsinglesbool.*(1:lastgood)>0); 
                        badsinglesbool=~goodsinglesbool;
                        badsinglesbool(singlespecfit(1,:)==0 | singlespecfit(3,:)>peak1+fwhm1/2)=0;    %single/double red peaks
                        badsinglesbool(alldbl)=0;
                        badsinglesnr=badsinglesbool(firstgood:lastgood).*(firstgood:lastgood); 
                        badsinglesnr=badsinglesnr(badsinglesnr>0);
                        dataminvib(:,badsinglesbool(firstgood:lastgood).*(firstgood:lastgood)>0)=mat(:,badsinglesbool(firstgood:lastgood).*(firstgood:lastgood)>0);
                        if size(goodsingles,2)>=1
                            if size(goodsingles,2)>1
                                meansinglespecs=mean(goodsingles,2);
                                stdsinglespecs=std(goodsingles,0,2);
                            elseif size(goodsingles,2)==1
                                meansinglespecs=goodsingles;
                                meansinglespecs(2)=(meansinglespecs(2)+fwhm1)/2; %due to potential broadening of Q peaks
                                stdsinglespecs=sqrt(abs(meansinglespecs));
                                stdsinglespecs(2)=abs(meansinglespecs(2)-fwhm1);
                                stdsinglespecs(3)=abs(meansinglespecs(3)-peak1);
                            end
                            if stdsinglespecs(3)<1
                                stdsinglespecs(3)=1;
                            end
                            newlbs=[SQthr meansinglespecs(2)-2*stdsinglespecs(2) meansinglespecs(3)-3*stdsinglespecs(3) meansinglespecs(4)-2*stdsinglespecs(4)];
                            newubs=[Sdbl meansinglespecs(2)+2*stdsinglespecs(2) meansinglespecs(3)+3*stdsinglespecs(3) meansinglespecs(4)+2*stdsinglespecs(4)];
                        else
                            newlbs=lbs;
                            newubs=ubs;
                        end
                    end

             %%% 4. Guess vib band for misfitted ones, based on good fits %%%
                    if (fitvib)&&(~isempty(badsinglesnr))
                        [singlespecfit,dataminvib]=guessvibband(mat,allvibfit,singlespecfit,dataminvib,badsinglesnr,fwhmdbl,meanvibs,meansinglespecs,newlbs,newubs,noise,mainvsvibfactor,fwhmvib,peakvib);
                    end

             %%% 5. Identify doubles in three different ways, apart from earlier and later methods %%%
                    if (finddbltemp)&&(~manualdbl)
                        if isempty(alldbl)
                            [alldbl alldbl3] = finddoubles(alldbl,dataminvib,singlespecfit,allSNR(1:lastgood),noise,fwhmjump,firstgood,lastgood,SNRdbl,Qthr,fwhmdbl,2*fwhmmax1+fwhmmax2,peak1,fwhm1,skewred);
                        else  % if doubles have already been found, make width, skewness, and intensity thresholds less strict
                            [alldbl alldbl3] = finddoubles(alldbl,dataminvib,singlespecfit,allSNR(1:lastgood),noise,fwhmjump,firstgood,lastgood,SNRdbl/2,Qthr,fwhmdouble*ones(1,lastgood),2*fwhmmax1+fwhmmax2,peak1,fwhm1,skewred/2);
                        end
                       %keep alldbl3 separate for fitting purposes. Intensity of alldbl3 and alldbl5 could be smaller than other doubles
                       [alldbl ~] = addQdbl(alldbl,allSNR(1:lastgood)/SNRdbl);
                    end

             %%% 6. Refit misfitted singles and vib bands, except alldbl, based on fit parameters of good fits %%%
             %%%    If still misfitted, it could be another double
              
                    if (constrainblue)&&(~isempty(goodsingles))
                        if (fitvib)&&(~isempty(goodvibs))
                             if stdvibs(3)<2
                                stdvibs(3)=2;
                            end
                            newlbvib=[noise meanvibs(2)-3*stdvibs(2) max(meanvibs(3)-2*stdvibs(3),peakvibmin)];
                            newubvib=[Sdbl meanvibs(2)+3*stdvibs(2) min(meanvibs(3)+2*stdvibs(3),peakvibmax)];
                        end
                        for spec=badsinglesnr
                            if spec==7
                            end
                            data=mat(:,spec+1);
                            newubs(2)=fwhmdbl(spec);
                            if allSNR(spec)>SNRsmall
                                data=adjavg_both(data,wavbin);
                            else
                                data=adjavg_both(data,wavbinsmall);
                            end
                            if (fitvib)&&(~isempty(goodvibs))&&(allSNR(spec)>2*Qthr)
    %                            dblspecfit=dblskewgaussfit2([wav data],[meansinglespecs(2),meansinglespecs(3),meansinglespecs(4),meanvibs(2),meanvibs(3),meanvibs(4)],[newlbs newlbvib],[newubs newubvib]);
    %                            if (dblspecfit(2)<newlbs(2)+.1)||(dblspecfit(5)<noise)||(dblspecfit(6)>newlbvib(2)-.1)||(dblspecfit(7)>newlbvib(3)-.1) %misfit, then try to fit to bands separately   
                                mainspecfit=skewgaussfit4([wav(1:singlemax) data(1:singlemax)],[meansinglespecs(2),meansinglespecs(3),meansinglespecs(4)],fittype,newlbs,newubs);
                                vibband=data-gaussian(wav,mainspecfit(1),mainspecfit(2),mainspecfit(3)); 
                                vibbandfit=gaussfit4([wav(vibmin:end) vibband(vibmin:end)],meanvibs(2),meanvibs(3),newlbvib,newubvib);
                                %For the sake of saving an extra for loop, I'm  taking this shortcut to approximate misfits
                                %Normally they're just too broad or skew
                                if vibbandfit(1)<newlbvib(1)+.1
                                    allvibfit(:,spec)=[0 0 0];
    %                                 singlepeak2=data-skewgaussian3([mainspecfit(1),mainspecfit(2),mainspecfit(3),mainspecfit(4)],wav);
    %                                 if sum(singlepeak2(singlemax:end))>intfarred
    %                                     alldbl=[alldbl specnr];
    %                                 end
    %                             elseif vibbandfit(2)>newubvib(2)-.1
    %                                 index=find(wav>=meansinglespecs(3),1,'first'); %- change: just leave as it is. Approx already calculated above!!!
    %                                 amplblue=sum(data(index-1:index+1))/3;
    %                                 vibpeakampl=amplblue*meanvibs(1)/meansinglespecs(1);
    %                                 vibbandfit=[vibpeakampl,meanvibs(2),meanvibs(3)]; %approximation
                                elseif vibbandfit(2)<newubvib(2)-.1
                                    goodvibs=[goodvibs vibbandfit'];
                                    meanvibs=mean(goodvibs,2);
                                    stdvibs=std(goodvibs,0,2);
                                    allvibfit(:,spec)=vibbandfit';
                                end
                                data=data-gaussian(wav,vibbandfit(1),vibbandfit(2),vibbandfit(3)); 
                                specfit=skewgaussfit4([wav,data],[],fittype,newlbs,newubs);
                                dataminvib(:,spec+1)=data;
                                singlepeak=skewgaussian3([specfit(1),specfit(2),specfit(3),specfit(4)],wav); 
                                if ~isempty(alldbl)   %if there is a big misfit, then prob double band
                                    bordermin=find(wav>=specfit(3)-specfit(2)*1.5,1,'first');
                                    bordermax=find(wav>=specfit(3)+specfit(2)*2,1,'first');
                                    rms=sqrt(sum((singlepeak(bordermin:bordermax)-data(bordermin:bordermax)).^2)/length(singlepeak));
                                    if rms>rmsthr  %not 100% accurate, but identifies some more doubles. Threshold should change for noisy data!! (Check this again)
                                        alldbl=[alldbl spec];  % Shift this to "findsmalldbl.m"!!!
                                    end
                                end
                            else
                                specfit=skewgaussfit4([wav,data],meansinglespecs(2:4)',fittype,[],[]);
                                if (sum(specfit(1:3)<lbs(1:3))>0)||(sum(specfit(1:3)>ubs(1:3))>0)
                                    specfit=skewgaussfit4([wav,data],meansinglespecs(2:4)',fittype,newlbs,newubs);
                                    if (sum(specfit(1:3)<lbs(1:3)+.1)>0)||(sum(specfit(1:3)>ubs(1:3)-.1)>0)
                                        specfit=[-3 -3 -3 -3];
                                    end
                                end
                            end
                            if specfit(3)>peak1+fwhm1/2
                                allvibfit(:,spec)=[0 0 0];
                            end
                            singlespecfit(:,spec)=specfit';
                        end
                        if size(goodsingles,2)<=3 % otherwise the mean will be more or less good
                            meansinglespecs=mean(goodsingles,2);
                            stdsinglespecs=std(goodsingles,0,2);
                        end
                    else
                        newlbs=lbs;
                        newubs=ubs;
                    end

                    allspecfit=-2*ones(8,lastgood);
                    allspecfit(1:4,:)=singlespecfit;

            %%% 7. Look for doubles in the region of the vib band and further into the red, after subtracting the better estimated vibrational bands
            %%%   (If very small (quenched) far-red spectra are missed, these Q spectra should be added to "allsingle"!!) 
                
                    alldbl5=[];
                    if (fitvib)&&(finddbltemp)&&(findfarreds)&&(~manualdbl)
                        if ~isempty(meansinglespecs)
                            alldblbool=zeros(1,lastgood);
                            alldblbool(alldbl)=1;       % check only spectra not already identified as doubles
                            allsingle=~alldblbool.*(1:lastgood); 
                            allsingle=allsingle(allsingle>0);
                            [alldbl4,alldbl5] = findsmalldbl(mat,alldbl,allvibfit,singlespecfit,allsingle,meanvibs,meansinglespecs,peak1,fwhm1,skewred,intfarred,farredminwav,allSNR/SNRdbl);
                            
                            %NEW!!
                            for spec=alldbl
                                if (singlespecfit(2,spec)>fwhmmax1+fwhmmax2)&&(singlespecfit(3,spec)>peak1+fwhm1*3/4)
                                    alldbl5=[alldbl5 spec];
                                end
                            end
                            
                            if ~isempty(alldbl4)
                                alldbl=sort([alldbl alldbl4]);
                            end
                            if ~isempty(alldbl5)
                                alldbl=sort([alldbl alldbl5]);
                                alldbl5=sort([alldbl3 alldbl5]);
                            else
                                alldbl5=alldbl3;
                            end 
                        end
                    end 

           %%% 8. Fit double-band spectra (I'll make separate function later) %%%
           
                     if (finddbltemp)&&(~isempty(alldbl))
                        if ~manualdbl
                            [alldbl Qdbl] = addQdbl(alldbl,allSNR(1:lastgood)/SNRdbl);
                            if (~isempty(alldbl5))&&(~isempty(Qdbl))
                                alldbl5=sort([alldbl5 Qdbl]);
                            end
                            allspecfit(:,Qdbl)=0;
                            allspecfit(:,alldbl)=-3;

                            [trans,transi] = transitdoubles(alldbl,singlespecfit,allspecfit,peak1,fwhm1); %use fewer variables??
                            alldbl=trimarray(alldbl,transi);

                            if ~isempty(alldbl5)
                                alldbl5=remvalarray(alldbl5,trans);  %remove trans from alldbl5
                                alldbl=remvalarray(alldbl,alldbl5);   %remove alldbl5 from alldbl
                            end
                        end

                        if (constrainblue)&&(~isempty(goodsingles))
                            if size(goodsingles,2)>1 %at least 3 contribute to std
                                olbds = [SQthr meansinglespecs(2)-2*stdsinglespecs(2),meansinglespecs(3)-2*stdsinglespecs(3),meansinglespecs(4)-stdsinglespecs(4) SQthr fwhmmin max(peak2min,peak1+fwhmmin*2/3) skewmin2]; % 3/4!!!
                                oubds = [Sdbl meansinglespecs(2)+2*stdsinglespecs(2),meansinglespecs(3)+2*stdsinglespecs(3),meansinglespecs(4)+stdsinglespecs(4) Sdbl fwhmmax2 peakmx skewmax2]; % double skew
                            else %only 1 good fit
                                olbds = [SQthr fwhmmin goodsingles(3)-1.5 skewmin1 SQthr fwhmmin max(peak2min,peak1+fwhmmin/2) skewmin2/2]; % double skew
                                oubds = [Sdbl fwhmmax1 goodsingles(3)+1.5 skewmax1/2 Sdbl fwhmmax2 peakmx skewmax2]; % double skew    
                            end
                        else 
                            olbds = [SQthr fwhmmin peakmin skewmin1 SQthr fwhmmin max(peak2min,peak1+fwhmmin*2/3) skewmin2]; % double skew
                            oubds = [Sdbl fwhmmax1 peak1+fwhmmin/2 skewmax1 Sdbl fwhmmax2 peakmx skewmax2]; % double skew
                        end
                        if isempty(meansinglespecs)
                            singlefitguess=[noise*10 fwhm1 peak1 skewness]; % guess starting values
                        else
                            singlefitguess=meansinglespecs;
                        end
                        alloverlaps=zeros(1,lastgood);

                        for spec=alldbl5 % special case for alldbl3 and alldbl5. (Later a distinction between the two can be made. Then I need to perform transitdoubles on both and separate loops spec=alldbl3 and spec=alldbl5)
                                % This means that in (hopefully) all cases the two peaks are well separated
                                % (Gaussian can be fitted to peaks with sufficiently small intensities)
                            ubds=oubds; lbds=olbds; lbds(8)=lbs(4);
                            if allSNR(spec)<Qthr
                                allspecfit(:,spec)=zeros(1,8);
                            else
                                data=mat(:,spec+1);
                                if fwhmdbl(spec)>0
                                    newubs(2)=fwhmdbl(spec);
                                elseif (scalewidth_int)&&(specfit(3)<peak1+fwhm1/2)
                                    newubs(2)=fwhmdouble+(scalef_int/allSNR(spec))^2; %scale fwhmdbl according to intensity 
                                else
                                    newubs(2)=fwhmdouble;
                                end
                                ubds(6)=newubs(2);
                                data=adjavg_both(data,wavbinsmall);
                                %subtract vib band - important for accurate fitting of reds
                                if (fitvib)&&(singlespecfit(3,spec)<peak1+fwhm1/2)
                                    if (~isempty(meansinglespecs))&&(~isempty(goodvibs))
                                        index=find(wav>=meansinglespecs(3),1,'first');
                                        amplblue=sum(data(index-1:index+1))/3;
                                        vibpeakampl=amplblue*meanvibs(1)/meansinglespecs(1);
                                        dblvib=[vibpeakampl,meanvibs(2),meanvibs(3)];
                                    else %guess
                                        index=find(wav>=peak1,1,'first');
                                        amplblue=sum(data(index-1:index+1))/3;
                                        vibpeakampl=amplblue/(mainvsvibfactor*1.5);
                                        dblvib=[vibpeakampl,fwhmvib,peakvib];
                                    end
                                    allvibfit(:,spec)=dblvib;
                                    vibband=gaussian(wav,dblvib(1),dblvib(2),dblvib(3));
                                    data=data-vibband;
                                    dataminvib(:,spec+1)=data;
                                end
        %                            dblspecfit=dblskewgaussfit2([wav data],[meansinglespecs(2),meansinglespecs(3),meansinglespecs(4),fwhm2,peak2,skewness],lbds,ubds);      %doublefit always worse than 2 separate single fits!
                                midwav=find(mat(:,1)>peak1+2*fwhm1,1,'first');
                                overlap=false;
                                specfit1=skewgaussfit4([wav(1:midwav),data(1:midwav)],[singlefitguess(2),singlefitguess(3),singlefitguess(4)],fittype,newlbs,newubs); % should restrict borders??
                                if specfit1(1)<newlbs(1)+.1  %only red band
                                    specfit1=skewgaussfit4([wav(midwav:end),data(midwav:end)],[],fittype,lbds(5:8),ubds(5:8));
                                    if (specfit1(1)<lbds(5)+.1) || (specfit1(2)>ubds(6)-.1) || (specfit1(2)<lbds(6)+.1) || (specfit1(4)<lbds(8)+.001) %misfit, prob too noisy to fit skewed gaussian
                                        specfit1=[gaussfit3([wav(midwav:end),data(midwav:end)]) 0];
                                        if (specfit1(1)<lbds(1))||(specfit1(2)<lbds(2))||(specfit1(2)>ubds(2))||(specfit1(3)<lbds(3))||(specfit1(3)>ubds(7))
                                            specfit1=[0 0 0 0];
                                        end
                                    end
                                    specfit2=[0 0 0 0];
                                    alldbl5=remvalarray(alldbl5,spec); 
                                    allvibfit(:,spec)=[0 0 0];
                                else
                                    if (specfit1(2)>newubs(2)-.1)||(specfit1(3)>newubs(3)-.1) %two peaks prob not well separated
                                        singlepeak2=data-skewgaussian3([specfit1(1),specfit1(2),specfit1(3),specfit1(4)],wav);
                                        [~,I]=max(singlepeak2);
                                        if (I>2)&&(I<length(singlepeak2)-2)
                                            peakamp=sum(singlepeak2(I-2:I+2))/5;
                                        else
                                            peakamp=SQthr;
                                        end
                                        if peakamp<SQthr+.1
                                            overlap=true;
                                        else
                                            specfit2=skewgaussfit4([wav,singlepeak2],[],fittype,[],[]);
                                            if (specfit2(3)-specfit1(3)<specfit1(2)+specfit2(2))||(specfit2(2)<fwhmmin+.1)||(specfit2(1)<SQthr+.1) % then it's very likely that this is a single-band spectrum
                                                overlap=true;
                                            elseif specfit(4)<-.1
                                                specfit2=gaussfit3([wav,singlepeak2]);
                                                specfit2(4)=0;
                                            end
                                        end
                                        if overlap %not single but too strong overlap
                                           alldbl=[alldbl spec];
                                           alldbl5=remvalarray(alldbl5,spec);
                                           alloverlaps(spec)=1;
                                        else
                                           alloverlaps(spec)=0;
                                        end
                                    elseif (specfit1(2)<newlbs(2)+.1) %misfit
                                        specfit1=skewgaussfit4([wav(1:midwav),data(1:midwav)],[singlefitguess(2),singlefitguess(3),singlefitguess(4)],fittype,newlbs,newubs);
                                        if (specfit1(2)<newlbs(2)+.1) || (specfit1(2)>newubs(2)-.1) %still misfit
                                            dblspecfit=[-1 -1 -1 -1 -1 -1 -1 -1];
                                            overlap=true;
                                        end
                                    end
                                    if ~overlap
            %                             int1=sum(data(1:midwav));
            %                             int2=allint1(spec)-int1;
                                        singlepeak2=data-skewgaussian3(specfit1,wav);
                                        specfit2=skewgaussfit4([wav(midwav:end),singlepeak2(midwav:end)],[],fittype,lbds(5:8),ubds(5:8));
                                        if (specfit2(2)>ubds(6)-.1)&&(scalewidth_wav)
                                            ubds(6)=max(ubds(6),oubds(6)+(specfit2(3)-specfit1(3))/scalef_wav);
                                            specfit2=skewgaussfit4([wav(midwav:end),singlepeak2(midwav:end)],[],fittype,lbds(5:8),ubds(5:8));
                                        end
                                        if specfit2(4)<lbds(8)+.001 %too negatively skewed, likely fitting part of blue band
                                            midwav=midwav+20;   %can multiply with spectral resolution
                                            specfit2=skewgaussfit4([wav(midwav:end),singlepeak2(midwav:end)],[],fittype,lbds(5:8),ubds(5:8));
                                        end
                                        if (specfit2(1)<SQthr+.1) || (specfit2(2)>ubds(6)-.1) || (specfit2(2)<lbds(6)+.1) || (specfit2(3)<lbds(7)+.1) || (specfit2(3)>ubds(7)-.1) || (specfit2(4)<-.1) %misfit, prob too noisy to fit skewed gaussian, and neg skewness doesn't make sense
                                            specfit2=gaussfit3([wav(midwav:end),singlepeak2(midwav:end)]);
                                            specfit2(4)=0;
                                            if (specfit2(1)<SQthr/2) || (specfit2(2)>ubds(6)-.1) || (specfit2(2)<lbds(6)+.1) || (specfit2(3)<lbds(7)+.1) || (specfit2(3)>ubds(7)-.1) %too noisy for double band
                                                allspecfit(:,spec)=[-1 -1 -1 -1 -1 -1 -1 -1];
                                                overlap=true;
                                            end
                                        end
                                    end
                                end
                                if ~overlap
                                    allspecfit(:,spec)=[specfit1 specfit2];
                                end
                            end
                        end

                        countdbl=1;
                        for spec=alldbl
                            data=mat(:,spec+1);
                            ubds=oubds; lbds=olbds;
                            if fwhmdbl(spec)>0
                                newubs(2)=fwhmdbl(spec);
                            elseif (scalewidth_int)&&(specfit(3)<peak1+fwhm1/2)
                                newubs(2)=fwhmdouble+(scalef_int/allSNR(spec))^2; %scale fwhmdbl according to intensity 
                            else
                                newubs(2)=fwhmdouble;
                            end
                            
                            if (allSNR(spec)<SNRdbl)&&(alloverlaps(spec)==0)
                                dblspecfit=zeros(1,8);
                            else
                                if allSNR(spec)<SNRsmall
                                    data=adjavg_both(data,wavbinsmall);
                                else
                                    data=adjavg_both(data,wavbin);
                                    %subtract vib band - important for accurate fitting of reds
                                    if fitvib
                                        if (~isempty(meansinglespecs))&&(~isempty(goodvibs))
                                            index=find(wav>=meansinglespecs(3),1,'first');
                                            amplblue=sum(data(index-1:index+1))/3;
                                            vibpeakampl=amplblue*meanvibs(1)/meansinglespecs(1);
                                            dblvib=[vibpeakampl,meanvibs(2),meanvibs(3)];
                                        else %guess
                                            index=find(wav>=peak1,1,'first');
                                            amplblue=sum(data(index-1:index+1))/3;
                                            vibpeakampl=amplblue/(mainvsvibfactor*1.5);
                                            dblvib=[vibpeakampl,fwhmvib,peakvib];
                                        end
                                        allvibfit(:,spec)=dblvib;
                                        vibband=gaussian(wav,dblvib(1),dblvib(2),dblvib(3));
                                        data=data-vibband;
                                        dataminvib(:,spec+1)=data;
                                    end
                                end

                                dblspecfit=dblskewgaussfit2([wav data],[singlefitguess(2),singlefitguess(3),singlefitguess(4),fwhm2,peak2,skewness],lbds,ubds);
    %                             if dblspecfit(7)-dblspecfit(3)>dblspecfit(2)+dblspecfit(6) % well separated peaks: use fitting procedure of alldbl5
    %                                 alldbl5=[alldbl5 spec];
    %                             end
                                if (dblspecfit(6)>ubds(6)-.1)&&(dblspecfit(7)>dblspecfit(3)+dblspecfit(2)+dblspecfit(6))&&(scalewidth_wav) %only for well separated peaks
                                    ubds(6)=oubds(6)+(dblspecfit(7)-dblspecfit(3))/scalef_wav; %scale fwhmdbl according to wavelength
    %                                ubs2(2)=ubds(6);
                                    dblspecfit=dblskewgaussfit2([wav data],[singlefitguess(2),singlefitguess(3),singlefitguess(4),fwhm2,peak2,skewness],lbds,ubds); %can also use starting values of dblspecsfit!!
                                end

    %                             if (testbluing)&&(~incfile)&&(allSNR(spec)>SNRdbl*2)&&(dblspecfit(4)<skewmin1*3/4)||((dblspecfit(4)<skewmin1/2)&&(dblspecfit(2)>ubds(2)-.1))
    %                                blue=true;
    %                                blues(i)=complex;
    %                                break; 
    %                             end
    %                             while (dblspecfit(8)>ubds(8)-.01)&&(dblspecfit(7)>dblspecfit(3)+dblspecfit(2))&&(ubds(6)<fwhm2*2)
    %                                 ubds(6)=ubds(6)+5;
    %                                 dblspecfit=dblskewgaussfit2([wav data],[singlefitguess(2),singlefitguess(3),singlefitguess(4),fwhm2,peak2,skewness],lbds,ubds);
    %                             end
                                if dblspecfit(1)<SQthr+.1
        %                           specfit=skewgaussfit4([wav,data],meansinglespecs(2:4),fittype,lbs,ubs);   - necessary to fit again??
        %                           if (specfit(2)>fwhmdbl(spec)-.1)   % single-band spectrum; scale fwhmdbl acc to wavelength and int
        %                                 int2=0; int1=int;
        %                                 specfit1=skewgaussfit4([wav,data],[fwhm1,peak1,skewness],fittype,lbs,ubs);
        %                                 if (specfit1(2)>fwhmmax1-.1)
        %                                    ub=[Inf fwhmmax1+5 peakmax skewmax2];
        %                                    specfit1=skewgaussfit4([wav,data],[fwhm1,peak1,skewness],fittype,lbs,ub);
        %                                 end

                                     alldbl=trimarray(alldbl,countdbl);
        %                             dblspecfit=[specfit1 specfit2];
        %                            end
                                elseif (sum(dblspecfit([2,3,6,7])<lbds([2,3,6,7])+.1)>0)||...
                                        (sum(dblspecfit([2,3,6,7])>ubds([2,3,6,7])-.1)>0)||...
                                        (dblspecfit(4)<lbds(4)+.01)||(dblspecfit(8)>ubds(8)-.01)   %misfit
                                    midwavx=min(dblspecfit(3)+dblspecfit(7),peak1+peak2)/2;
                                    midwav=find(wav>midwavx,1,'first');
                                    specfit1=skewgaussfit4([wav(1:midwav),data(1:midwav)],[dblspecfit(2),dblspecfit(3),skewness],fittype,[],[]);  %fit not good if restricted
                                    if (specfit1(2)>2*fwhmmax1+fwhmmax2)||(specfit1(3)<peakmin)||(specfit1(3)>dblspecfit(7)) %misfit
                                        specfit1=skewgaussfit4([wav(1:midwav),data(1:midwav)],[dblspecfit(2),dblspecfit(3),skewness],fittype,lbs,ubs); 
                                    end
                                    %fit 2 bands separately if overlap not too strong
                                    if (testbluing)&&(~incfile)&&(allSNR(spec)>SNRdbl*2)&&(specfit1(4)<skewmin1/2)&&(specfit1(3)<peak1)&&(specfit1(2)>fwhmmin)
                                       blue=true;
                                       blues(i)=complex;
                                       break; 
                                    elseif (specfit1(2)<newubs(2)-.01)%||(dblspecfit(7)>peak1+fwhm1/2) %fit is good (or width at least is not to large)
                                        border=midwav;
                                        borderx=midwavx; 
                                        if (specfit1(2)<newlbs(2)+.01)&&(dblspecfit(6)<ubds(6)-.1) %both peaks likely resolved, but blue one too narrow
                                            inc=3;
                                        else
                                            inc=10; 
                                        end
                                        while (borderx<dblspecfit(7))&&((specfit1(2)<fwhmmin+.1)||(specfit1(3)<peakmin)||(specfit1(3)>peakmx)) %misfit
                                            border=border+inc;
                                            borderx=wav(border);
                                            ub = [Sdbl newubs(2) borderx skewmax2];
                                            specfit1=skewgaussfit4([wav(1:border),data(1:border)],[dblspecfit(2),dblspecfit(3),skewness],fittype,newlbs,ub); 
                                            if specfit1(2)>ub(2)
                                                specfit2=[-1 -1 -1 -1];
                                                break
                                            end
                                        end
                                        singlepeak2=data-skewgaussian3([specfit1(1),specfit1(2),specfit1(3),specfit1(4)],wav);                                        
%                                        specfit2=skewgaussfit4([wav,singlepeak2],[dblspecfit(6),dblspecfit(7),dblspecfit(8)],fittype,lbds(5:8),ubds(5:8));  %contraints cause more misfits than without
                                        specfit2=skewgaussfit4([wav,singlepeak2],[dblspecfit(6),dblspecfit(7),dblspecfit(8)],fittype,[],[]);
                                        if (specfit2(2)<lbds(6)+.01)||(specfit2(3)<specfit1(3)+specfit1(2)/2)||(specfit2(1)<lbds(5)+.1)  % too much taken by spec1. This part called by external function? Same as above. Put in loop?
                                            border=midwav-5;            %2nd condition=>too strong overlap
                                            borderx=wav(border);
                                            inc=2;
                                            ub = [Sdbl newubs(2) borderx skewmax2];
                                            specfit1=skewgaussfit4([wav(1:border),data(1:border)],[specfit1(2)-2,specfit1(3)-1,skewness],fittype,newlbs,ub); 
                                            while (border<midwav+3)&&((specfit1(2)<newlbs(2)+.1)||(specfit1(2)>ub(2)-.1)||(specfit1(3)<newlbs(3)+.1)||(specfit1(3)>ub(3)-.1)||(specfit1(4)>ub(4)-.01)) %misfit
                                                border=border+inc;
                                                borderx=wav(border);
                                                ub = [Sdbl newubs(2) borderx skewmax2];
                                                specfit1=skewgaussfit4([wav(1:border),data(1:border)],[dblspecfit(2),dblspecfit(3),skewness],fittype,newlbs,ub); 
                                                if specfit1(2)>ub(2)
                                                    specfit2=[-1 -1 -1 -1];
                                                    break
                                                end
                                            end
                                            singlepeak2=data-skewgaussian3([specfit1(1),specfit1(2),specfit1(3),specfit1(4)],wav);
                                            %check first SNR before refit -faster, but not always effective
%                                             [~,I]=max(singlepeak2);
%                                             if (I>2)&&(I<length(singlepeak2)-2)
%                                                 peakamp=sum(singlepeak2(I-2:I+2))/5;
%                                             else
%                                                 peakamp=SQthr;
%                                             end
%                                             if peakamp<SQthr+.1
%                                                 specfit1=singlespecfit(:,spec)';
%                                                 specfit2=[-2 -2 -2 -2];
%                                                 alldbl=trimarray(alldbl,countdbl);
%                                                 countdbl=countdbl-1;
%                                                 disp('increase fwhmdoubles or scalewidth_int!!');
%                                             else
                                                specfit2=skewgaussfit4([wav,singlepeak2],[dblspecfit(6),dblspecfit(7),dblspecfit(8)],fittype,[],[]);
                                                if (specfit2(2)<fwhmmin+.01)||(specfit2(1)<SQthr+.1)||(specfit2(3)<specfit1(3)+specfit1(2)/2) % misfit
                                                    specfit2=[-1 -1 -1 -1];
                                                end
%                                             end
                                        else  %prob well separated peaks
                                            border=midwav;
                                            borderx=wav(border);
                                            if scalewidth_wav
                                                ubds(6)=max(ubds(6),oubds(6)+(specfit2(3)-specfit1(3))/scalef_wav);
                                            end
                                            while (border>midwav/2)&&((specfit2(1)<SQthr+.1)||(specfit2(2)<fwhmmin+.01)||(specfit2(2)>ubds(6)-.1)||(specfit2(3)>peakmx-1)||(specfit2(3)<borderx+1)||(specfit2(4)<lbs(4)+.01))
                                                lb=[SQthr fwhmmin borderx lbs(4)];
                                                specfit2=skewgaussfit4([wav(border:end),singlepeak2(border:length(wav))],[dblspecfit(2),dblspecfit(3),dblspecfit(4)],fittype,lb,ubds(5:8));
                                                border=border-20;
                                                borderx=wav(border);
                                            end

        %                                    if (border<midwav/2)||(abs(specfit2(3)-dblspecfit(7))>10)     % misfit
                                            if border<midwav/2     % misfit
                                                specfit2=[-1 -1 -1 -1];
%                                                ub=[Sdbl fwhmmax2 midwavx skewmax2];
%                                                specfit1=skewgaussfit4([wav,data],[fwhm1,peak1,skewness],fittype,lbs,ub);
                                            elseif specfit2(4)<0
                                                specfit2=gaussfit3([wav(border:end),singlepeak2(border:length(wav))]);
                                                specfit2(4)=0;
                                            end
                                        end
                                    else %overlap of 2 bands too strong: use good single- or double-band fitting data in misfit.m
                                        specfit2=[-1 -1 -1 -1];
                                    end
                                    dblspecfit=[specfit1 specfit2];
                                % If an additional small red peak is fitted in the tail, fit with single peak only
                                elseif (singlespecfit(2,spec)<fwhmdbl(spec)+5)...
                                      &&((dblspecfit(1)/dblspecfit(5)>3)||((dblspecfit(1)/dblspecfit(5)>2.5)&&(dblspecfit(7)-dblspecfit(3)<20))) % check this again!!
                                    dblspecfit=[singlespecfit(:,spec)' -2 -2 -2 -2];
                                    alldbl=trimarray(alldbl,countdbl);
                                    countdbl=countdbl-1;
                                end
                            end
                            % check misfits
                            if (dblspecfit(5)>0)&&(dblspecfit(1)/dblspecfit(5)>6)&&(dblspecfit(7)<dblspecfit(3)+dblspecfit(2)+dblspecfit(6))       % small red band fitted in vib tail. check factor!
                                allspecfit(:,spec)=[singlespecfit(:,spec)' 0 0 0 0];
                            elseif (dblspecfit(4)<lbds(4)/3)&&(dblspecfit(8)>ubds(8)/2)&&(dblspecfit(2)+dblspecfit(6)>fwhmdouble*2)  
                                allspecfit(:,spec)=zeros(1,8);                          % broad spectrum misfitted
                            else
                                allspecfit(:,spec)=dblspecfit;
                            end
                            countdbl=countdbl+1;
                        end
                        
               %%% 9. Refit misfitted doubles %%%
                        
                        if ~blue
                            if ~isempty(alldbl5)
                                [allspecfit alldbl5] = refitbaddbl(dataminvib,allspecfit,alldbl5,lbds,ubds,[peak1 peak2 fwhm2 skewmax2],fwhmdbl,constrainblue,SQthr,mainvsvibfactor,scalewidth_wav,scalef_wav);
                                for spec=alldbl5 %estimate intensities of double-band spectra
                                    int1=sum(skewgaussian3([allspecfit(1,spec),allspecfit(2,spec),allspecfit(3,spec),allspecfit(4,spec)],wav));
                                    int2=sum(skewgaussian3([allspecfit(5,spec),allspecfit(6,spec),allspecfit(7,spec),allspecfit(8,spec)],wav));
                                    allint1(spec)=int1;
                                    allint2(spec)=int2;
                                end
                            end
                            if ~isempty(alldbl)
                                [allspecfit alldbl] = refitbaddbl(dataminvib,allspecfit,alldbl,olbds,oubds,[peak1 peak2 fwhm2 skewmax2],fwhmdbl,constrainblue,SQthr,mainvsvibfactor,scalewidth_wav,scalef_wav);
                                for spec=alldbl %estimate intensities of double-band spectra
                                    % check misfits
                                    if (allspecfit(5,spec)>0)
                                        if (allspecfit(1,spec)/allspecfit(5,spec)>6)&&(allspecfit(7,spec)<allspecfit(3,spec)+allspecfit(2,spec)+allspecfit(6,spec))
                                            allspecfit(:,spec)=[singlespecfit(:,spec)' 0 0 0 0];     % small red band fitted in vib tail or small blue band. 
                                        elseif (allspecfit(4,spec)<lbds(4)/3)&&(allspecfit(8,spec)>ubds(8)/2)&&(dblspecfit(2)+dblspecfit(6)>fwhmdouble*2)   % broad/dbl spectrum misfitted
                                            allspecfit(:,spec)=zeros(1,8); 
                                        elseif allspecfit(3,spec)<lbds(3)+.1
                                            if singlespecfit(3,spec)>peak1+fwhm1/2
                                                if allspecfit(7,spec)>peak1+fwhm1/2
                                                    allspecfit(:,spec)=zeros(1,8); 
                                                else
                                                    allspecfit(:,spec)=[allspecfit(5:8) 0 0 0 0];
                                                end
                                            else
                                                allspecfit(:,spec)=[singlespecfit(:,spec)' 0 0 0 0];
                                            end
                                        end
                                    end
                                    int1=sum(skewgaussian3([allspecfit(1,spec),allspecfit(2,spec),allspecfit(3,spec),allspecfit(4,spec)],wav));
                                    int2=sum(skewgaussian3([allspecfit(5,spec),allspecfit(6,spec),allspecfit(7,spec),allspecfit(8,spec)],wav));
                                    allint1(spec)=int1;
                                    allint2(spec)=int2;
                                end
                            end
                            if ~isempty([alldbl alldbl5])
                                doubles(i)=complex;
                            else
                                allspecfit(5:8,:)=-2;  % for the one figure
                            end
                        end
                     end
                    
              %%% 10. Average fitting parameters and jumps are calculated   %%%

                    allpeaks1=zeros(maxspecs,1);
                    allpeaks2=zeros(maxspecs,1);
                    allfwhm1=zeros(maxspecs,1);
                    allfwhm2=zeros(maxspecs,1);
                    allskew1=zeros(maxspecs,1);
                    allskew2=zeros(maxspecs,1);
                    allfi1=zeros(maxspecs,1);
                    alljumps1=zeros(maxspecs,1);
                    alljumps2=zeros(maxspecs,1);
                    alljumps12=zeros(maxspecs,1); %difference between peaks1 and peaks2
                    allvibampl=zeros(maxspecs,1);
                    allvibpeaks=zeros(maxspecs,1);
                    allvibfwhm=zeros(maxspecs,1);
                    allvibjumps=zeros(maxspecs,1);
                    allpeakjumps4vib=zeros(maxspecs,1); %jumps in peak1 to correlate with jumps in vib band
                    allrms=zeros(maxspecs,1);

                    allpeaks1(firstgood:lastgood)=allspecfit(3,firstgood:lastgood);
                    allpeaks2(firstgood:lastgood)=allspecfit(7,firstgood:lastgood);
                    allfwhm1(firstgood:lastgood)=allspecfit(2,firstgood:lastgood);
                    allfwhm2(firstgood:lastgood)=allspecfit(6,firstgood:lastgood);
                    allskew1(firstgood:lastgood)=allspecfit(4,firstgood:lastgood);
                    allskew2(firstgood:lastgood)=allspecfit(8,firstgood:lastgood);

                    %clean up fits - check why there are still misfits sometimes! 
    %                badspecs=(allpeaks2>0 & allpeaks2<peakmin+fwhmmin) | (allpeaks1>0 & allpeaks1<peakmin) | (allfwhm1>0 & allfwhm1<fwhmmin+.1) | allfwhm1>fwhmmax2-.1 | allskew1>=skewmax1-.01 | allskew2>=skewmax2-.01 | allint1<=0 | (allint2>0 & allint2<SNRdbl/2); %| (allfwhm2>0 & allfwhm2<allfwhm1)
    %                allpeaks1(badspecs)=0;
    %                allpeaks2(badspecs)=0;
    %                allfwhm1(badspecs)=0;
    %                allfwhm2(badspecs)=0;
    %                allskew1(badspecs)=0;
    %                allskew2(badspecs)=0;

                    %only values of fitted peaks
                    goodpeaksbool1=allint1>0 & allpeaks1>0 & (allpeaks2<=0 | allpeaks2>allpeaks1+fwhm1/2);
                    goodpeaksbool2=allint2>0 & allpeaks2>allpeaks1+fwhm1/2;
                    numgood1=length(goodpeaksbool1(goodpeaksbool1>0));
                    numgood2=length(goodpeaksbool2(goodpeaksbool2>0));

                    avgint1(i)=sum(allint1(goodpeaksbool1))/numgood1; %faster than built-in "mean"
                    avgfwhm1(i)=sum(allfwhm1(goodpeaksbool1))/numgood1;
                    avgpeaks1(i)=sum(allpeaks1(goodpeaksbool1))/numgood1;
                    avgskew1(i)=sum(allskew1(goodpeaksbool1))/numgood1;
                    stdint1(i)=std(allint1(goodpeaksbool1));
                    stdfwhm1(i)=std(allfwhm1(goodpeaksbool1));
                    stdpeaks1(i)=std(allpeaks1(goodpeaksbool1));
                    stdskew1(i)=std(allskew1(goodpeaksbool1));
                    if numgood2>0  % don't bother with other values (-2 or -1). Default = 0.
                        avgint2(i)=sum(allint2(goodpeaksbool2))/numgood2;
                        avgfwhm2(i)=sum(allfwhm2(goodpeaksbool2))/numgood2;
                        avgpeaks2(i)=sum(allpeaks2(goodpeaksbool2))/numgood2;
                        avgskew2(i)=sum(allskew2(goodpeaksbool2))/numgood2;
                        stdint2(i)=std(allint2(goodpeaksbool2));
                        stdfwhm2(i)=std(allfwhm2(goodpeaksbool2));
                        stdpeaks2(i)=std(allpeaks2(goodpeaksbool2));
                        stdskew2(i)=std(allskew2(goodpeaksbool2));
                    end

                    if stdfwhm2(i)>3  % arbitrary number: red fluctuation should be large enough for such a comparison
                        allfi1(firstgood:lastgood)=allint1(firstgood:lastgood)./(allint1(firstgood:lastgood)+allint2(firstgood:lastgood))*100;
                    end

                    %calculate spectral jumps
                    goodpeaks1=allpeaks1(goodpeaksbool1);
                    goodpeaks2=allpeaks2(goodpeaksbool2);
                    alljumps1(1:length(goodpeaks1)-1)=diff(goodpeaks1);
                    alljumps2(1:length(goodpeaks2)-1)=diff(goodpeaks2);
                    goodpeaks12=allpeaks2(allpeaks2>0)-allpeaks1(allpeaks2>0);
                    alljumps12(1:length(goodpeaks12))=goodpeaks12;

                    %calculate jumps in vib peak
                    if fitvib
                        allvibampl(1:lastgood)=allvibfit(1,:);
                        allvibfwhm(1:lastgood)=allvibfit(2,:);
                        allvibpeaks(1:lastgood)=allvibfit(3,:);
                        allgoodvibpeaks=allvibfit(3,allvibfit(3,:)>0);
                        allvibjumps(1:length(allgoodvibpeaks)-1)=diff(allgoodvibpeaks);
                        allgoodpeaks4vib=allpeaks1(allvibfit(3,:)>0);
                        allgoodpeaks4vib=allgoodpeaks4vib(allgoodpeaks4vib>0);
                        allpeakjumps4vib(1:length(allgoodpeaks4vib)-1)=diff(allgoodpeaks4vib);
                    end

                    %figure with all fits
                    l=ceil(sqrt(min(150,lastgood)*3/2));
                    w=ceil(sqrt(min(150,lastgood)*2/3));
                    m=ceil(lastgood/(l*w));
                    for k=1:m
                        h=figure;
                        set(h,'visible','off') 
                        set(gcf, 'Units', 'centimeters');
                        afFigurePosition = [2 2 42 22]; % [pos_x pos_y width_x width_y]
                        set(gcf, 'Position', afFigurePosition);
                        %axes('FontSize',floor(15/n)+5);
                        for j=firstgood+(k-1)*100:min(k*100,lastgood)
                            hh=subplot(l,w,j-(k-1)*100);
                            hold on;
                            h=plot(wav,mat(:,j+1),'k');
                            axis tight;
                            set(hh,'FontSize',floor(10/w)+5);
                            if allpeaks1(j)>0;
                                switch fittype
                                    case 1
                                        fit1=skewgaussian1(allspecfit(1:4,j),wav);
                                        fit2=skewgaussian1(allspecfit(5:8,j),wav);
                                    case 2
                                        fit1=skewgaussian2(allspecfit(1:4,j),wav);
                                        fit2=skewgaussian2(allspecfit(5:8,j),wav);
                                    case 3
                                        fit1=skewgaussian3(allspecfit(1:4,j),wav);
                                        fit2=skewgaussian3(allspecfit(5:8,j),wav);
                                end
                                fit3=gaussian(wav,allvibfit(1,j),allvibfit(2,j),allvibfit(3,j));
                                plot(wav,fit1,'b');
                                if allspecfit(5,j)>0
                                    plot(wav,fit2,'b');
                                end
                                if allvibfit(1,j)>0
                                    plot(wav,fit3,'g');
                                end
                                fit=fit1+fit2+fit3;
                                plot(wav,fit,'r');
                                allrms(j)=sqrt(sum((fit-mat(:,j+1)).^2)/length(mat));
                            end
                        end
                        if m>1
                            saveas(h,fullfile(writedir,['spec' int2str(i) '_' int2str(k) '.jpg']));
                            saveas(h,fullfile(writedir,['spec' int2str(i) '_' int2str(k) '.pdf']));
                        else
                            saveas(h,fullfile(writedir,['spec' int2str(i) '.jpg']));
                            saveas(h,fullfile(writedir,['spec' int2str(i) '.pdf']));
                        end
                        close all;
                    end
                    dlmwrite(fullfile(writedir,['specminbg' int2str(i)]),mat,' '); 
                    dlmwrite(fullfile(writedir,['dataspec' int2str(i)]),[allint1(firstgood:lastgood) allpeaks1(firstgood:lastgood) allfwhm1(firstgood:lastgood) allskew1(firstgood:lastgood) allint2(firstgood:lastgood) allpeaks2(firstgood:lastgood) allfwhm2(firstgood:lastgood) allskew2(firstgood:lastgood)],'\t');
                    if isempty(alldbl)
                        goods(i)=complex;
                    end
                    
                    allcpeaks1(:,i)=allpeaks1;
                    allcpeaks2(:,i)=allpeaks2;
                    allcfwhm1(:,i)=allfwhm1;
                    allcfwhm2(:,i)=allfwhm2;
                    allcskew1(:,i)=allskew1;
                    allcskew2(:,i)=allskew2;
                    allcfi1(:,i)=allfi1;
                    allcjumps1(:,i)=alljumps1;
                    allcjumps2(:,i)=alljumps2;
                    allcjumps12(:,i)=alljumps12;
                    allcvibampl(:,i)=allvibampl;
                    allcvibpeaks(:,i)=allvibpeaks;
                    allcvibfwhm(:,i)=allvibfwhm;
                    allcvibjumps(:,i)=allvibjumps;
                    allcpeakjumps4vib(:,i)=allpeakjumps4vib;
                    allcrms(:,i)=allrms;
                    allcint1(:,i)=allint1;
                    allcint2(:,i)=allint2;
                    allcSNR(:,i)=allSNR;
                end
            end
        end
    end
    disp(i);
end

%% Figures

if sum(goods+doubles)>0
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

 %% write to files
    % write all fit data to matrices
    dlmwrite(fullfile(writedir,'allints1'),allcint1,'\t');    
    dlmwrite(fullfile(writedir,'allpeaks1'),allcpeaks1,'\t');   
    dlmwrite(fullfile(writedir,'allfwhm1'),allcfwhm1,'\t');    
    dlmwrite(fullfile(writedir,'allskew1'),allcskew1,'\t');   
    dlmwrite(fullfile(writedir,'allints2'),allcint2,'\t');             
    dlmwrite(fullfile(writedir,'allpeaks2'),allcpeaks2,'\t');      
    dlmwrite(fullfile(writedir,'allfwhm2'),allcfwhm2,'\t');      
    dlmwrite(fullfile(writedir,'allskew2'),allcskew2,'\t');        
    dlmwrite(fullfile(writedir,'allskew2'),allcskew2,'\t');       
    dlmwrite(fullfile(writedir,'allfi1'),allcfi1,'\t');
    dlmwrite(fullfile(writedir,'avgspec'),[avgint1' avgpeaks1' avgfwhm1' avgskew1' avgint2' avgpeaks2' avgfwhm2' avgskew2'],'\t');
    dlmwrite(fullfile(writedir,'stdspec'),[stdint1' stdpeaks1' stdfwhm1' stdskew1' stdint2' stdpeaks2' stdfwhm2' stdskew2'],'\t');
    dlmwrite(fullfile(writedir,'allrms'),allcrms,'\t');

    %write good fit data to single column
    dlmwrite(fullfile(writedir,'col_allints1'),allcint1(allcpeaks1>0));    
    dlmwrite(fullfile(writedir,'col_allpeaks1'),allcpeaks1(allcpeaks1>0));    
    dlmwrite(fullfile(writedir,'col_allfwhm1'),allcfwhm1(allcpeaks1>0));    
    dlmwrite(fullfile(writedir,'col_allskew1'),allcskew1(allcpeaks1>0));   
    dlmwrite(fullfile(writedir,'col_allints2'),allcint2(allcpeaks2>0));            
    dlmwrite(fullfile(writedir,'col_allpeaks2'),allcpeaks2(allcpeaks2>0));
    dlmwrite(fullfile(writedir,'col_allfwhm2'),allcfwhm2(allcpeaks2>0));
    dlmwrite(fullfile(writedir,'col_allskew2'),allcskew2(allcpeaks2>0));
    dlmwrite(fullfile(writedir,'col_allfi1'),allcfi1(allcpeaks2>0));

    dlmwrite(fullfile(writedir,'intvspeak_blues'),[allcpeaks1(goodspecs2)' allcint1(goodspecs2)'],'\t');
    dlmwrite(fullfile(writedir,'fwhmvspeak_blues'),[allcpeaks1(goodspecs2)' allcfwhm1(goodspecs2)'],'\t');
    dlmwrite(fullfile(writedir,'skewvspeak_blues'),[allcpeaks1(goodspecs2)' allcskew1(goodspecs2)'],'\t');
    dlmwrite(fullfile(writedir,'intvsfwhm_blues'),[allcfwhm1(goodspecs2)' allcint1(goodspecs2)'],'\t');
    dlmwrite(fullfile(writedir,'intvspeak_reds'),[allcpeaks2(goodspecs2)' allcint2(goodspecs2)'],'\t');
    dlmwrite(fullfile(writedir,'fwhmvspeak_reds'),[allcpeaks2(goodspecs2)' allcfwhm2(goodspecs2)'],'\t');
    dlmwrite(fullfile(writedir,'skewvspeak_reds'),[allcpeaks2(goodspecs2)' allcskew2(goodspecs2)'],'\t');
    dlmwrite(fullfile(writedir,'intvsfwhm_reds'),[allcfwhm2(goodspecs2)' allcint2(goodspecs2)'],'\t');
    dlmwrite(fullfile(writedir,'intvspeak_singles'),[allcpeaks1(allcpeaks2==-2&goodspecs1)' allcint1(allcpeaks2==-2&goodspecs1)'],'\t');
    dlmwrite(fullfile(writedir,'fwhmvspeak_singles'),[allcpeaks1(allcpeaks2==-2&goodspecs1)' allcfwhm1(allcpeaks2==-2&goodspecs1)'],'\t');
    dlmwrite(fullfile(writedir,'skewvspeak_singles'),[allcpeaks1(allcpeaks2==-2&goodspecs1)' allcskew1(allcpeaks2==-2&goodspecs1)'],'\t');
    dlmwrite(fullfile(writedir,'intvsfwhm_singles'),[allcfwhm1(allcpeaks2==-2&goodspecs1)' allcint1(allcpeaks2==-2&goodspecs1)'],'\t');

    saveas(h1,fullfile(writedir,'FLP and jumps.jpg'));
    saveas(h2,fullfile(writedir,'Blue yield vs red FLP.jpg'));
    saveas(h3,fullfile(writedir,'spec par.jpg'));

    dlmwrite(fullfile(writedir,'histavgpeaks1'),[markavgpeak' histpeaks1'],'\t');
    dlmwrite(fullfile(writedir,'histavgpeaks2'),[markavgpeak' histpeaks2'],'\t');
    dlmwrite(fullfile(writedir,'histstdpeaks1'),[markstdpeak' histstdpeaks1'],'\t');
    dlmwrite(fullfile(writedir,'histstdpeaks2'),[markstdpeak' histstdpeaks2'],'\t');
    dlmwrite(fullfile(writedir,'histallpeaks1'),[markallpeaks' histallpeaks1'],'\t');
    dlmwrite(fullfile(writedir,'histallpeaks2'),[markallpeaks' histallpeaks2'],'\t');
    dlmwrite(fullfile(writedir,'histalljumps1'),[markjumps' jumps1'],'\t');    
    dlmwrite(fullfile(writedir,'histalljumps2'),[markjumps' jumps2'],'\t');
    dlmwrite(fullfile(writedir,'histalljumps12'),[markjumps' jumps12'],'\t');
end

deads=deads(deads>0); blues=blues(blues>0); doubles=doubles(doubles>0); broads=broads(broads>0); goods=goods(goods>0);
l1=length(deads); l2=length(blues); l3=length(doubles); l4=length(broads); l5=length(goods);
l=max(max(max(max(l1,l2),l3),l4),l5);
deads=[deads zeros(1,l-l1)];
blues=[blues zeros(1,l-l2)];
doubles=[doubles zeros(1,l-l3)];
broads=[broads zeros(1,l-l4)];
goods=[goods zeros(1,l-l5)];
dlmwrite(fullfile(writedir,'sortedpeaks'),'normal red blue noisy broad','');
dlmwrite(fullfile(writedir,'sortedpeaks'),[goods' doubles' blues' deads' broads'],'-append','delimiter','\t');

save(strcat(writedir,'\peakanalysis.mat'));
toc