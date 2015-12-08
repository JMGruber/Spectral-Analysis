function curve = dblskewgaussian2(bb,x)
%  Usage:  curve = dblskewgaussian2([A1,FWHM1,Offset1,A2,FWHM2,Offset2],x);

A1=bb(1); FWHM1=bb(2); Offset1=bb(3); Skew1=bb(4); A2=bb(5); FWHM2=bb(6); Offset2=bb(7); Skew2=bb(8);

if Skew1==0
    Skew1=0.001;
elseif Skew2==0
    Skew2=0.001;
else
    sigma1=Skew1*FWHM1/sinh(Skew1);
    sigma2=Skew2*FWHM2/sinh(Skew2);
    %curve=A1.*exp((-log(2)/Skew1^2).*(log(1+2*Skew1/FWHM1*(x-Offset1))).^2)+A2.*exp((-log(2)/Skew2^2).*(log(1+2*Skew2/FWHM2*(x-Offset2))).^2);
    curve1=A1.*exp((-log(2)/Skew1^2).*(log(1+2*Skew1/sigma1*(x-Offset1))).^2);
    curve2=A2.*exp((-log(2)/Skew2^2).*(log(1+2*Skew2/sigma2*(x-Offset2))).^2);
    if Skew1>0
        curve11=curve1(x>Offset1-sigma1/2/Skew1);
        curve1=[zeros(1,length(curve1)-length(curve11)) curve11']';
    else
        curve11=curve1(x<Offset1-sigma1/2/Skew1);
        curve1=[curve11' zeros(1,length(curve1)-length(curve11))]';
    end
    if Skew2>0
        curve22=curve2(x>Offset2-sigma2/2/Skew2);
        nullen=zeros(1,length(curve2)-length(curve22));
        if isempty(nullen)
            curve2=curve22;
        else
            curve2=[nullen curve22']';
        end
    else
        curve22=curve2(x<Offset2-sigma2/2/Skew2);
        curve2=[curve22' zeros(1,length(curve2)-length(curve22))]';
    end
    curve=curve1+curve2;
end