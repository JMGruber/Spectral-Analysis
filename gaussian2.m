function curve = gaussian2(v,x)

%  Usage:  curve = gaussian(v,x);
%     x is a vector of x data values
%     v = 1x3 matrix where
%     v(1)=A is the amplitude of the gaussian
%     v(2)=FWHM is the full width at half maximum
%     v(3)=Offset is the center of the gaussian
%     curve is a vector of y values corresponding 
%         to the provided x values.

curve=v(1).*exp((-4*log(2)/v(2).^2).*(x-v(3)).^2);


