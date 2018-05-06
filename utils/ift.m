%
%	function im = ift(dat)
%
%	Function does a centered inverse 2DFT of the data:
%
%	im = fftshift(ifft2(fftshift(dat)));

function im = ift(dat)

im = fftshift(ifft2(fftshift(dat)));

