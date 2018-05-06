function imlow = HP_hamm(L, imsize, im)

w = hamming(L);

m = w(:)*w(:).'; % Create 2D window

mask = zeros(imsize);

mask((-L/2:(L/2-1))+imsize(1)/2,(-L/2:(L/2-1))+imsize(2)/2,:)=repmat(m,[1 1 imsize(3)]);

imlow = real(ift3d(fft3d(im).*(1-mask)));

end