function y = dipole_fit_backward(D,Mask,xx)

x = reshape(xx,size(Mask));

Ax = real(ifftn(D.* fftn(Mask.*x)));
y = Ax( Mask(:) == 0);
end