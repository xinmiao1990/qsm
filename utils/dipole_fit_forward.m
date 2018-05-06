
function y = dipole_fit_forward(D,Mask,xx)

x = zeros(size(D));
x(Mask(:) == 0) = xx(1:end);

y = Mask.*real(ifftn(D.* fftn(x) ));
y = y(:);
end