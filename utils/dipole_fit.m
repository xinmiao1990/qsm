function y = dipole_fit(D,Mask,xx)

% xx is a colum vector
x = zeros(size(D));
x(Mask(:) == 0) = xx(1:end);

Ax = real(ifftn(D.* fftn(Mask.*real(ifftn(D.* fftn(x) )))));
y = Ax( Mask(:) == 0);
end

