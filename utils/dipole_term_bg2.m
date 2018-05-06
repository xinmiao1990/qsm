% Normal Equation of the Forward Calculation used in PDF
%   y = dipole_term(W,D,Mask,xx)
% 
%   output
%   y - a background field 
% 
%   input
%   W - noise covariance matrix
%   D - dipole kernel
%   xx - the background dipoles
%
%   Created by Tian Liu in 2009
%   Last modified by Tian Liu on 2013.07.24

function y = dipole_term_bg2(D,Masks,Maskl,xx)

% xx is a colum vector
x = zeros(size(D));
x(Maskl(:) == 0) = xx(1:end);

Ax = real(ifftn(D.* fftn(Masks.*real(ifftn(D.* fftn(x) )))));
y = Ax( Maskl(:) == 0);






