% Projection onto Dipole Fields (PDF)
%   [RDF shim] = PDF(iFreq, N_std, Mask,matrix_size,voxel_size, B0_dir,tol)
% 
%   output
%   RDF - the relative difference field, or local field
%   shim (optional) - the cropped background dipole distribution
% 
%   input
%   iFreq - the unwrapped field map
%   N_std - the noise standard deviation on the field map. (1 over SNR for single echo)
%   Mask - a binary 3D matrix denoting the Region Of Interest
%   matrix_size - the size of the 3D matrix
%   voxel_size - the size of the voxel in mm
%   B0_dir - the direction of the B0 field
%   tol(optional) - tolerance level
%
%   When using the code, please cite 
%   T. Liu et al. NMR Biomed 2011;24(9):1129-36
%   de Rochefort et al. MRM 2010;63(1):194-206
%
%   Created by Tian Liu in 2009
%   Modified by Tian Liu on 2011.02.01
%   Last modified by Tian Liu on 2013.07.24


function [RDF, back_PDF, shim] = PDF_xm(iFreq, Mask,matrix_size,voxel_size, B0_dir)

if (nargin<7)
    tol = 0.1;
end


% zero pad
matrix_size0 = matrix_size;
d1 = max(max(Mask,[],2),[],3);
d1first = find(d1,1,'first');
d1last = find(d1,1,'last');

d2 = max(max(Mask,[],1),[],3);
d2first = find(d2,1,'first');
d2last = find(d2,1,'last');
        
d3 = max(max(Mask,[],1),[],2);
d3first = find(d3,1,'first');
d3last = find(d3,1,'last');

matrix_size = [ floor((d1last - d1first+40)/2)*2,...
                floor((d2last - d2first+40)/2)*2,...
                floor((d3last - d3first+40)/2)*2];
iFreq=iFreq(d1first:d1last,d2first:d2last, d3first:d3last);

Mask=Mask(d1first:d1last,d2first:d2last, d3first:d3last);
padsize = [matrix_size(1)-size(iFreq,1) matrix_size(2)-size(iFreq,2) matrix_size(3)-size(iFreq,3)];
iFreq = padarray(iFreq, padsize, 0,'post');

Mask = padarray(Mask, padsize, 0,'post');


% generate the weighting
% W = 1./N_std;
% W(isinf(W)) =0;
% W = W.*(Mask>0);
% W_std = W;
% W_var = W.^2;
% W_var = ones(size(W_var));

%%%%%% start the PDF method %%%%%
D = dipole_kernel(matrix_size,voxel_size,B0_dir);

% generating the RHS vector in Eq. 6 in the PDF paper
p_temp = real(ifftn(D.*fftn(Mask.*(iFreq) )));
b = p_temp( Mask(:) == 0);

% set erance level and maximum iteration allowed
% E_noise_level = real(ifftn(D.*fftn(W_std.*ones(size(N_std)))));
itermax= 100;
A=@(xx)(dipole_term_bg(D,Mask,xx) );
cg_tol = 0;%tol*norm(E_noise_level( Mask(:) == 1))/norm(b(:));

[x res num_iter] = cgsolve(A, b, cg_tol, itermax, 1);
 
xx = zeros(size(D));
xx(Mask(:) == 0) = x(1:end);
xx(Mask(:) > 0) = 0;

% background dipole field
p_dipole = real(ifftn(D.*fftn(xx)));

p_final = (iFreq-p_dipole).*Mask;  
p_back = p_dipole.*Mask;

% remove zero pad
RDF = zeros(matrix_size0);
RDF(d1first:d1last,d2first:d2last, d3first:d3last) = ...
    p_final(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1);

back_PDF = zeros(matrix_size0);
back_PDF(d1first:d1last,d2first:d2last, d3first:d3last) = ...
    p_back(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1);

shim = zeros(matrix_size0);
shim(d1first:d1last,d2first:d2last, d3first:d3last) = ...
    xx(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1);

