
function QSM = admm_qsm_withMask_withZeroCSF(phase_use, N, mask_use, mask_CSF, spatial_res, lambda1, lambda2)

% This function solves this optimization problem:
% || M Fh D F X - M f ||^2_2 + lambda1 * || TV(X) ||^1 + ...
% lambda2 * || M_csf X ||^2_2 + ...
% lambda3 * || M_notissue X ||^2_2
%
% lambda1: regularization weighting for the total variation constraint
% lambda2: regularization weighting for the zero CSF constraint
% lambda3: regularization weighting for the finite spatial support
% constraint
% 
% Algorithm: ADMM
% || M Fh D F X - M f ||^2_2 + ...
% lambda1 * || zdx ||^1 + lambda1 * || zdy ||^1 + lambda1 * || zdz ||^1 +
% mu * || FDx(X) - zdx ||^2_2 + mu * || FDy(X) - zdy ||^2_2 + mu * || FDz(X) - zdz ||^2_2
% lambda2 * || Mcsf X ||^2_2 + ...
% lambda3 * || Mnon-tissue X ||^2_2
% 
% Author: Xin Miao
% Email: xinm@usc.edu
%% create dipole kernel 
kernel = dipole_kernel(N, spatial_res); 

%% create finite difference operator (FDx) and its Hermician form in Fourier space
[kx, ky, kz] = ndgrid(0:N(1)-1, 0:N(2)-1, 0:N(3)-1);
FDx = 1 - exp(2i .* pi .* kx / N(1));
FDy = 1 - exp(2i .* pi .* ky / N(2));
FDz = 1 - exp(2i .* pi .* kz / N(3));

FDxH = conj(FDx);
FDyH = conj(FDy);
FDzH = conj(FDz);

FD2 = FDxH .* FDx + FDyH .* FDy + FDzH .* FDz; clear kx ky kz;

%% ADMM recon
mu = 1e-2;              % gradient consistency
outer_iter = 50;
tol_update = 1;

z_dx = zeros(N, 'single');
z_dy = zeros(N, 'single');
z_dz = zeros(N, 'single');

e_dx = zeros(N, 'single');
e_dy = zeros(N, 'single');
e_dz = zeros(N, 'single');

x = zeros(N, 'single');

tic
for t = 1:outer_iter
    % update x : susceptibility estimate
    tx = FDxH .* fftn(z_dx - e_dx);
    ty = FDyH .* fftn(z_dy - e_dy);
    tz = FDzH .* fftn(z_dz - e_dz);
    
    x_prev = x;
    
    % solve L2 terms
    kspace = fftn(mask_use.*phase_use);
    b = ifftn(conj(kernel) .* kspace);
    b = b+ ifftn(mu * (tx + ty + tz));
    b = b(:);

    A=@(xx)( L2_terms(kernel, mask_use, mask_CSF, FD2, mu, lambda2, xx) );
    cg_tol = 1e-4;
    inner_iter = 20;
    [x,~,~] = cgsolve(A, b, cg_tol, inner_iter, 1);
    x = reshape(x,N);
    
    %
    x_update = 100 * norm(x(:)-x_prev(:)) / norm(x(:));
    disp(['Iter: ', num2str(t), '   Update: ', num2str(x_update)])
    
    if x_update < tol_update
        break
    end
    
    if t < outer_iter
        % update z : gradient varible
        Fx = fftn(x);
        x_dx = ifftn(FDx .* Fx);
        x_dy = ifftn(FDy .* Fx);
        x_dz = ifftn(FDz .* Fx);

        z_dx = max(abs(x_dx + e_dx) - lambda1 / mu, 0) .* sign(x_dx + e_dx);
        z_dy = max(abs(x_dy + e_dy) - lambda1 / mu, 0) .* sign(x_dy + e_dy);
        z_dz = max(abs(x_dz + e_dz) - lambda1 / mu, 0) .* sign(x_dz + e_dz);

        % update s : Lagrange multiplier
        e_dx = e_dx + x_dx - z_dx;
        e_dy = e_dy + x_dy - z_dy;            
        e_dz = e_dz + x_dz - z_dz;    
        view3dgui(real(x),[0.45 0.45 1.3]);
        
        if(mod(t,2)==0)
            save(sprintf('QSM_iter%d.mat',t),'x');
            disp(t);
        end
    end
    
end
toc

QSM = real(x);

end

function kernel = dipole_kernel(N, spatial_res)
    [ky,kx,kz] = meshgrid(-N(2)/2:N(2)/2-1, -N(1)/2:N(1)/2-1, -N(3)/2:N(3)/2-1);

    kx = (kx / max(abs(kx(:)))) / spatial_res(1);
    ky = (ky / max(abs(ky(:)))) / spatial_res(2);
    kz = (kz / max(abs(kz(:)))) / spatial_res(3);

    k2 = kx.^2 + ky.^2 + kz.^2;

    R_tot = eye(3);

    kernel = fftshift( 1/3 - (kx * R_tot(3,1) + ky * R_tot(3,2) + kz * R_tot(3,3)).^2 ./ (k2 + eps) ); 
end

function y = L2_terms(kernel, mask_use, mask_CSF, FD2, mu, lambda2, xx)
    [Nx, Ny, Nz] = size(mask_use);
    x = reshape(xx, [Nx, Ny, Nz]);

    y1 = ifftn( kernel .* fftn(x) );
    y1 = y1 .* mask_use;
    y1 = ifftn( conj(kernel) .* fftn(y1) );

    if mu~=0
        Fx = fftn(x);
        y2 = mu *ifftn( FD2 .* Fx );
    else 
        y2=0;
    end
    
    y3 = lambda2 * (mask_CSF .* x);

    y = y1 + y2 + y3;
    y = y(:);
end
