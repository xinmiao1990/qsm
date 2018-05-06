
function QSM = admm_qsm_withMask_withZeroCSF_v2(phase_use, N, mask_use, mask_tissue, spatial_res, lambda1)

% This function solves this optimization problem:
% || M Fh D F E X - M f ||^2_2 + lambda1 * || TV(E X) ||^1 
%
% E: embeding X into regions that are non-CSF tissue
% lambda1: regularization weighting for the total variation constraint
% lambda2: regularization weighting for the zero CSF constraint
% lambda3: regularization weighting for the finite spatial support
% constraint
% 
% Algorithm: ADMM
% || M Fh D F E X - M f ||^2_2 + ...
% lambda1 * || zdx ||^1 + lambda1 * || zdy ||^1 + lambda1 * || zdz ||^1 +
% mu * || FDx(V) - zdx - edx ||^2_2 + mu * || FDy(V) - zdy - edy ||^2_2 
% + mu * || FDz(V) - zdz - edz ||^2_2
% + mu * || E X - V - ev ||^2_2
%
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

% set initial value
x = zeros([sum(mask_tissue(:)), 1], 'single');

z_dx = zeros(N, 'single');
z_dy = zeros(N, 'single');
z_dz = zeros(N, 'single');
v = zeros(N, 'single');

e_dx = zeros(N, 'single');
e_dy = zeros(N, 'single');
e_dz = zeros(N, 'single');
e_v = zeros(N, 'single');

tic
for t = 1:outer_iter
    % update x : susceptibility estimate    
    x_prev = x;
    
    % solve for x
    kspace = fftn(mask_use.*phase_use);
    b = ifftn(conj(kernel) .* kspace);
    b = b + mu * (v + e_v);
    b = b(mask_tissue(:));

    A=@(xx)( L2_terms_x(kernel, mask_use, mask_tissue, mu, xx) );
    cg_tol = 1e-3;
    inner_iter = 20;
    [x,~,~] = cgsolve(A, b, cg_tol, inner_iter, 1);
    
    %
    x_update = 100 * norm(x(:)-x_prev(:)) / norm(x(:));
    disp(['Iter: ', num2str(t), '   Update: ', num2str(x_update)])
    
    if x_update < tol_update
        break
    end
    
    if t < outer_iter
        % solve for v
        tx = FDxH .* fftn(z_dx + e_dx);
        ty = FDyH .* fftn(z_dy + e_dy);
        tz = FDzH .* fftn(z_dz + e_dz);
        v = ifftn( ( tx + ty + tz + fftn(embed(x,mask_tissue)) ) ./ (FD2+1) );
    
        % solve for z
        Fv = fftn(v);
        v_dx = ifftn(FDx .* Fv);
        v_dy = ifftn(FDy .* Fv);
        v_dz = ifftn(FDz .* Fv);

        z_dx = max(abs(v_dx - e_dx) - lambda1 / mu, 0) .* sign(v_dx - e_dx);
        z_dy = max(abs(v_dy - e_dy) - lambda1 / mu, 0) .* sign(v_dy - e_dy);
        z_dz = max(abs(v_dz - e_dz) - lambda1 / mu, 0) .* sign(v_dz - e_dz);

        % update s : Lagrange multiplier
        e_dx = e_dx - (v_dx - z_dx);
        e_dy = e_dy - (v_dy - z_dy);            
        e_dz = e_dz - (v_dz - z_dz);  
        e_v = e_v - ( embed(x,mask_tissue) - v );
        QSM = embed(real(x),mask_tissue);
        view3dgui( QSM,[0.45 0.45 1.3]);
        
        if(mod(t,2)==0)
            save(sprintf('QSM_iter%d.mat',t),'QSM'); clear QSM;
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

function y = L2_terms_x(kernel, mask_use, mask_tissue, mu, xx)
    x = embed(xx,mask_tissue);

    y1 = ifftn( kernel .* fftn(x) );
    y1 = y1 .* mask_use;
    y1 = ifftn( conj(kernel) .* fftn(y1) );

    y = y1 + mu * x;
    y = y(mask_tissue(:));
end

