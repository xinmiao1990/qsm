
function QSM= admm_qsm(phase_use, N,mask_use, spatial_res,lambda)


% N = matrix_size; clear matrix_size;
% mask_use = Mask; clear Mask;
% phase_use = RDF; clear RDF iFreq iFreq_raw;
% spatial_res = voxel_size; clear voxel_size;


%% create dipole kernel and noisy phase

[ky,kx,kz] = meshgrid(-N(2)/2:N(2)/2-1, -N(1)/2:N(1)/2-1, -N(3)/2:N(3)/2-1);

kx = (kx / max(abs(kx(:)))) / spatial_res(1);
ky = (ky / max(abs(ky(:)))) / spatial_res(2);
kz = (kz / max(abs(kz(:)))) / spatial_res(3);

k2 = kx.^2 + ky.^2 + kz.^2;

R_tot = eye(3);

kernel = fftshift( 1/3 - (kx * R_tot(3,1) + ky * R_tot(3,2) + kz * R_tot(3,3)).^2 ./ (k2 + eps) );    
 

%% Closed-form L2 recon

[kx, ky, kz] = ndgrid(0:N(1)-1, 0:N(2)-1, 0:N(3)-1);
Ex = 1 - exp(2i .* pi .* kx / N(1));
Ey = 1 - exp(2i .* pi .* ky / N(2));
Ez = 1 - exp(2i .* pi .* kz / N(3));

Ext = conj(Ex);
Eyt = conj(Ey);
Ezt = conj(Ez);

E2 = Ext .* Ex + Eyt .* Ey + Ezt .* Ez;
K2 = abs(kernel).^2;

%% TV ADMM recon

mu = 1e-2;              % gradient consistency
% lambda = 4e-4;          % gradient L1 penalty
 
num_iter = 50;
tol_update = 1;

z_dx = zeros(N, 'single');
z_dy = zeros(N, 'single');
z_dz = zeros(N, 'single');

s_dx = zeros(N, 'single');
s_dy = zeros(N, 'single');
s_dz = zeros(N, 'single');

x = zeros(N, 'single');

kspace = fftn(phase_use);
Dt_kspace = conj(kernel) .* kspace;


tic
for t = 1:num_iter
    % update x : susceptibility estimate
    tx = Ext .* fftn(z_dx - s_dx);
    ty = Eyt .* fftn(z_dy - s_dy);
    tz = Ezt .* fftn(z_dz - s_dz);
    
    x_prev = x;
    x = ifftn( (mu * (tx + ty + tz) + Dt_kspace) ./ (eps + K2 + mu * E2) );

    x_update = 100 * norm(x(:)-x_prev(:)) / norm(x(:));
    disp(['Iter: ', num2str(t), '   Update: ', num2str(x_update)])
    
    if x_update < tol_update
        break
    end
    
    if t < num_iter
        % update z : gradient varible
        Fx = fftn(x);
        x_dx = ifftn(Ex .* Fx);
        x_dy = ifftn(Ey .* Fx);
        x_dz = ifftn(Ez .* Fx);

        z_dx = max(abs(x_dx + s_dx) - lambda / mu, 0) .* sign(x_dx + s_dx);
        z_dy = max(abs(x_dy + s_dy) - lambda / mu, 0) .* sign(x_dy + s_dy);
        z_dz = max(abs(x_dz + s_dz) - lambda / mu, 0) .* sign(x_dz + s_dz);

        % update s : Lagrange multiplier
        s_dx = s_dx + x_dx - z_dx;
        s_dy = s_dy + x_dy - z_dy;            
        s_dz = s_dz + x_dz - z_dz;            
    end
end
toc

QSM=real(x) .* mask_use;

end