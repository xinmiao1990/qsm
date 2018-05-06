function iField = eddy_current_cor(iField,Mask,Mask_large)

% Mask is smaller than Mask_large, Mask is for fitting and correction is applied to the range of Mask_large 
iField = permute(iField,[2 1 3 4]);
Mask = permute(Mask,[2 1 3]); 
Mask_large = permute(Mask_large,[2 1 3]); 

mask = Mask; 
mask_large = Mask_large; 

Mask = repmat(Mask,[1 1 1 3]);
Mask_large = repmat(Mask_large,[1 1 1 3]);

nx = size(iField,1);
ny = size(iField,2);
nz = size(iField,3);

% apply mask

y1 = angle(iField(:,:,:,1).*conj(iField(:,:,:,2)));
for i =1:nz
    y1(:,:,i) = medfilt2(y1(:,:,i));
end

y2 = angle(iField(:,:,:,2).*conj(iField(:,:,:,3)));
for i =1:nz
    y2(:,:,i) = medfilt2(y2(:,:,i));
end

y3 = angle(iField(:,:,:,3).*conj(iField(:,:,:,4)));
for i =1:nz
    y3(:,:,i) = medfilt2(y3(:,:,i));
end

b=cat(4,y1,y2,y3);

%---------
line1 = b(:,nx/2,40,1);
line2 = b(:,nx/2,40,2);
line3 = b(:,nx/2,40,3);
figure(1);
subplot(1,3,1);plot(line1-line2);
subplot(1,3,2);plot(line1-line3);
subplot(1,3,3);plot(line2-line3);

%---------

b = b.*Mask;

b = b(Mask==1);

N = length(b)/3;

clear y1 y2 y3 y1l y2l y3l;

grad1 = repmat(((-nx/2:nx/2-1)/nx)',[1 ny nz])/sqrt(nx*ny*nz);
grad1 = grad1(mask==1);

grad0 = ones(nx,ny,nz)/sqrt(nx*ny*nz);
grad0 = grad0(mask==1);


%% Direct inverse
% x = (A'*A)\(A'*y);
phi0 = zeros(3,1);
phi1 = zeros(3,1);
x = zeros(N,1);

% A1 operator:
A1 = @ (x) A1_v2(x, grad1);

% A0 operator:
A0 = @ (x) A0_v2(x, grad0);

% A operator:
A = @ (x) repmat(x,[3 1]);

max_iter = 1;

for iter = 1:max_iter
%% CG 1: solve phi1: A1*phi1 = s - A0*phi0 - A*x
%A'A:
Q1 = eye(3)* (sum(grad1(:).^2));

% y
y = b(:) - A0(phi0) - A(x);

%A'y
y = reshape(y,[N,3]);
Y = y.*repmat(grad1,[1 3]); %s: nx,ny,nz,3
clear y;
Y=Y';
Y = sum(Y,2); 

phi1 = CG_solver(Q1,Y,phi1, 20);

%% CG 2: solve phi0: A0*phi0 = s - A1*phi1 - A*x
%A'A:
Q0 = eye(3)* (sum(grad0(:).^2));

% y
y = b(:) - A1(phi1) - A(x);

%A'y
y = reshape(y,[N,3]);
Y = y.*repmat(grad0,[1 3]); %s: nx,ny,nz,3
clear y;
Y=Y';
Y = sum(Y,2); 

phi0 = CG_solver(Q0,Y,phi0, 20);

%% CG 3: solve x: A*x = s - A1*phi1 - A0*phi0

% y
y = b(:) - A1(phi1) - A0(phi0);
y = reshape(y,[N,3]);
x_old = x;
x = sum(y,2)/3;
x = x(:);

Diff = x - x_old;

disp('Error:')
sum(Diff.^2)

error = A1(phi1) + A0(phi0);
error = embed(error, Mask);
% figure;imshow(error(:,:,40,1),[-1.5 1.5])
% figure;imshow(error(:,:,40,2),[-1.5 1.5])
% figure;imshow(error(:,:,40,3),[-1.5 1.5])

line1 = error(:,nx/2,40,1);
line2 = error(:,nx/2,40,2);
line3 = error(:,nx/2,40,3);
figure(1);
subplot(1,3,1);hold on; plot(line1-line2);
subplot(1,3,2);hold on; plot(line1-line3);
subplot(1,3,3);hold on; plot(line2-line3);
drawnow;

end

%% apply to larger range
grad1 = repmat(((-nx/2:nx/2-1)/nx)',[1 ny nz])/sqrt(nx*ny*nz);
grad1 = grad1(mask_large==1);

grad0 = ones(nx,ny,nz)/sqrt(nx*ny*nz);
grad0 = grad0(mask_large==1);

error = A1_v2(phi1,grad1) + A0_v2(phi0,grad0);
error = embed(error, Mask_large);

iField(:,:,:,2)=iField(:,:,:,2).*exp(1i*error(:,:,:,1));
iField(:,:,:,3)=iField(:,:,:,3).*exp(1i*error(:,:,:,2)).*exp(1i*error(:,:,:,1));
iField(:,:,:,4)=iField(:,:,:,4).*exp(1i*error(:,:,:,3)).*exp(1i*error(:,:,:,2)).*exp(1i*error(:,:,:,1));

y1 = angle(iField(:,:,:,1).*conj(iField(:,:,:,2)));
y2 = angle(iField(:,:,:,2).*conj(iField(:,:,:,3)));
y3 = angle(iField(:,:,:,3).*conj(iField(:,:,:,4)));

figure(1);
subplot(1,3,1);hold on; plot(y1(:,nx/2,40)-y2(:,nx/2,40));
subplot(1,3,2);hold on; plot(y1(:,nx/2,40)-y3(:,nx/2,40));
subplot(1,3,3);hold on; plot(y2(:,nx/2,40)-y3(:,nx/2,40));
drawnow;
iField = permute(iField,[2 1 3 4]);
end