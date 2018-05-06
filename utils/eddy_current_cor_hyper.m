function iField = eddy_current_cor_hyper(iField)

% Mask is manually chosen, region within Mask should have no phase wraps 

% Mask is for fitting and correction is applied to the range of Mask_large 
iField = permute(iField,[2 1 3 4]);
nx = size(iField,1);
ny = size(iField,2);
nz = size(iField,3);

%---------
center=30;
Mask=zeros(nx,ny,nz);
Mask(nx/2+1+(-40:39)+center,ny/2+1+(-60:59),41:60)=1;
Mask = logical(Mask);
mask = Mask; 

Mask = repmat(Mask,[1 1 1 3]);

% apply mask

y1 = angle(iField(:,:,:,1).*conj(iField(:,:,:,2))).*mask;
% for i =1:nz
%     y1(:,:,i) = medfilt2(y1(:,:,i));
% end
figure(1);montage(reshape(y1,[nx ny 1 nz]),'DisplayRange',[]);

y2 = angle(iField(:,:,:,2).*conj(iField(:,:,:,3))).*mask;
% for i =1:nz
%     y2(:,:,i) = medfilt2(y2(:,:,i));
% end
figure(2);montage(reshape(y2,[nx ny 1 nz]),'DisplayRange',[]);

y3 = angle(iField(:,:,:,3).*conj(iField(:,:,:,4))).*mask;
% for i =1:nz
%     y3(:,:,i) = medfilt2(y3(:,:,i));
% end
figure(3);montage(reshape(y3,[nx ny 1 nz]),'DisplayRange',[]);



b=cat(4,y1,y2,y3);

%---------
line1 = b(:,nx/2,50,1);
line2 = b(:,nx/2,50,2);
line3 = b(:,nx/2,50,3);
figure(4);
subplot(2,3,1);hold on; plot(line1-line2,'r');
subplot(2,3,2);hold on; plot(line1-line3,'r');
subplot(2,3,3);hold on; plot(line2-line3,'r');
subplot(2,3,4);hold on; plot(line1,'r');
subplot(2,3,5);hold on; plot(line2,'r');
subplot(2,3,6);hold on; plot(line3,'r');

keyboard;
%---------

b = b.*Mask;

b = b(Mask==1);

N = length(b)/3;

clear y1 y2 y3 y1l y2l y3l;

grad1 = repmat(((-nx/2:nx/2-1)-center)',[1 ny nz]);
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

model = repmat(x,[3,1])+A1(phi1) + A0(phi0);
model = embed(model, Mask);


line1 = model(:,nx/2,50,1);
line2 = model(:,nx/2,50,2);
line3 = model(:,nx/2,50,3);
figure(4);
subplot(2,3,1);hold on; plot(line1-line2,'b');
subplot(2,3,2);hold on; plot(line1-line3,'b');
subplot(2,3,3);hold on; plot(line2-line3,'b');
subplot(2,3,4);hold on; plot(line1,'b');
subplot(2,3,5);hold on; plot(line2,'b');
subplot(2,3,6);hold on; plot(line3,'b');
drawnow;

end

%% apply to larger range
grad1 = repmat(((-nx/2:nx/2-1)-center)',[1 ny nz]);
grad1 = grad1(:);

grad0 = ones(nx,ny,nz)/sqrt(nx*ny*nz);
grad0 = grad0(:);

error = A1_v2(phi1,grad1) + A0_v2(phi0,grad0);
error = reshape(error,[nx,ny,nz,3]);

iField(:,:,:,2)=iField(:,:,:,2).*exp(1i*error(:,:,:,1));
iField(:,:,:,3)=iField(:,:,:,3).*exp(1i*error(:,:,:,2)).*exp(1i*error(:,:,:,1));
iField(:,:,:,4)=iField(:,:,:,4).*exp(1i*error(:,:,:,3)).*exp(1i*error(:,:,:,2)).*exp(1i*error(:,:,:,1));

y1 = angle(iField(:,:,:,1).*conj(iField(:,:,:,2)));
y2 = angle(iField(:,:,:,2).*conj(iField(:,:,:,3)));
y3 = angle(iField(:,:,:,3).*conj(iField(:,:,:,4)));

figure(4);
subplot(2,3,1);hold on; plot(y1(:,nx/2,50)-y2(:,nx/2,50),'g');
subplot(2,3,2);hold on; plot(y1(:,nx/2,50)-y3(:,nx/2,50),'g');
subplot(2,3,3);hold on; plot(y2(:,nx/2,50)-y3(:,nx/2,50),'g');

subplot(2,3,4);hold on; plot(y1(:,nx/2,50),'g');
subplot(2,3,5);hold on; plot(y2(:,nx/2,50),'g');
subplot(2,3,6);hold on; plot(y3(:,nx/2,50),'g');
drawnow;
iField = permute(iField,[2 1 3 4]);
end