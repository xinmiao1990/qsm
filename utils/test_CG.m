
s1 = 1;
s2 = 2;
s3 = 3;

nx = 4;
ny = 4;
nz = 2;

grad = repmat((-nx/2:nx/2-1)',[1 ny*nz]);
grad = grad(:);
Zero = zeros(nx*ny*nz,1);

A = [grad, Zero, Zero; Zero, grad, Zero; Zero, Zero, grad];

slope = [s1;s2;s3];

y0 = A*slope;

y = y0 + randn(nx*ny*nz*3,1) * 0.1;

%% Direct inverse
% x = (A'*A)\(A'*y);

%% CG 
Q = eye(3)* (sum(grad.^2));
y = [sum(grad.*y(1:nx*ny*nz)); sum(grad.*y((1:nx*ny*nz)+nx*ny*nz));sum(grad.*y((1:nx*ny*nz)+2*nx*ny*nz))]; %A'*y
x = CG_solver(Q,y,[0;0;0], 20);