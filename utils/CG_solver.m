function x = CG_solver(Q,y,x0, max_iter)

% initial image
if (nargin<3)
    max_iter = 20;
end

x = x0;
% r0 = A' * (y - A * x); % residual
r0 = y - Q * x; % residual

% p0 = A' * (y - A * x); % the negative of the gradient of f
p0 = y - Q * x;
n_iter = 0;


r = r0;
p = p0;
err = 1;


while (1)
%     q = A' * (A *p);
    q = Q *p;
    alpha = r(:)' * r(:) / (p(:)' * q(:) );
    x = x + alpha * p;
    r_new = r - alpha * q;
    
    if sqrt(r_new(:)' * r_new(:))<1e-10
        break;
    end
    
    beta = r_new(:)' * r_new(:) / (r(:)' * r(:));
    p = r_new + beta * p;
    
    r = r_new;
    err_old = err;
    err= r(:)' * r(:) / (r0(:)' * r0(:));
    n_iter = n_iter + 1;
    

    fprintf('\n iteration #%d, alpha = %f, beta = %f, error = %f', n_iter, alpha, beta, abs(err));

    if (n_iter > max_iter) % || (abs(err) > abs (err_old))
        break;
    end
end