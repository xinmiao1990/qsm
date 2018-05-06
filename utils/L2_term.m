function y = L2_term(kernel,mask_use,E2,mu, xx)
    [Nx,Ny,Nz] = size(mask_use);
    x = reshape(xx,[Nx,Ny,Nz]);

    y1 = ifftn(kernel.*fftn(x));
    y1 = y1.*mask_use;
    y1 = ifftn(conj(kernel).*fftn(y1));

    if mu~=0
        Fx = fftn(x);
        y2 = mu*ifftn(E2.*Fx);
    else 
        y2=0;
    end

    y = y1+y2;
    y = y(:);
end