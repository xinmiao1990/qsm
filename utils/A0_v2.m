function out = A0_v2(in, grad0)

out = cat(1, in(1)*grad0, in(2)*grad0, in(3)*grad0);

end