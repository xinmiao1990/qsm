function out = A1_v2(in, grad1)

out = cat(1, in(1)*grad1, in(2)*grad1, in(3)*grad1);

end