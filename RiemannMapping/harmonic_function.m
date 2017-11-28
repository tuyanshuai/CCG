function f = harmonic_function(mesh, f0)
% given initial value on boundary (specified by f0), compute corresponding
% harmonic function

L = laplace_beltrami(mesh);
f = f0;
for i = 1:size(f0,2)
    ind = isnan(f0(:,i)) | isinf(f0(:,i));
    L2 = L(ind,ind);
    b = -L(ind,~ind)*f0(~ind,i);
    f(ind,i) = L2\b;
end
