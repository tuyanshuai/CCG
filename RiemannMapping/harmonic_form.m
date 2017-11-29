function dh = harmonic_form(mesh, w)
% given a closed one form w, compute the harmonic one form in the
% homological class of w
% find function h such that delta(w + dh) = 0, delta(dh) = -delta_w
% since d(w + dh) = dw + d(dh) = 0, so w+dh is harmonic form
delta_w = exterior_co_derivative(mesh,w,1);
L = laplace_beltrami(mesh);
L(1,1) = L(1,1)-1;
h = -L\delta_w;
dh = exterior_derivative(mesh,h,0);
% dh is closed harmonic one form
dh = w + dh;
