function uv = annulus_riemann_map(mesh)
% compute riemann map from annulus surface to canonical annulus in unit
% disk

% compute harmonic function f with fixed value on two boundaries of annulus
% f(b1) = 0, f(b2) = 1
L = laplace_beltrami(mesh.face,mesh.vert);
bds = boundary(mesh.face);
b1 = bds{1};
b2 = bds{2};
f = zeros(mesh.nv,1);
f(b1) = 0;
f(b2) = 1;
ind = true(mesh.nv,1);
ind(b1) = false;
ind(b2) = false;
L2 = L(ind,ind);
b = -L(ind,~ind)*f(~ind);
f(ind) = L2\b;
% find a shortest path between two boundaries b1,b2
path = shortest_path(mesh,[b1(1),b2(1)]);
