function iv = center(mesh)
% find the (rough) central vertex of input mesh, return index of central
% vertex. input mesh must be simply-connected with one boundary
L = laplace_beltrami(mesh.face, mesh.vert);
[V,~] = eigs(L,3,0);
[~,iv] = min(dot(V,V,2));
