function iv = center(mesh)
% find the (rough) central vertex of input mesh, return index of central
% vertex. input mesh must be simply-connected with one boundary
L = laplace_beltrami(mesh);
% the smallest eigenvalue is 0, eigenvector is constant vector
% second smallest eigenvector (Fiedler vector)'s signs partition mesh into
% two approximately even parts, the cut is approximately the central
% foliation along that direction.
% third smallest eigenvector has similar property (as observed, not 
% confirmed), the cut foliates in orthogonal direction.
% thus the intersection of second and third smallest eigenvector's central
% foliation gives an rough center of mesh.
[V,~] = eigs(L,3,0);
[~,iv] = min(dot(V,V,2));
