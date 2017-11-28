function uv = riemann_map(mesh)
% compute riemann map from simply-connected surface to unit disk

% first find a central triangle and remove it, remaining mesh becomes 
% topological annulus, with same vertex set
iv = center(mesh);
vfr = vert_face_ring(mesh,iv);
ivf = vfr{1};
ind = true(mesh.nf,1);
ind(ivf(1)) = false;
face2 = mesh.face(ind,:);
mesh2 = make_mesh(face2,mesh.vert);
uv = annulus_riemann_map(mesh2);
