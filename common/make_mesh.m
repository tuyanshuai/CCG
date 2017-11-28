function mesh = make_mesh(face, vert)

[edge,eif] = compute_edge(face);
[he,heif] = halfedge(face);

nf = size(face,1);
nv = size(vert,1);
ne = size(edge,1);
nh = size(he,1);

mesh.nf = nf;
mesh.nv = nv;
mesh.ne = ne;
mesh.nh = nh;

mesh.face = face;
mesh.vert = vert;
mesh.edge = edge;
mesh.eif = eif;
mesh.halfedge = he;
mesh.heif = heif;
