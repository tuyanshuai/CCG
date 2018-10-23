function vf = vector_representation(mesh, w)
% vector valued 2-form of closed 1-form w

edge = mesh.edge;
vert = mesh.vert;
face = mesh.face;
normal = face_normal(face, vert);

nv = mesh.nv;
es = sparse(edge(:,1),edge(:,2),w,nv,nv);
es = es - conj(es');

vf = - full(diag(es(face(:, 1), face(:, 2)))) .* cross(vert(face(:, 3), :), normal) ...
    - full(diag(es(face(:, 3), face(:, 1)))) .* cross(vert(face(:, 2), :), normal) ...
    - full(diag(es(face(:, 2), face(:, 3)))) .* cross(vert(face(:, 1), :), normal);

end
