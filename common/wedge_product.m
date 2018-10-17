function wp = wedge_product(mesh, w1, w2)
% compute wedge product of two one-forms
w1_vector = vector_representation(mesh, w1);
w2_vector = vector_representation(mesh, w2);
f_normal = face_normal(mesh.face, mesh.vert);
wp = dot(cross(w1_vector, w2_vector), f_normal, 2);
end