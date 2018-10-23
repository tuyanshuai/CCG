function cw = conjugate_harmonic_one_form(mesh, w, hf)
% compute conjugate harmonic 1-form of w

if ~exist('hf','var') || isempty(hf)
    if ~exist('cf', 'var')
        cf = characteristic_one_form(mesh);
    end
    hf = cell(size(cf));
    for i = 1:length(hf)
        hf{i} = harmonic_form(mesh, cf{i});
    end
end

hf_vector = cell(size(hf));
for i = 1:length(hf)
    hf_vector{i} = vector_representation(mesh, hf{i});
end

fnormal = face_normal(mesh.face, mesh.vert);
farea = face_area(mesh.face, mesh.vert);

w_vector = vector_representation(mesh, w);
dual_w = cross(fnormal, w_vector);

% build linear system
s = ones(length(hf), 1);
for i = 1:length(hf)  
    s(i) = sum(dot(cross(hf_vector{i}, dual_w), fnormal, 2).*farea);
end

A = zeros(length(hf), length(hf));
for i = 1 : length(hf)
    for j = 1 : length(hf)
        wp = wedge_product(mesh, hf{i}, hf{j});
        %integrate this 2-form
        A(i,j) = sum(wp.*farea);
    end
end
mu = A\s;

cw = zeros(size(w));
for i = 1 : length(hf)
    cw = cw + mu(i) * hf{i};
end

end