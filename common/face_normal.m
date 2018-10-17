function normal = face_normal(face, vert)
normal = cross(vert(face(:, 2), :) - vert(face(:, 1), :), ...
    vert(face(:, 3), :) - vert(face(:, 1), :));
normal = bsxfun(@rdivide,normal,sqrt(sum(normal.^2,2)));
end