function dw = exterior_derivative(mesh, w, order)
% exterior derivative of differential form w
edge = mesh.edge;
if order == 0
    if size(w,1) ~= mesh.nv
        error('differential form w has incorrect shape');
    end
    dw = w(edge(:,2)) - w(edge(:,1));
end
if order == 1
    if size(w,1) ~= mesh.ne
        error('differential form w has incorrect shape');
    end
    face = mesh.face;
    nv = mesh.nv;
    ws = sparse(edge(:,1),edge(:,2),w,nv,nv);
    ws = ws - ws';
    dw = zeros(mesh.nf,1);
    dw = dw + ws(face(:,1)+(face(:,2)-1)*nv);
    dw = dw + ws(face(:,2)+(face(:,3)-1)*nv);
    dw = dw + ws(face(:,3)+(face(:,1)-1)*nv);
end
