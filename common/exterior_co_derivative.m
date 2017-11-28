function delta_w = exterior_co_derivative(mesh, w, order)
% exterior co-derivative of differential form w
edge = mesh.edge;
nv = mesh.nv;
if ~isfield(mesh,'ew')
    ew = edge_weight(mesh);
else
    ew = mesh.ew;
end
if order == 1
    if size(w) ~= mesh.ne
        error('differential form w has incorrect shape');
    end
    delta_w = accumarray(edge(:,1),ew.*w,[nv,1]) + accumarray(edge(:,2),-ew.*w,[nv,1]);
end
if order == 2
    if size(w) ~= mesh.nf
        error('differential form w has incorrect shape');
    end
end
    