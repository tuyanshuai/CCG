function E = harmonic_energy(mesh, f)
% compute harmonic energy for a map f
edge = mesh.edge;
% if edge weight is computed already, reuse it
if ~isfield(mesh,'ew')
    ew = edge_weight(mesh);
else
    ew = mesh.ew;
end
% harmonic energy
df = f(edge(:,2),:)-f(edge(:,1),:);
E = sum(ew.*dot(df,df,2))/2;
