function f = integration(mesh, dw)
% integration of differential one form
edge = mesh.edge;
nv = mesh.nv;
es = sparse(edge(:,1),edge(:,2),dw,nv,nv);
es = es - conj(es');
vvr = vert_vert_ring(mesh);
f = zeros(mesh.nv,1);
f(1) = 0;
ind = false(nv,1);
ind(1) = true;
qe = 1;
while ~isempty(qe)
    i = qe(end);
    qe(end) = [];
    vri = vvr{i};
    for j = vri
        if ind(j)
            continue
        end
        f(j) = f(i) + es(i,j);
        ind(j) = true;
        qe = [qe,j];
    end
end
