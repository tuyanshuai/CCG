function path = shortestpath(mesh,cc)
nv = mesh.nv;
vert = mesh.vert;
edge = mesh.edge;
de = vert(edge(:,1),:)-vert(edge(:,2),:);
el = sqrt(dot(de,de,2));
G = sparse(edge(:,1),edge(:,2),el,nv,nv);
G = G+G';
path = [];
for i = 2:length(cc)
    [~,pathi] = graphshortestpath(G,cc(i-1),cc(i));
    path = [path,pathi(1:end-1)];
end
path = [path,cc(end)];
