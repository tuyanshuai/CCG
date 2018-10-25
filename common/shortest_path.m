function path = shortest_path(mesh,cc)
vert = mesh.vert;
edge = mesh.edge;
de = vert(edge(:,1),:)-vert(edge(:,2),:);
el = sqrt(dot(de,de,2));
G = graph([edge(:,1);edge(:,2)],[edge(:,2);edge(:,1)],[el;el]);
path = [];
for i = 2:length(cc)
    pathi = shortestpath(G,cc(i-1),cc(i));
    path = [path,pathi(1:end-1)];
end
path = [path,cc(end)];
path = path(:);
