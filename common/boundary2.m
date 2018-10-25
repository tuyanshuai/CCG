function bds = boundary2(mesh)
nv = mesh.nv;
ind = false(nv,1);
ind(mesh.bd) = true;
% find boundary halfedge
indhe = ind(mesh.halfedge(:,1))&ind(mesh.halfedge(:,2));
he = mesh.halfedge(indhe,:);
G = sparse(he(:,1),he(:,2),ones(size(he,1),1),nv,nv);
bds = {};
ind = false(nv,1);
ind(he) = true;
k = 1;
while true
    % find a boundary vert
    b = find(ind,1,'first');
    if isempty(b)
        break;
    end
    bd = b;
    ind(b) = false;
    while true
        bs = find(G(b,:));
        if isempty(bs)
            error(['unable to find halfedge with source vert ' num2str(b)]);
        end
        
        if length(bs)>1
            ib = ind(bs);
            i = find(ib,1);
            if isempty(i)
                break
            end
            bs = bs(i);
        end
        if ~ind(bs)
            break
        end

        ind(bs) = false;
        b = bs;
        bd = [bd;b];
    end
    bds{k} = bd;
    k = k+1;
end
