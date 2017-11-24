function bds = boundary(face)
nv = max(max(face));
vert = zeros(nv,3);
tri = triangulation(face,vert);
ff = freeBoundary(tri);
if isempty(ff)
    bds = {};
    return;
end
S = sparse(ff(:,1),ff(:,2),ones(size(ff,1),1));
[I,~] = find(S);
bds = cell(1);

c = I(1);
b = c;
k = 1;
while true
    n = find(S(c,:));
    if length(n)>1
        n;
    end
    if isempty(n) || n(1) == b(1)
        bds{k} = b;
        k = k+1;
        I = setdiff(I,b);
        if isempty(I)
            break;
        else
            c = I(1);
            b = c;
        end
    else
        c = n(1);
        b = [b;c];        
    end
end

if length(bds) == 1
    bds = bds{1};
end
