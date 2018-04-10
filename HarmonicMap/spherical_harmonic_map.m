function uvw = spherical_harmonic_map(face,vert)
nv = size(vert,1);
A = laplace_beltrami(face,vert);
[V,D] = eigs(A,2,0);
fd = V(:,2);
ind = fd>0;
indf = sum(ind(face),2)>=2;
f1 = face(indf,:);
ind = false(nv,1);
ind(f1) = true;
uv1 = disk_harmonic_map(f1,vert);
bd1 = compute_bd(f1);
uv = disk_harmonic_map(face,vert,bd1,uv1(bd1,:));
uvw = zeros(nv,3);
uvw(:,1:2) = uv*2;
uvw(:,3) = dot(uv,uv,2)-1;
uvw = uvw./(dot(uv,uv,2)+1);
uvw(~ind,3) = -uvw(~ind,3);

vr = compute_vertex_ring(face,uvw,bd1);
for i = 1:length(vr)
    ri = vr{i};
    x = mean(uvw(ri,:));
    uvw(bd1(i),:) = x./norm(x);
end

dt = 0.2;
k = 0;
while k <= 3000
    d = A*uvw;
    dd = d-dot(d,uvw,2).*uvw;
    if mod(k,100) == 0
        fprintf('#%d |df/dt| = %.10f\n',k,norm(dd));
    end
    if norm(dd) < eps
        break;
    end
    while true
        uvw1 = uvw + dt*(dd);
        l = sqrt(dot(uvw1,uvw1,2));
        uvw1 = uvw1./l;
        d = A*uvw1;
        dd1 = d-dot(d,uvw1,2).*uvw1;
        if norm(dd1)>norm(dd)
            dt = dt/2;
        else
            break;
        end
    end
    uvw = uvw + dt*(dd);
    l = sqrt(dot(uvw,uvw,2));
    uvw = uvw./l;
    k = k+1;
end
