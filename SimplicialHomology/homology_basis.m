%% homology basis 
% Compute a basis for the homology group H_1(M,Z), based on the algorithm
% 6 in book [1].
%  
% # Gu, Xianfeng David, and Shing-Tung Yau, eds. Computational conformal
%   geometry. Vol. 3. Somerville: International Press, 2008.
%
%% Syntax
%   hb = homology_basis(mesh)
%
%% Description
%  mesh: mesh structure
%
%  hb: cell array, n x 1, a basis of homology group, each cell is a closed 
%      loop based. Return empty for genus zero surface. Two loops on each
%      handle. If there is boundary on surface, each boundary will be an
%      element of hb
%
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/13
%  Revised: 2014/03/24 by Wen, add doc
%  Revised: 2018/09/16 by Wen, simplify code
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui
function hb = homology_basis(mesh)
face = mesh.face;
ee = cut_graph(face);
G = graph(ee(:,1),ee(:,2));
[tree,pred] = minspantree(G,'Type','forest','Root',ee(1));
v = ee(1);
tree_edges = table2array(tree.Edges);
I = [tree_edges(:,1);tree_edges(:,2)];
J = [tree_edges(:,2);tree_edges(:,1)];
eh = setdiff(ee,[I,J],'rows');
hb = cell(size(eh,1),1);
for i = 1:size(eh,1)
    p1 = trace_path(pred,eh(i,1),v);
    p2 = trace_path(pred,eh(i,2),v);
    loop = [flipud(p1);eh(i,1);eh(i,2);p2];
    hb{i} = prune_path(loop);
end

function path = trace_path(pred,v,root)
path = [];
while true
    path = [path;pred(v)];
    v = pred(v);
    if v == root
        break;
    end
end

function path_new = prune_path(path)
ind = path ~= flipud(path);
i = find(ind,1)-1;
path_new = path(i:end-i+1);
