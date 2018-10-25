%% cut graph 
% Compute a cut graph of mesh, such that surface becomes simply-connected
% if slice mesh along the cut-graph. There are two versions: if both face 
% and vertex provided, invoke version 1; if only face provided, invoke 
% version 2. Version 1 is exact implementation of algorithm 3 in book [1].
% Version 2 is translated from David Gu's C++ code of cut graph, which is
% much faster than version 1.
% 
% Though version 1 takes vertex into consideration, both algorithms do not
% generated optimal cut graph (shortest one). In fact this problem seems
% to be open until now.
%
%% Syntax
%   ee = cut_graph(face)
%
%% Description
%  face  : double array, nf x 3, connectivity of mesh
% 
%  ee: double array, n x 2, edges in the cut graph, each row is an edge on 
%      mesh, may not be in consecutive order. ee's mainly purpose is been 
%      passed to slice_mesh, which will slice the mesh open along edges in
%      ee, to form a simply-connected surface
%
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/13
%  Revised: 2014/03/13 by Wen, implement another cut graph algorithm
%  Revised: 2014/03/17 by Wen, merge two cut graph algorithm
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function ee = cut_graph(face)
[am,amd] = adjacency_matrix(face);
nf = size(face,1);
% use array to emulate queue
queue = zeros(nf,1);
queue(1) = 1;
qs = 1; % point to queue start
qe = 2; % point to queue end

ft = false(nf,1);
ft(1) = true;
face4 = face(:,[1 2 3 1]);

while qe > qs
    fi = queue(qs);
    qs = qs+1;
    for i = 1:3
        he = face4(fi,[i i+1]);
        sf = amd(he(2),he(1));
        if sf       
            if ~ft(sf)
                queue(qe) = sf;
                qe = qe+1;
                ft(sf) = true;
                am(he(1),he(2)) = -1;
            end
        end
    end
end
am((am<0)') = 0;
G = triu(am>0);

% prune the graph cut
while true
    Gs = full(sum(G,2))+full(sum(G,1))';
    ind = (Gs == 1);
    if sum(ind) ==0
        break;
    end
    G(ind,:) = 0;
    G(:,ind) = 0;
end

[I,J,~] = find(G);
ee = [I,J];
