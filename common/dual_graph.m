%% compute dual graph 
% Dual graph of a triangle mesh, regarded as graph. 
% Each face in original mesh corresponds to a vertex in dual graph, vertex
% position be the centroid of the original face.
%
%% Syntax
%   [amf] = compute_dual_graph(face);
%   [amf,dual_vertex] = compute_dual_graph(face,vertex);
%
%% Description
%  face  : double array, nf x 3, connectivity of mesh
%  vertex: double array, nv x 3, vertex of mesh
% 
%  amf: sparse matrix, nf x nf, connectivity of dual graph
%  dual_vertex: nf x 3, dual vertex in dual graph, if vertex is not
%               supplied, will return []
%
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/14
%  Revised: 2014/03/18 by Wen, add doc
%  Revised: 2014/03/23 by Wen, revise doc
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function G = dual_graph(mesh)
eif = mesh.eif;
ind1 = eif(:,1)>0;
eif(~ind1,1) = inf;
ind2 = eif(:,2)>0;
eif(~ind2,2) = inf;
ind = ind1&ind2;
G = graph(eif(ind,1),eif(ind,2));
