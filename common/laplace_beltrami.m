%% laplace_beltrami 
% Laplace Beltrami operator on the mesh.
% 
% Cotangent formula is used, while there are some variants:
% 
% * 'Polthier', see paper [1]
% * 'Meyer', see paper [2]
% * 'Desbrun', see paper [3]
% 
% For comparison and convergence analysis, see paper [4]
% 
% # K. Polthier. Computational Aspects of Discrete Minimal Surfaces. In 
%   Proc. of the Clay Summer School on Global Theory of Minimal Surfaces, 
%   J. Hass, D. Hoffman, A. Jaffe, H. Rosenberg, R. Schoen, M. Wolf (Eds.), 
%   to appear, 2002.
% # M. Meyer, M. Desbrun, P. Schröder, and A. Barr. Discrete 
%   Differential-Geometry Operator for Triangulated 2-manifolds. In Proc. 
%   VisMath'02, Berlin, Germany, 2002.
% # M. Desbrun, M. Meyer, P. Schröder, and A. H. Barr. Implicit Fairing 
%   of Irregular Meshes using Diffusion and Curvature Flow. SIGGRAPH99, 
%   pages 317-324, 1999.
% # Xu, Guoliang. "Convergent discrete laplace-beltrami operators over 
%   triangular surfaces." Geometric Modeling and Processing, 2004. 
%   Proceedings. IEEE, 2004.
%
%% Syntax
%   A = laplace_beltrami(mesh)
%   A = laplace_beltrami(mesh, method)
%
%% Description
%  face  : double array, nf x 3, connectivity of mesh
%  vertex: double array, nv x 3, vertex of mesh
%  method: string, optional, method of cotangent formula, can be one of
%          three: 'Polthier', 'Meyer', 'Desbrun'. Default is 'Polthier'.
% 
%  A: sparse matrix, nv x nv, Laplace Beltrami operator
%
%% Example
%   A = laplace_beltrami(mesh)
%   A = laplace_beltrami(mesh,'Polthier') % same as last 
%   A = laplace_beltrami(mesh,'Meyer')
%   A = laplace_beltrami(mesh,'Desbrun')
%
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/03
%  Revised: 2014/03/03 by Wen, add more cotangent formula variants, not
%           implemented
%  Revised: 2014/03/23 by Wen, add doc
%  Revised: 2014/03/28 by Wen, add code for 'Meyer' and 'Desbrun' methods
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function A = laplace_beltrami(mesh, method)
% default method is 'Polthier'
if nargin == 1
    method = 'Polthier';
end
face = mesh.face;
edge = mesh.edge;
vert = mesh.vert;
if ~isfield(mesh,'ew')
    ew = edge_weight(mesh);
else
    ew = mesh.ew;
end

switch method
    case 'Polthier'
        A = sparse([edge(:,1);edge(:,2)],[edge(:,2);edge(:,1)],[ew;ew]);
        sA = sum(A,2);
        A = A - diag(sA);
    case 'Meyer'
        va = vertex_area(face,vert,'mixed');
        ew = (ew./va(edge(:,1))+ew./va(edge(:,2)))/2;
        A = sparse([edge(:,1);edge(:,2)],[edge(:,2);edge(:,1)],[ew;ew]);
        sA = sum(A,2);
        A = A - diag(sA);
    case 'Desbrun'
        va = vertex_area(face,vert,'one_ring');
        ew = (ew./va(edge(:,1))+ew./va(edge(:,2)))/2*3;
        A = sparse([edge(:,1);edge(:,2)],[edge(:,2);edge(:,1)],[ew;ew]);
        sA = sum(A,2);
        A = A - diag(sA);
    otherwise
        error('Wrong method. Available methods are: Polthier, Meyer, Desbrun')
end
