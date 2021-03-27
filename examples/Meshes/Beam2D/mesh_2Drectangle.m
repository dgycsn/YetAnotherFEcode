% mesh_2Drectangle
%
% syntax:
% [nodes, elements, nset] = ...
%           mesh_2Drectangle(Lx,Ly,nx,ny)
%
% Description: create the mesh for a rectangle using quadratic
%              quadrilaterals (Quad8).
%
% INPUTS    Lx, Ly: length, width
%           nx,ny: number of elements for Lx and Ly, respectively.
%           
% OUTPUTS   nodes --> [x, y]
%           elements --> [node_1, ... , node_8]
%           nset --> struct containing node sets (IDs) of the 4 external
%           edges. They are ordered as: x=0, y=0, x=l, y=w.

function [nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny)

nel = nx*ny;

fprintf(' Meshing %d elements (Q8) ... ', nel)
tic

elnodes = [         % element coordinates in the natural space (/2)
     0     0
     1     0
     1     1
     0     1
     0.5   0
     1     0.5
     0.5   1
     0     0.5];
lx = Lx/nx;
ly = Ly/ny;
elnodes(:,1) = elnodes(:,1)*lx; % element coordinates in physical space
elnodes(:,2) = elnodes(:,2)*ly;
% create nodes
nodes = zeros(nel*8,2);
nn = 1;
for ii = 1:nx
    for jj = 1:ny   
        elnodes_temp = elnodes;
        elnodes_temp(:,1) = elnodes_temp(:,1)+lx*(ii-1);
        elnodes_temp(:,2) = elnodes_temp(:,2)+ly*(jj-1);
        nodes(nn:nn+7,:)  = elnodes_temp;
        nn = nn+8;
    end
end

tol = 1e12;
nodes = round(nodes*tol) / tol;

% remove duplicate nodes from 'nodes' matrix
[nodes, ~, ic] = unique(nodes, 'rows', 'stable');
elements = ic;
elements = reshape(elements, 8, nel)';
idb = reshape(1 : length(nodes(:)), 2, size(nodes,1))';

conn = zeros(nel,8*2);
for ii = 1:nel
    conn(ii,:) = reshape(idb(elements(ii,:),:)',1,8*2);
end

% external surfaces
tol = 1e-15;
node_IDs = 1 : size(nodes, 1);
edge1_X0_nodes = node_IDs( nodes(:,1)==0 );
edge2_Y0_nodes = node_IDs( nodes(:,2)==0 );
edge3_XL_nodes = node_IDs( abs(nodes(:,1)-Lx)<tol );
edge4_YW_nodes = node_IDs( abs(nodes(:,2)-Ly)<tol );

% edge1_X0_dofs=reshape(idb(edge1_X0_nodes,2:end)',1,length(edge1_X0_nodes)*2);
% edge2_Y0_dofs=reshape(idb(edge2_Y0_nodes,2:end)',1,length(edge2_Y0_nodes)*2);
% edge3_XL_dofs=reshape(idb(edge3_XL_nodes,2:end)',1,length(edge3_XL_nodes)*2);
% edge4_YW_dofs=reshape(idb(edge4_YW_nodes,2:end)',1,length(edge4_YW_nodes)*2);

nset = {edge1_X0_nodes, edge2_Y0_nodes, edge3_XL_nodes, edge4_YW_nodes};

fprintf(' %.2f s\n',toc)
fprintf('  Nodes: %d \n',size(nodes,1))
fprintf('  Dofs:  %d \n\n',size(nodes,1)*2)

