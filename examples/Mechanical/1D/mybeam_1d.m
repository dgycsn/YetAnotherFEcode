clear
clc
close
%% parameters
% geometry
l = 1; 
h = 1e-2;
b = 1e-2; 

% mesh
nElements = 30;
dx = l/nElements;

% material properties
E       = 210e9;  % Young's modulus
rho     = 8750;   % density
nu      = 0.3;    % nu

%% Structural
% Material
myBeamMaterial  = KirchoffMaterial();
set(myBeamMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
D               = myBeamMaterial.get_stress_strain_matrix_2D();

% Element
myElementConstructor = @()BeamElement(b, h, myBeamMaterial); % same element all across the domain

% Mesh
x = (0:dx:l).';
nNodes = size(x,1);
Nodes = [x, zeros(nNodes,1)];

Elements = [1:nNodes-1;2:nNodes].';
BeamMesh = Mesh(Nodes);
BeamMesh.create_elements_table(Elements,myElementConstructor);

BeamMesh.set_essential_boundary_condition([1 31],1:2,0);

% ASSEMBLY ________________________________________________________________
BeamAssembly = Assembly(BeamMesh);
M = BeamAssembly.mass_matrix();
u0 = zeros( BeamMesh.nDOFs, 1);
[K,~] = BeamAssembly.tangent_stiffness_and_force(u0);

% store matrices
BeamAssembly.DATA.K = K;
BeamAssembly.DATA.M = M;

% Eigenvalue problem_______________________________________________________
n_VMs = 2; % first n_VMs modes with lowest frequency calculated 
Kc = BeamAssembly.constrain_matrix(K);
Mc = BeamAssembly.constrain_matrix(M);
[V0,om] = eigs(Kc, Mc, n_VMs, 'SM');
[f0,ind] = sort(sqrt(diag(om))/2/pi);
V0 = V0(:,ind);
for ii = 1:n_VMs
    V0(:,ii) = V0(:,ii)/max(sqrt(sum(V0(:,ii).^2,2)));
end
V0 = BeamAssembly.unconstrain_vector(V0);

% PLOT
mod = 1;
elementPlot = Elements;
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(Nodes, elementPlot, 0);
v1 = reshape(V0(:,mod), 2, []).';
PlotFieldonDeformedMesh(Nodes, elementPlot, v1, 'factor', Ly*1.1);
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ' Hz'])
