%% EXAMPLE 1D-BEAM using NLvib
% Example based on:
% Sombroek et al. (2018). Numerical computation of nonlinear normal modes 
% in a modal derivative subspace. Computers and Structures, 195, 34–46. 
% https://doi.org/10.1016/j.compstruc.2017.08.016
%
% Author: Jacopo Marconi, Politecnico di Milano
% Created: 21 April 2021
% Last modified: 27 April 2021

clear
clc
close all

imod = 1; % mode analyze
% NOTE: you can load "BeamNLvib.mat" which contains the results for the
% beam meshed with 8 elements (all other parameters set as in this script).
% Run the PLOT sections to inspect the results.

%% Parameters                                                       
% geometry
len = 1;        	% length
height = 1e-2;    	% height in the bending direction
thickness = 1e-2;	% thickness in the third dimension

% mesh
nElements = 30;
dx = len / nElements;

% material properties
E       = 210e9;    % Young's modulus
rho     = 8750;     % density
nu      = 0.3;      % nu

%% Structural Model                                                 
% Material
myBeamMaterial  = KirchoffMaterial();
set(myBeamMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);

% Element
myElementConstructor = @()BeamElement(thickness, height, myBeamMaterial); % same element all across the domain

% Mesh
x = (0:dx:len).';
nNodes = size(x,1);
Nodes = [x, zeros(nNodes,1)];

Elements = [1:nNodes-1;2:nNodes].';
BeamMesh = Mesh(Nodes);
BeamMesh.create_elements_table(Elements,myElementConstructor);

nNodes = BeamMesh.nNodes;
BeamMesh.set_essential_boundary_condition([1 nNodes],1:3,0);

ndofs = length( BeamMesh.EBC.unconstrainedDOFs );
ntot = BeamMesh.nDOFs;

% ASSEMBLY ________________________________________________________________
BeamAssembly = Assembly(BeamMesh);
M = BeamAssembly.mass_matrix();
u0 = zeros( BeamMesh.nDOFs, 1);
[K,~] = BeamAssembly.tangent_stiffness_and_force(u0);

% store matrices
BeamAssembly.DATA.K = K;
BeamAssembly.DATA.M = M;

% External force __________________________________________________________
forced_dof = round(nNodes/2) * 3 - 1; % middle node, vertical dof
Fext = zeros(ntot, 1);
Fext( forced_dof ) = 1;
Fextc = BeamAssembly.constrain_vector( Fext );

% Let us also define the index of the forced dof in the constrained vector:
forced_dof_c = BeamAssembly.free2constrained_index( forced_dof );

%% Eigenvalue problem                                               
n_VMs = ndofs; % first n_VMs modes with lowest frequency calculated 
Kc = BeamAssembly.constrain_matrix(K);
Mc = BeamAssembly.constrain_matrix(M);
[Phi,om2] = eigs(Kc, Mc, n_VMs, 'SM'); %Phi-->V0 in fe code
[om, ind] = sort(sqrt(diag(om2)));
f0 = om/2/pi;
Phi = Phi(:,ind);
for ii = 1:n_VMs
    Phi(:,ii) = Phi(:,ii)/max(sqrt(Phi(1:3:end,ii).^2+Phi(2:3:end,ii).^2));
end
Phi = BeamAssembly.unconstrain_vector(Phi);

u = Phi(1:3:end, imod);
v = Phi(2:3:end, imod);
x = Nodes(:, 1);
y = Nodes(:, 2);

% Damping _________________________________________________________________
Qfactor = 100;
csi = 1./(2*Qfactor);   % dimensionless damping
om0 = 2*pi*f0(1);       % first eigenfrequency
alfa = 2 * om0 * csi;
D = alfa*M;             % Mass-proportinal damping: D = a*M
BeamAssembly.DATA.D = D;
Dc = BeamAssembly.constrain_matrix(D);

%% reduced basis
m = 3; % use the first m VMs in reduction
V = Phi(:,1:m);
Vc = BeamAssembly.constrain_vector(Phi(:,1:m));
Kc = BeamAssembly.constrain_matrix(K);
Mc = BeamAssembly.constrain_matrix(M);
figure
for i=1:m
    LHS_matrix=[full(Kc-om(i)^2*Mc) -Mc*Vc(:,i);
                (-Mc*Vc(:,i))' 0];
    for j=1:m
        
        mystiffness_derivative=full(BeamAssembly.stiffness_derivative(u0,V(:,j)));
        mystiffness_derivative=BeamAssembly.constrain_matrix(mystiffness_derivative);
        RHS_vector=[-mystiffness_derivative*Vc(:,i);0];
        myresult=LHS_matrix^-1*RHS_vector;
        modal_der=myresult(1:end-1);
        modal_der = BeamAssembly.unconstrain_vector(modal_der);
        subplot(m,(m+1)/2+1,(i-1)*m+j)
        plot(modal_der(1:3:end))
        title(append('dphi_',num2str(i),'/dq_',num2str(j)))
        if i<=j
            V=[V modal_der];
        end
%     figure
%     plot(modal_der)
    end
end
sgtitle('modal derivatives in plane') 
%% Reduced assembly 
BeamReducedAssembly = ReducedAssembly(BeamMesh,V);

BeamReducedAssembly.DATA.M = BeamReducedAssembly.mass_matrix();
BeamReducedAssembly.DATA.D = V'*D*V;
BeamReducedAssembly.DATA.K =  BeamReducedAssembly.stiffness_matrix();


ndofs = length(BeamReducedAssembly.DATA.M);
Fext_red=V'*Fext ;
% Fextc_red = BeamReducedAssembly.constrain_vector( Fext_red );
%how to constrain reduced vectors (is it even possible?)
%% (2) NLvib:  FRF - Harmonic Balance reduced                            	
BeamSystem_red = FE_system( BeamReducedAssembly, Fext_red );
% BeamSystem_red.V=V;

PHI_lin = BeamReducedAssembly.constrain_vector(Phi);

omi = om(imod);         % linear eigenfrequency

% Analysis parameters
H = 7;                  % harmonic order
N = 3*H+1;              % number of time samples per period
Om_s = omi * 0.95;   	% start frequency
Om_e = omi * 1.1;    	% end frequency
ds = 1;                 % Path continuation step size
exc_lev = [1];       

fprintf('\n\n FRF from %.2f to %.2f rad/s \n\n', Om_s, Om_e)

% Excitation levels
r2 = cell(length(exc_lev),1);
for iex = 1 : length(exc_lev)
    % Set excitation level
    BeamSystem_red.Fex1 = Fext_red * exc_lev(iex);
    
    % Initial guess (solution of underlying linear system)
    Q1 = (-Om_s^2*BeamSystem_red.M + 1i*Om_s*BeamSystem_red.D +...
        BeamSystem_red.K) \ Fext_red;
    y0 = zeros( (2*H+1)*ndofs , 1);
    y0( ndofs + (1:2*ndofs) ) = [real(Q1);-imag(Q1)];
    
    % stuff for scaling
    qscl = max(abs((-omi^2*BeamSystem_red.M + ...
        1i*omi*BeamSystem_red.D + BeamSystem_red.K) \ Fext_red));
    dscale = [y0*0+qscl; omi];
    Sopt = struct('Dscale', dscale, 'dynamicDscale', 1);
    
    % Solve and continue w.r.t. Om	
    [X2, Solinfo, Sol] = solve_and_continue(y0, ...
        @(X) HB_residual_reduced(X, BeamSystem_red, H, N, 'FRF',V), ...
        Om_s, Om_e, ds, Sopt);
    
    % Interpret solver output
    r2{iex} = nlvib_decode(X2, Solinfo, Sol, 'FRF', 'HB', ndofs, H);
    
    results.FRF.HB{iex} = r2{iex};
end
%% (2) PLOT                                                         

r2 = results.FRF.HB;

figure;hold on
h = 1;
for iex = 1 : length(exc_lev)
    % 1st harmonic amplitude of the forced dof (use force_dof_c!)
    Qre_full=pagemtimes( V,r2{iex}.Qre);
    Qim_full=pagemtimes(V,r2{iex}.Qim);
    A = Qre_full(forced_dof_c, :, h);
    B = Qim_full(forced_dof_c, :, h);
    W = r2{iex}.omega;
    a1 = squeeze(sqrt( A.^2 + B.^2 ));
    plot(W, a1 / height, 'linewidth', 2,"Color","r","DisplayName",...
        "nlvib nonlin reduced"); hold on
end
grid on
axis tight
xlabel('\omega [rad/s]')
ylabel('|Q_1| / height [-]')
title('FRF with Harmonic Balance')


load("nlvib_freq_resp_30_8750.mat")
plot(nlvib_lin_freq,nlvib_lin_amp,"DisplayName","nlvib lin",...
    "LineStyle","--","LineWidth",2,"Color","k")
plot(nlvib_nonlin_freq,nlvib_nonlin_amp,"DisplayName","nlvib nonlin full",...
    "LineStyle",":","LineWidth",2,"Color","k")
legend