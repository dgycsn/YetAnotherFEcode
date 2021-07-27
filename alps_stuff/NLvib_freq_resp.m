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
tic
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
[Phi,om2] = eigs(Kc, Mc, n_VMs, 'SM');
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
% Qfactor = 100;
% csi = 1./(2*Qfactor);   % dimensionless damping
% om0 = 2*pi*f0(1);       % first eigenfrequency
% alfa = 2 * om0 * csi;
% D = alfa*M;             % Mass-proportinal damping: D = a*M
% BeamAssembly.DATA.D = D;
% Dc = BeamAssembly.constrain_matrix(D);
disp(' Rayleigh Damping')
Qfactors = [100 200]';
frequenz = [1 2]';
csi = 1./(2*Qfactors);
om0 = 2*pi*f0(frequenz);
AA = [1./(2*om0) om0/2];
if size(AA,1)==size(AA,2)
    XX  = AA\csi;
else
    XX = (AA'*AA)\(AA'*csi);
end
a = XX(1);
b = XX(2);
fprintf(' alpha = %.2i, beta = %.2i\n\n',a,b)
% 
D=a*M+b*K;
BeamAssembly.DATA.D=D;
Dc=BeamAssembly.constrain_matrix(D);

%% (1) NLvib: NMA - Harmonic Balance                                
% Example adapted from "09_beam_cubicSpring_NM" from NLvib

BeamSystem = FE_system( BeamAssembly, Fext );

PHI_lin = BeamAssembly.constrain_vector(Phi);

%% (2) NLvib:  FRF - Harmonic Balance                             	

omi = om(imod);         % linear eigenfrequency

% Analysis parameters
H = 7;                  % harmonic order
N = 3*H+1;              % number of time samples per period
Om_s = 390;   	% start frequency
Om_e = 410;    	% end frequency
ds = 1;                 % Path continuation step size
exc_lev = [5];       

fprintf('\n\n FRF from %.2f to %.2f rad/s \n\n', Om_s, Om_e)

% Excitation levels
r2 = cell(length(exc_lev),1);
for iex = 1 : length(exc_lev)
    % Set excitation level
    BeamSystem.Fex1 = Fextc * exc_lev(iex);
    
    % Initial guess (solution of underlying linear system)
    Q1 = (-Om_s^2*Mc + 1i*Om_s*Dc + Kc) \ Fextc;
    y0 = zeros( (2*H+1)*ndofs , 1);
    y0( ndofs + (1:2*ndofs) ) = [real(Q1);-imag(Q1)];
    
    % stuff for scaling
    qscl = max(abs((-omi^2*Mc + 1i*omi*Dc + Kc) \ Fextc));
    dscale = [y0*0+qscl; omi];
    Sopt = struct('Dscale', dscale, 'dynamicDscale', 1);
    
    % Solve and continue w.r.t. Om	
    [X2, Solinfo, Sol] = solve_and_continue(y0, ...
        @(X) HB_residual(X, BeamSystem, H, N, 'FRF'), ...
        Om_s, Om_e, ds, Sopt);
    
    % Interpret solver output
    r2{iex} = nlvib_decode(X2, Solinfo, Sol, 'FRF', 'HB', ndofs, H);
    
    results.FRF.HB{iex} = r2{iex};
end

%% (2) PLOT                                                         

r2 = results.FRF.HB;

figure
h = 1;
for iex = 1 : length(exc_lev)
    % 1st harmonic amplitude of the forced dof (use force_dof_c!)
    A = r2{iex}.Qre(forced_dof_c, :, h);
    B = r2{iex}.Qim(forced_dof_c, :, h);
    W = r2{iex}.omega;
    a1 = squeeze(sqrt( A.^2 + B.^2 ));
    plot(W, a1 / height, 'linewidth', 1); hold on
end
grid on
axis tight
xlabel('\omega [rad/s]')
ylabel('|Q_1| / height [-]')
title('FRF with Harmonic Balance')

try
    % LINEAR RESPONSE
    % compute the linear FRF for comparison
    nw = 201;
    w_linear = linspace(Om_s, Om_e, nw);
    for iex = 1 : length(exc_lev)
        fr_linear = zeros(nw, ndofs);
        for ii = 1:nw
            w = w_linear(ii);
            fr_linear(ii,:) = (-w^2*Mc + 1i*w*Dc + Kc) \ Fextc * exc_lev(iex);
        end
        plot(w_linear, abs(fr_linear(:, forced_dof_c))/height, 'k--')
    end
    drawnow
end
toc
nlvib_lin_freq=w_linear;
nlvib_lin_amp=abs(fr_linear(:, forced_dof_c))/height;
nlvib_nonlin_freq=W;
nlvib_nonlin_amp=a1/height;
% save("nlvib_freq_resp_30_8750","nlvib_lin_amp","nlvib_lin_freq","nlvib_nonlin_amp","nlvib_nonlin_freq")
%%
s=sprintf("freq_resp/NLvib_"+...
    "Q1_"+num2str(Qfactors(1))+...
    "Q2_"+num2str(Qfactors(2))+...
    "_famp"+num2str(exc_lev)+".fig");
saveas(gcf,s); %saves the figure generated with a timestamp
s=sprintf("freq_resp/NLvib"+...
    "Q1_"+num2str(Qfactors(1))+...
    "Q2_"+num2str(Qfactors(2))+...
    "_famp"+num2str(exc_lev)+".png");
saveas(gcf,s); %saves the figure generated with a timestamp