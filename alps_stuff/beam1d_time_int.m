clear
clc
close all
%% parameters

% geometry
len = 1;        	% length
height = 1e-2;    	% height in the bending direction
thickness = 1e-2;	% thickness in the third dimension

% mesh
nElements = 30;
dx = len/nElements;

% material properties
E       = 210e9;  % Young's modulus
rho     = 7800;   % density
nu      = 0.3;    % nu

%% Structural
% Material
myBeamMaterial  = KirchoffMaterial();
set(myBeamMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
% D               = myBeamMaterial.get_stress_strain_matrix_2D();

% Element
myElementConstructor = @()BeamElement(thickness, height, myBeamMaterial); % same element all across the domain

% Mesh
x = (0:dx:len).';
nNodes = size(x,1);
Nodes = [x, zeros(nNodes,1)];

Elements = [1:nNodes-1;2:nNodes].';
BeamMesh = Mesh(Nodes);
BeamMesh.create_elements_table(Elements,myElementConstructor);

BeamMesh.set_essential_boundary_condition([1 nElements+1],1:3,0);

% ASSEMBLY ________________________________________________________________
BeamAssembly = Assembly(BeamMesh);
M = BeamAssembly.mass_matrix();
u0 = zeros( BeamMesh.nDOFs, 1);
[K,~] = BeamAssembly.tangent_stiffness_and_force(u0);
C = BeamAssembly.damping_matrix();

% store matrices
BeamAssembly.DATA.K = K;
BeamAssembly.DATA.M = M;
BeamAssembly.DATA.C = C;


%% Dynamic response using Implicit Newmark
n_vm=1;
Kc = BeamAssembly.constrain_matrix(K);
Mc = BeamAssembly.constrain_matrix(M);
evals=sqrt(eigs(Kc/Mc,n_vm,"smallestabs"));
omega_ext=evals(n_vm);
% omega_ext=500;
% omega_ext=921;
% omega_ext=100;

T =  2*pi/omega_ext; % time period of forcing

% % load amplification factor
% amplification_factor = 1;
% 
% % forcing function
% Pressure = 1; % in Pascals 
% F = Pressure*BeamAssembly.uniform_body_force();
% F_ext = @(t) amplification_factor * F * sin(omega_ext * t);

% External force __________________________________________________________
forced_dof = round(nNodes/2) * 3 - 1; % middle node, vertical dof
F = zeros(BeamMesh.nDOFs, 1);
F( forced_dof ) = 1;
F_ext = @(t) F * sin(omega_ext * t);

TI_lin.Solution.u=zeros(3*length(x),1);

% Initial condition: equilibrium
u0 = TI_lin.Solution.u;
v0 = TI_lin.Solution.u;
a0 = TI_lin.Solution.u; % a0 = M\(F_ext(0)-C*v0-F(u0)) 

q0 = BeamAssembly.constrain_vector(u0);
qd0 = BeamAssembly.constrain_vector(v0);
qdd0 = BeamAssembly.constrain_vector(a0);

% time step for integration
dt = T/50;

% Instantiate object for linear time integration
TI_lin = ImplicitNewmark('timestep',dt,'alpha',0.005,'linear',true);

% Linear Residual evaluation function handle
residual_lin = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,BeamAssembly,F_ext);

% Linearized Time Integration
tmax = 40*T; 
TI_lin.Integrate(q0,qd0,qdd0,tmax,residual_lin);

% obtain full solution
TI_lin.Solution.u = BeamAssembly.unconstrain_vector(TI_lin.Solution.q);

% Animate solution on Mesh (very slow)
% AnimateFieldonDeformedMesh(BeamMesh.nodes,BeamMesh.Elements,TI_lin.Solution.u ,'factor',1,'index',1:3,'filename','lineardisp')
% AnimateFieldonDeformedMesh(Nodes,Elements,TI_lin.Solution.u ,'factor',10,'index',1:3,'filename','lineardisp')

%% nonlin integration 
% Instantiate object for nonlinear time integration
TI_NL = ImplicitNewmark('timestep',dt,'alpha',0.005);

% Linear Residual evaluation function handle
residual = @(q,qd,qdd,t)residual_nonlinear(q,qd,qdd,t,BeamAssembly,F_ext);

% Nonlinear Time Integration
% tmax = 50*T; 
TI_NL.Integrate(q0,qd0,qdd0,tmax,residual);
TI_NL.Solution.u = BeamAssembly.unconstrain_vector(TI_NL.Solution.q);

%% plotting


figure;hold on
plot(TI_lin.Solution.time, TI_lin.Solution.u(forced_dof,:)/height,"DisplayName","linear")
plot(TI_NL.Solution.time, TI_NL.Solution.u(forced_dof,:)/height,"DisplayName","nonlinear")
xlabel('time'); 
% ylabel('q'); 
ylabel('|Q_1| / height [-]')
grid on;
legend

%% FFT on last k cycles
k=2;
% last_k_cycles_time=TI_NL.Solution.time(end-k*T/dt:end);
last_k_cycles_sol=TI_NL.Solution.u(forced_dof,end-k*T/dt:end);
% L=length(last_k_cycles_time)-1;
% 
Fs=1/dt;
% Y = fft(last_k_cycles_sol);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = 2*pi*Fs*(0:(L/2))/L;
[freq, spectrummag, ~] = fft_yFs(last_k_cycles_sol, Fs);
figure; hold on
plot(freq,spectrummag,"DisplayName","nonlin")
xlabel('f (Hz)')
ylabel('|P1(f)|')

% last_k_cycles_time=TI_lin.Solution.time(end-k*T/dt:end);
last_k_cycles_sol=TI_lin.Solution.u(forced_dof,end-k*T/dt:end);
% L=length(last_k_cycles_time)-1;

Fs=1/dt;
% Y = fft(last_k_cycles_sol);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = 2*pi*Fs*(0:(L/2))/L;

[freq, spectrummag, ~] = fft_yFs(last_k_cycles_sol, Fs);
plot(freq,spectrummag,"DisplayName","lin")
legend