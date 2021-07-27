C% clear
% clc
% close all
%% parameters

% geometry
l = 1; 
h = 1e-2;
b = 1e-2; 

% mesh
nElements = 4;
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

freq=250:10:400;
lin_resp=zeros(length(freq),1);
nonlin_resp=zeros(length(freq),1);

%first initial condition
init_cond=zeros(BeamAssembly.Mesh.nDOFs, 1);

TI_lin.Solution.u=init_cond;
TI_lin.Solution.qd=BeamAssembly.constrain_vector(init_cond);

TI_NL.Solution.u=init_cond;
TI_NL.Solution.qd=BeamAssembly.constrain_vector(init_cond);

for i=1:length(freq)
%% Dynamic response using Implicit Newmark
% n_vm=3;
% Kc = BeamAssembly.constrain_matrix(K);
% Mc = BeamAssembly.constrain_matrix(M);
% evals=sqrt(eigs(Kc/Mc,n_vm,"smallestabs"));
% omega_ext=evals(n_vm);
% omega_ext=921;
% omega_ext=100;

omega_ext=freq(i)
pause(1)

T =  2*pi/omega_ext; % time period of forcing

% load amplification factor
amplification_factor = 1e-4;

% forcing function
Pressure = 1e6; % in Pascals 
F = Pressure*BeamAssembly.uniform_body_force();
F_ext = @(t) amplification_factor * F * sin(omega_ext * t);

% Initial condition: previous step
u0 = TI_lin.Solution.u(:,end);
% v0 = zeros(BeamAssembly.Mesh.nDOFs, 1);
a0 = zeros(BeamAssembly.Mesh.nDOFs, 1); % a0 = M\(F_ext(0)-C*v0-F(u0)) 

q0 = BeamAssembly.constrain_vector(u0);
qd0 = TI_lin.Solution.qd(:,end);
qdd0 = BeamAssembly.constrain_vector(a0);

% time step for integration
dt = T/50;

% Instantiate object for linear time integration
TI_lin = ImplicitNewmark('timestep',dt,'alpha',0.005,'linear',true);

% Linear Residual evaluation function handle
residual_lin = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,BeamAssembly,F_ext);

% Linearized Time Integration
tmax = 10*T; 
TI_lin.Integrate(q0,qd0,qdd0,tmax,residual_lin);

% obtain full solution
TI_lin.Solution.u = BeamAssembly.unconstrain_vector(TI_lin.Solution.q);

% Animate solution on Mesh (very slow)
% AnimateFieldonDeformedMesh(BeamMesh.nodes,BeamMesh.Elements,TI_lin.Solution.u ,'factor',1,'index',1:3,'filename','lineardisp')
% AnimateFieldonDeformedMesh(Nodes,Elements,TI_lin.Solution.u ,'factor',10,'index',1:3,'filename','lineardisp')

%% nonlin integration 

% Initial condition: equilibrium
u0 = TI_NL.Solution.u(:,end);
% v0 = zeros(BeamAssembly.Mesh.nDOFs, 1);
a0 = zeros(BeamAssembly.Mesh.nDOFs, 1); % a0 = M\(F_ext(0)-C*v0-F(u0)) 

q0 = BeamAssembly.constrain_vector(u0);
qd0 = TI_NL.Solution.qd(:,end);
qdd0 = BeamAssembly.constrain_vector(a0);

% Instantiate object for nonlinear time integration
TI_NL = ImplicitNewmark('timestep',dt,'alpha',0.005);

% Linear Residual evaluation function handle
residual = @(q,qd,qdd,t)residual_nonlinear(q,qd,qdd,t,BeamAssembly,F_ext);

% Nonlinear Time Integration
% tmax = 50*T; 
TI_NL.Integrate(q0,qd0,qdd0,tmax,residual);
TI_NL.Solution.u = BeamAssembly.unconstrain_vector(TI_NL.Solution.q);


%% plotting

dof = round(nNodes/2) * 3 - 1; % random degree of freedom at which time response is compared
% 
% 
% figure;hold on
% plot(TI_lin.Solution.time, TI_lin.Solution.u(dof,:),"DisplayName","linear")
% plot(TI_NL.Solution.time, TI_NL.Solution.u(dof,:),"DisplayName","nonlinear")
% xlabel('time'); ylabel('q'); grid on;
% legend

%% linear regression on last k cycles

k=2;

last_k_cycles_sol=TI_NL.Solution.u(dof,end-k*T/dt:end);
Y = fft(last_k_cycles_sol);
w=max(Y);
nonlin_resp(i)=norm(w);

last_k_cycles_sol=TI_lin.Solution.u(dof,end-k*T/dt:end);
Y = fft(last_k_cycles_sol);
w=max(Y);
lin_resp(i)=norm(w);
end

%% plot
figure;hold on
% plot(fliplr(freq),lin_resp,"DisplayName","linear")
% plot(fliplr(freq),nonlin_resp,"DisplayName","nonlin")

plot(freq,lin_resp,"DisplayName","linear")
plot(freq,nonlin_resp,"DisplayName","nonlin")

legend


