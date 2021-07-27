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
ntot = BeamMesh.nDOFs;

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
% n_vm=1;
% Kc = BeamAssembly.constrain_matrix(K);
% Mc = BeamAssembly.constrain_matrix(M);
% evals=sqrt(eigs(Kc/Mc,n_vm,"smallestabs"));
% omega_ext=evals(n_vm);
% omega_ext=500;
% omega_ext=921;
% omega_ext=100;
% Eigenvalue problem
n_VMs = 5; % first n_VMs modes with lowest frequency calculated 
K_red = BeamAssembly.constrain_matrix(K);
M_red = BeamAssembly.constrain_matrix(M);
[V0,omega2] = eigs(K_red,M_red,n_VMs,'SM');
omega = sqrt(diag(omega2));

V0 = BeamAssembly.unconstrain_vector(V0);

omega_ext=omega(1);
% omega_ext=325

T =  2*pi/omega_ext; % time period of forcing

ndofs = length( BeamMesh.EBC.unconstrainedDOFs );
n_VMs = ndofs; % first n_VMs modes with lowest frequency calculated 
Kc = BeamAssembly.constrain_matrix(K);
Mc = BeamAssembly.constrain_matrix(M);
[V0,om2] = eigs(Kc, Mc, n_VMs, 'SM');
V0 = BeamAssembly.unconstrain_vector(V0);
[om, ~] = sort(sqrt(diag(om2)));
f0 = om/2/pi;
% Damping _____________________________________________________________
Qfactor = 100;
csi = 1./(2*Qfactor);   % dimensionless damping
om0 = 2*pi*f0(1);       % first eigenfrequency
alfa = 2 * om0 * csi;
D = alfa*M;             % Mass-proportinal damping: D = a*M
BeamAssembly.DATA.C = D;
BeamAssembly.DATA.D = D;

% % load amplification factor
% amplification_factor = 1;
% 
% % forcing function
% Pressure = 1; % in Pascals 
% F = Pressure*BeamAssembly.uniform_body_force();
% F_ext = @(t) amplification_factor * F * sin(omega_ext * t);
%% Reduced basis
m = 3; % use the first five VMs in reduction
V = V0(:,1:m);
Vc = BeamAssembly.constrain_vector(V0(:,1:m));
Kc = BeamAssembly.constrain_matrix(K);
Mc = BeamAssembly.constrain_matrix(M);
figure
for i=1:m
    LHS_matrix=[full(Kc-omega(i)^2*Mc) -Mc*Vc(:,i);
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

%%
% External force __________________________________________________________
forced_dof = round(nNodes/2) * 3 - 1; % middle node, vertical dof
F = zeros(BeamMesh.nDOFs, 1);
F( forced_dof ) = 1;
Famp = 1;
F_ext = @(t) Famp * F * sin(omega_ext * t);

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
tmax = 260*T; 
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
plot(TI_lin.Solution.time, TI_lin.Solution.u(forced_dof,:)/height,"DisplayName","linear full")
plot(TI_NL.Solution.time, TI_NL.Solution.u(forced_dof,:)/height,"DisplayName","nonlinear full")
xlabel('time'); 
% ylabel('q'); 
ylabel('|Q_1| / height [-]')
grid on;
legend


%% Reduced solution Linear
BeamReducedAssembly = ReducedAssembly(BeamMesh,V);

BeamReducedAssembly.DATA.M = BeamReducedAssembly.mass_matrix();
% BeamReducedAssembly.DATA.C = BeamReducedAssembly.damping_matrix();
BeamReducedAssembly.DATA.C = V'*D*V;
BeamReducedAssembly.DATA.K =  BeamReducedAssembly.stiffness_matrix();

q0 = zeros(m+m*(m+1)/2,1);
qd0 = zeros(m+m*(m+1)/2,1);
qdd0 = zeros(m+m*(m+1)/2,1);

TI_lin_red = ImplicitNewmark('timestep',dt,'alpha',0.005,'linear',true);

% Modal linear Residual evaluation function handle
Residual_lin_red = @(q,qd,qdd,t)residual_reduced_linear(q,qd,qdd,t,BeamReducedAssembly,F_ext);

% time integration
TI_lin_red.Integrate(q0,qd0,qdd0,tmax,Residual_lin_red);
TI_lin_red.Solution.u = V * TI_lin_red.Solution.q;

%% Reduced solution Noninear
% For demonstration purposes, we simply reduce the nonlinear system using
% out-of-plane bending modes. This is expected to produce bad results when 
% in-plane stretching is involved in the response.


TI_NL_alpha_red = GeneralizedAlpha('timestep',dt,'rho_inf',0.7);

% Modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear(q,qd,qdd,t,BeamReducedAssembly,F_ext);

% time integration
TI_NL_alpha_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_alpha_red.Solution.u = V * TI_NL_alpha_red.Solution.q;

%% plotting


figure;hold on
plot(TI_lin_red.Solution.time, TI_lin_red.Solution.u(forced_dof,:)/height,"DisplayName","linear red")
plot(TI_NL_alpha_red.Solution.time, TI_NL_alpha_red.Solution.u(forced_dof,:)/height,"DisplayName","nonlinear red")
xlabel('time'); 
% ylabel('q'); 
ylabel('|Q_1| / height [-]')
grid on;
legend
