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

%% Eigenvalue problem
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

%% Reduced basis
m = 3; % use the first m VMs in reduction
V = V0(:,1:m);
Vc = BeamAssembly.constrain_vector(V0(:,1:m));
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

%% External force 

omega_ext=om0;

forced_dof = round(nNodes/2) * 3 - 1; % middle node, vertical dof
F = zeros(BeamMesh.nDOFs, 1);
F( forced_dof ) = 1;
Famp = 1;
F_ext = @(t) Famp * F * sin(omega_ext * t);

% time step for integration
T =  2*pi/omega_ext; % time period of forcing
dt = T/50/4;
tmax = 160*T; 

%% nonlin full integration 
% Instantiate object for nonlinear time integration
TI_NL = ImplicitNewmark('timestep',dt,'alpha',0.005);

u0=zeros( BeamMesh.nDOFs, 1);
q0 = BeamAssembly.constrain_vector(u0);
qd0 = BeamAssembly.constrain_vector(u0);
qdd0 = BeamAssembly.constrain_vector(u0);

% Linear Residual evaluation function handle
residual_NL_full = @(q,qd,qdd,t)residual_nonlinear(q,qd,qdd,t,BeamAssembly,F_ext);

% Nonlinear Time Integration
tic
TI_NL.Integrate(q0,qd0,qdd0,tmax,residual_NL_full);
TI_NL.Solution.u = BeamAssembly.unconstrain_vector(TI_NL.Solution.q);
time_nonlin_full=toc;

%% Reduced basis
% For demonstration purposes, we simply reduce the nonlinear system using
% out-of-plane bending modes. This is expected to produce bad results when 
% in-plane stretching is involved in the response.
BeamReducedAssembly = ReducedAssembly(BeamMesh,V);

BeamReducedAssembly.DATA.M = BeamReducedAssembly.mass_matrix();
% BeamReducedAssembly.DATA.C = BeamReducedAssembly.damping_matrix();
BeamReducedAssembly.DATA.C = V'*D*V;
BeamReducedAssembly.DATA.K =  BeamReducedAssembly.stiffness_matrix();

q0 = zeros(m+m*(m+1)/2,1);
qd0 = zeros(m+m*(m+1)/2,1);
qdd0 = zeros(m+m*(m+1)/2,1);
%% Reduced solution Nonlinear
TI_NL_red = ImplicitNewmark('timestep',dt,'alpha',0.005);

% Modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear(q,qd,qdd,t,BeamReducedAssembly,F_ext);
tic
% time integration
TI_NL_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_red.Solution.u = V * TI_NL_red.Solution.q;
time_nonlin_red=toc;

%% tensors
%different notations T2-->Q3 T3-->Q4
F2 = BeamAssembly.vector('F2',u0);
T2 = BeamAssembly.tensor('T2',[BeamMesh.nDOFs, BeamMesh.nDOFs, BeamMesh.nDOFs], [2,3]);

F3 = BeamAssembly.vector('F3',u0);
T3 = BeamAssembly.tensor('T3',[BeamMesh.nDOFs, BeamMesh.nDOFs, BeamMesh.nDOFs, BeamMesh.nDOFs], [2,3,4]);

T2_red=ttt(ttt(tensor(V),ttt(tensor(T2),tensor(V),3,1),1,1),tensor(V),2,1);
T3_red=ttt(ttt(ttt(tensor(V),ttt(tensor(T3),tensor(V),4,1),1,1),tensor(V),3,1),tensor(V),2,1);

T2t = T2_red + permute(T2_red, [1 3 2]); 
T3t = T3_red + permute(T3_red, [1 3 2 4]) + permute(T3_red, [1 4 2 3]);

tensors.M=BeamReducedAssembly.DATA.M;
tensors.C=BeamReducedAssembly.DATA.C;
tensors.K=BeamReducedAssembly.DATA.K;
tensors.T2=T2_red;
tensors.T3=T3_red;
tensors.T2t=T2t;
tensors.T3t=T3t;
%% nonlin integration tensors
% Instantiate object for nonlinear time integration
TI_NL_t = ImplicitNewmark('timestep',dt,'alpha',0.005);

% Linear Residual evaluation function handle
F_ext_red = @(t) Famp * V'*F * sin(omega_ext * t);
residual_NL_t = @(q,qd,qdd,t)residual_nonlinear_tensor(q,qd,qdd,t,tensors,F_ext_red);

% Nonlinear Time Integration
tic
TI_NL_t.Integrate(q0,qd0,qdd0,tmax,residual_NL_t);
TI_NL_t.Solution.u = V*TI_NL_t.Solution.q;
time_nonlin_tensor=toc;
%% plotting
figure;hold on
plot(TI_NL.Solution.time, TI_NL.Solution.u(forced_dof,:)/height,"DisplayName","nonlinear full")
plot(TI_NL_red.Solution.time, TI_NL_red.Solution.u(forced_dof,:)/height,"DisplayName","nonlinear red")
plot(TI_NL_t.Solution.time, TI_NL_t.Solution.u(forced_dof,:)/height,"DisplayName","nonlinear tensor")
xlabel('time'); 
ylabel('|Q_1| / height [-]')
grid on;
legend
ax=gca;
text(tmax/4,ax.YLim(2),"time full = " + num2str(time_nonlin_full))
text(tmax/4,ax.YLim(2)/2,"time red = " + num2str(time_nonlin_red))
text(tmax/4,0,"time tensor = " + num2str(time_nonlin_tensor))
text(tmax/2,ax.YLim(1)*3/4,"tensor error wrt full= " + ...
    num2str(norm(TI_NL.Solution.u(forced_dof,:)-TI_NL_t.Solution.u(forced_dof,:))))
text(tmax/2,ax.YLim(1)/2,"red error wrt full= " + ...
    num2str(norm(TI_NL.Solution.u(forced_dof,:)-TI_NL_red.Solution.u(forced_dof,:))))
