function [lin_full,nonlin_full,lin_red,nonlin_red,nonlin_tensor]=...
    freq_resp_1d_new(freq,amplification)  
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
rho     = 8750;   % density
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
% ntot = BeamMesh.nDOFs;

%% Assembly full
BeamAssembly = Assembly(BeamMesh);
M = BeamAssembly.mass_matrix();
u0 = zeros( BeamMesh.nDOFs, 1);
[K,~] = BeamAssembly.tangent_stiffness_and_force(u0);
% store matrices
BeamAssembly.DATA.K = K;
BeamAssembly.DATA.M = M;

%% Eigenvalue problem
ndofs = length( BeamMesh.EBC.unconstrainedDOFs );

n_VMs = ndofs; % first n_VMs modes with lowest frequency calculated 
Kc = BeamAssembly.constrain_matrix(K);
Mc = BeamAssembly.constrain_matrix(M);
[V0,om2] = eigs(Kc, Mc, n_VMs, 'SM');
V0 = BeamAssembly.unconstrain_vector(V0);
[om, ~] = sort(sqrt(diag(om2)));
f0 = om/2/pi;

%% Damping
Qfactor = 100;
csi = 1./(2*Qfactor);   % dimensionless damping
om0 = 2*pi*f0(1);       % first eigenfrequency
alfa = 2 * om0 * csi;
D = alfa*M;             % Mass-proportinal damping: D = a*M
BeamAssembly.DATA.C = D;

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
%% Reduced assembly 
BeamReducedAssembly = ReducedAssembly(BeamMesh,V);

BeamReducedAssembly.DATA.M = BeamReducedAssembly.mass_matrix();
BeamReducedAssembly.DATA.C = V'*D*V;
BeamReducedAssembly.DATA.K =  BeamReducedAssembly.stiffness_matrix();
%% tensors
%different notations T2-->Q3 T3-->Q4
T2 = BeamAssembly.tensor('T2',[BeamMesh.nDOFs, BeamMesh.nDOFs, BeamMesh.nDOFs], [2,3]);
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
tensors.V=V;
%% variables for IC and return

lin_full=zeros(length(freq),length(amplification));
lin_red=zeros(length(freq),length(amplification));
nonlin_full=zeros(length(freq),length(amplification));
nonlin_red=zeros(length(freq),length(amplification));
nonlin_tensor=zeros(length(freq),length(amplification));

%first initial condition full
init_cond_full=zeros(BeamAssembly.Mesh.nDOFs, 1);

TI_lin_full.Solution.q=BeamAssembly.constrain_vector(init_cond_full);
TI_lin_full.Solution.qd=BeamAssembly.constrain_vector(init_cond_full);

TI_NL_full.Solution.q=BeamAssembly.constrain_vector(init_cond_full);
TI_NL_full.Solution.qd=BeamAssembly.constrain_vector(init_cond_full);

%first initial condition red
init_cond_red=zeros(m+m*(m+1)/2,1);

TI_lin_red.Solution.q=init_cond_red;
TI_lin_red.Solution.qd=init_cond_red;

TI_NL_red.Solution.q=init_cond_red;
TI_NL_red.Solution.qd=init_cond_red;

TI_NL_t.Solution.q=init_cond_red;
TI_NL_t.Solution.qd=init_cond_red;
    
%% main loop
    for curr_freq_step=1:length(freq)
        for curr_amp_step=1:length(amplification)
        %% External force 

        omega_ext=freq(curr_freq_step)
        pause(5)
        
        forced_dof = round(nNodes/2) * 3 - 1; % middle node, vertical dof
        F = zeros(BeamMesh.nDOFs, 1);
        F( forced_dof ) = 1;
        Famp = amplification(curr_amp_step);
        F_ext = @(t) Famp * F * sin(omega_ext * t);

        % time step for integration
        T =  2*pi/omega_ext; % time period of forcing
        dt = T/50/4;
        tmax = 160*T; 
        Fs=1/dt;
        %% Linear full
        % Initial condition: previous step
        q0 = TI_lin_full.Solution.q(:,end);
        qd0 = TI_lin_full.Solution.qd(:,end);
        qdd0 = BeamAssembly.constrain_vector(init_cond_full);

        % Instantiate object for linear time integration
        TI_lin_full = ImplicitNewmark('timestep',dt,'alpha',0.005,'linear',true);

        % Linear Residual evaluation function handle
        residual_lin = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,BeamAssembly,F_ext);
        TI_lin_full.Integrate(q0,qd0,qdd0,tmax,residual_lin);

        % obtain full solution
        TI_lin_full.Solution.u = BeamAssembly.unconstrain_vector(TI_lin_full.Solution.q);

        %% nonlin full integration 
        % Initial condition: previous step
        q0 = TI_NL_full.Solution.q(:,end);
        qd0 = TI_NL_full.Solution.qd(:,end);
        qdd0 = BeamAssembly.constrain_vector(init_cond_full);
        
        % Instantiate object for nonlinear time integration
        TI_NL_full = ImplicitNewmark('timestep',dt,'alpha',0.005);

        % Linear Residual evaluation function handle
        residual_NL_full = @(q,qd,qdd,t)residual_nonlinear(q,qd,qdd,t,BeamAssembly,F_ext);

        % Nonlinear Time Integration
        TI_NL_full.Integrate(q0,qd0,qdd0,tmax,residual_NL_full);
        TI_NL_full.Solution.u = BeamAssembly.unconstrain_vector(TI_NL_full.Solution.q);

        %% Linear reduced
        q0 = TI_lin_red.Solution.q(:,end);
        qd0 = TI_lin_red.Solution.qd(:,end);
        qdd0 = init_cond_red;

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

        q0 = TI_NL_red.Solution.q(:,end);
        qd0 = TI_NL_red.Solution.qd(:,end);
        qdd0 = init_cond_red;
        
        TI_NL_red = ImplicitNewmark('timestep',dt,'alpha',0.005);
        % Modal nonlinear Residual evaluation function handle
        Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear(q,qd,qdd,t,BeamReducedAssembly,F_ext);
   
        % time integration
        TI_NL_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
        TI_NL_red.Solution.u = V * TI_NL_red.Solution.q;

        %% nonlin integration tensors
        q0 = TI_NL_t.Solution.q(:,end);
        qd0 = TI_NL_t.Solution.qd(:,end);
        qdd0 = init_cond_red;
        % Instantiate object for nonlinear time integration
        TI_NL_t = ImplicitNewmark('timestep',dt,'alpha',0.005);
        % Linear Residual evaluation function handle
        residual_NL_t = @(q,qd,qdd,t)residual_nonlinear_tensor(q,qd,qdd,t,tensors,F_ext);

        % Nonlinear Time Integration
        TI_NL_t.Integrate(q0,qd0,qdd0,tmax,residual_NL_t);
        TI_NL_t.Solution.u = V*TI_NL_t.Solution.q;
        
        %% fft on last k cycles
        k=2;
        %lin full
        last_k_cycles_sol=TI_lin_full.Solution.u(forced_dof,end-k*T/dt:end);
        [~, Y, ~] = fft_yFs(last_k_cycles_sol, Fs);
        w=max(Y);
        lin_full(curr_freq_step,curr_amp_step)=norm(w)/height;
        %nl full
        last_k_cycles_sol=TI_NL_full.Solution.u(forced_dof,end-k*T/dt:end);
        [~, Y, ~] = fft_yFs(last_k_cycles_sol, Fs);
        w=max(Y);
        nonlin_full(curr_freq_step,curr_amp_step)=norm(w)/height;
        %lin red
        last_k_cycles_sol=TI_lin_red.Solution.u(forced_dof,end-k*T/dt:end);
        [~, Y, ~] = fft_yFs(last_k_cycles_sol, Fs);
        w=max(Y);
        lin_red(curr_freq_step,curr_amp_step)=norm(w)/height;
        %nl red
        last_k_cycles_sol=TI_NL_red.Solution.u(forced_dof,end-k*T/dt:end);
        [~, Y, ~] = fft_yFs(last_k_cycles_sol, Fs);
        w=max(Y);
        nonlin_red(curr_freq_step,curr_amp_step)=norm(w)/height;
        %nl tensor
        last_k_cycles_sol=TI_NL_t.Solution.u(forced_dof,end-k*T/dt:end);
        [~, Y, ~] = fft_yFs(last_k_cycles_sol, Fs);
        w=max(Y);
        nonlin_tensor(curr_freq_step,curr_amp_step)=norm(w)/height;
        end
    end


end