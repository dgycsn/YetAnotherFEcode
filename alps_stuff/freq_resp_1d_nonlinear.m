function [nonlin_full,nonlin_red,nonlin_tensor]=...
    freq_resp_1d_nonlinear(freq,amplification,BeamAssembly,BeamReducedAssembly,tensors,varargin)
%% variables for IC and return

fft_cycles=2;
nonlin_full=zeros(length(freq),length(amplification));
nonlin_red=zeros(length(freq),length(amplification));
nonlin_tensor=zeros(length(freq),length(amplification));

%first initial condition full
init_cond_full=zeros(BeamAssembly.Mesh.nDOFs, 1);

TI_NL_full.Solution.q=BeamAssembly.constrain_vector(init_cond_full);
TI_NL_full.Solution.qd=BeamAssembly.constrain_vector(init_cond_full);

%first initial condition red
V_size=size(BeamReducedAssembly.V);
init_cond_red=zeros(V_size(2),1);

TI_NL_red.Solution.q=init_cond_red;
TI_NL_red.Solution.qd=init_cond_red;

TI_NL_t.Solution.q=init_cond_red;
TI_NL_t.Solution.qd=init_cond_red;
%% freq response type

if isempty(varargin)
   istensor=true; 
   isfull=true;
   isreduced=true;
else
    varargin=string(varargin);
    if ismember("full",varargin)
        isfull=true;
    else 
        isfull=false;
    end
    if ismember("reduced",varargin)
        isreduced=true;
    else 
        isreduced=false;
    end
    if ismember("tensor",varargin)
        istensor=true;
    else 
        istensor=false;
    end
    if ismember("none",varargin)
        return
    end    
end
%% main loop
    for curr_freq_step=1:length(freq)
        for curr_amp_step=1:length(amplification)
        %% External force 

        omega_ext=freq(curr_freq_step)
        pause(5)
        
        forced_dof = round(BeamAssembly.Mesh.nNodes/2) * 3 - 1; % middle node, vertical dof
        F = zeros(BeamAssembly.Mesh.nDOFs, 1);
        F( forced_dof ) = 1;
        Famp = amplification(curr_amp_step);
        F_ext = @(t) Famp * F * sin(omega_ext * t);

        % time step for integration
        T =  2*pi/omega_ext; % time period of forcing
        dt = T/50/4;
        tmax = 160*T; 
        Fs=1/dt;

        %% nonlin full integration 
        if isfull
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
        
        last_cycles_sol=TI_NL_full.Solution.u(forced_dof,end-fft_cycles*T/dt:end);
        [~, Y, ~] = fft_yFs(last_cycles_sol, Fs);
        w=max(Y);
        nonlin_full(curr_freq_step,curr_amp_step)=norm(w);
        end
       
        %% Reduced solution Noninear
        % For demonstration purposes, we simply reduce the nonlinear system using
        % out-of-plane bending modes. This is expected to produce bad results when 
        % in-plane stretching is involved in the response.
        if isreduced
        q0 = TI_NL_red.Solution.q(:,end);
        qd0 = TI_NL_red.Solution.qd(:,end);
        qdd0 = init_cond_red;
        
        TI_NL_red = ImplicitNewmark('timestep',dt,'alpha',0.005);
        % Modal nonlinear Residual evaluation function handle
        Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear(q,qd,qdd,t,BeamReducedAssembly,F_ext);
   
        % time integration
        TI_NL_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
        TI_NL_red.Solution.u = BeamReducedAssembly.V * TI_NL_red.Solution.q;
        
        last_cycles_sol=TI_NL_red.Solution.u(forced_dof,end-fft_cycles*T/dt:end);
        [~, Y, ~] = fft_yFs(last_cycles_sol, Fs);
        w=max(Y);
        nonlin_red(curr_freq_step,curr_amp_step)=norm(w);

        end
        %% nonlin integration tensors
        if istensor
        q0 = TI_NL_t.Solution.q(:,end);
        qd0 = TI_NL_t.Solution.qd(:,end);
        qdd0 = init_cond_red;
        % Instantiate object for nonlinear time integration
        TI_NL_t = ImplicitNewmark('timestep',dt,'alpha',0.005);
        % Linear Residual evaluation function handle
        residual_NL_t = @(q,qd,qdd,t)residual_nonlinear_tensor(q,qd,qdd,t,tensors,F_ext);

        % Nonlinear Time Integration
        TI_NL_t.Integrate(q0,qd0,qdd0,tmax,residual_NL_t);
        TI_NL_t.Solution.u = BeamReducedAssembly.V*TI_NL_t.Solution.q;
        
        last_cycles_sol=TI_NL_t.Solution.u(forced_dof,end-fft_cycles*T/dt:end);
        [~, Y, ~] = fft_yFs(last_cycles_sol, Fs);
        w=max(Y);
        nonlin_tensor(curr_freq_step,curr_amp_step)=norm(w);
        end

        end
    end


end