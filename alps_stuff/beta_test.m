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
% Qfactor = 100;
% csi = 1./(2*Qfactor);   % dimensionless damping
% om0 = 2*pi*f0(1);       % first eigenfrequency
% alfa = 2 * om0 * csi;
% 
% om1 = 2*pi*f0(2);       % second eigenfrequency
% beta = 2 * om1 * csi;
% 
% D = alfa*M+beta*K;             % Mass-proportinal damping: D = a*M
% BeamAssembly.DATA.C = D;

% RAYLEIGH DAMPING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% C = a*M + b*Kt
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
% show dimensionless damping
freq_range=0:2:500;
xi_a=a./(2*freq_range);
xi_b=b*freq_range/2;
% figure
% hold on
% plot(freq_range,xi_a)
% plot(freq_range,xi_b)
% plot(freq_range,xi_b+xi_a)
s=sprintf("freq_resp/"+...
    "bfreq"+num2str(frequenz(2)+...
    "Q1_"+num2str(Qfactors(1))+...
    "Q2_"+num2str(Qfactors(2))...
    ));
% saveas(gcf,s+".fig"); 
% saveas(gcf,s+".png");
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
D=a*M+b*K;
BeamAssembly.DATA.C=D;
Dc=BeamAssembly.constrain_matrix(D);
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
        title(append('d\phi_',num2str(i),'/dq_',num2str(j)))
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
%%
% f1=300:3:315;
% f2=316:1:322;
% f3=322.25:0.25:323;
% f4=323:0.1:323.5;
% f5=323.75;
% f6=324;
% f7=325:5:350;
% 
% freq_forward_nonlin=[f1 f2 f3 f4 f5 f6 f7];
% 
% f1=300:3:312;
% f2=313:1:319;
% f3=319.5:0.5:320.5;
% f4=320.6:0.1:321;
% f5=322:1:324;
% f6=325:5:350;
% 
% freq_backward_nonlin=flip([f1 f2 f3 f4 f5 f6]);

freq_forward_nonlin=300:3:450;
freq_backward_nonlin=flip(freq_forward_nonlin);

f1=300:3:310;
f2=310:1:315;
f3=315.25:0.25:318;
f4=319:1:324;
f5=325:5:350;
f6=360:10:450;

freq_forward_lin=[f1 f2 f3 f4 f5 f6];
freq_backward_lin=flip(freq_forward_lin);

force_amplitude=5;
%%
[lin_forward_full, lin_forward_red,lin_forward_analytical]=...
    freq_resp_1d_linear(freq_forward_lin,force_amplitude,...
    BeamAssembly,BeamReducedAssembly,"full", "reduced");
[lin_backward_full, lin_backward_red,lin_backward_analytical]=...
    freq_resp_1d_linear(freq_backward_lin,force_amplitude,...
    BeamAssembly,BeamReducedAssembly,"full" ,"reduced");
%%
[nl_backward_full,nl_backward_red,nl_backward_t]=...
    freq_resp_1d_nonlinear(freq_backward_nonlin,force_amplitude,...
    BeamAssembly,BeamReducedAssembly,tensors,"full" ,"reduced");
[nl_forward_full,nl_forward_red,nl_forward_t]=...
    freq_resp_1d_nonlinear(freq_forward_nonlin,force_amplitude,...
    BeamAssembly,BeamReducedAssembly,tensors,"full" ,"reduced");
total_time=toc;
%% plotting
figure;hold on

plot(freq_forward_lin,lin_forward_full/height,"DisplayName","lin forward full",...
    "LineStyle","--","LineWidth",2,"Color","r","Marker","o")
plot(freq_forward_nonlin,nl_forward_full/height,"DisplayName","nonlin forward full",...
    "LineStyle",":","LineWidth",2,"Color","r","Marker","s")
plot(freq_forward_lin,lin_forward_red/height,"DisplayName","lin forward red",...
    "LineStyle","--","LineWidth",2,"Color","r","Marker","x")
plot(freq_forward_nonlin,nl_forward_red/height,"DisplayName","nonlin forward red",...
    "LineStyle",":","LineWidth",2,"Color","r","Marker","d")
plot(freq_forward_nonlin,nl_forward_t/height,"DisplayName","nonlin forward t",...
    "LineStyle","-.","LineWidth",2,"Color","r","Marker","*")
plot(freq_forward_lin,lin_forward_analytical/height,"DisplayName","lin forward anlyt",...
    "LineStyle","--","LineWidth",2,"Color","r","Marker","h")


plot(freq_backward_lin,lin_backward_full/height,"DisplayName","lin backward full",...
    "LineStyle","--","LineWidth",2,"Color","b","Marker","o")
plot(freq_backward_nonlin,nl_backward_full/height,"DisplayName","nonlin backward full",...
    "LineStyle",":","LineWidth",2,"Color","b","Marker","s")
plot(freq_backward_lin,lin_backward_red/height,"DisplayName","lin backward red",...
    "LineStyle","--","LineWidth",2,"Color","b","Marker","x")
plot(freq_backward_nonlin,nl_backward_red/height,"DisplayName","nonlin backward red",...
    "LineStyle",":","LineWidth",2,"Color","b","Marker","d")
plot(freq_backward_nonlin,nl_backward_t/height,"DisplayName","nonlin backward t",...
    "LineStyle","-.","LineWidth",2,"Color","b","Marker","*")
plot(freq_backward_lin,lin_backward_analytical/height,"DisplayName","lin backward anlyt",...
    "LineStyle","--","LineWidth",2,"Color","b","Marker","h")

% load("nlvib_freq_resp_30_8750.mat")
% plot(nlvib_lin_freq,nlvib_lin_amp,"DisplayName","nlvib lin",...
%     "LineStyle","--","LineWidth",2,"Color","k")
% plot(nlvib_nonlin_freq,nlvib_nonlin_amp,"DisplayName","nlvib nonlin",...
%     "LineStyle",":","LineWidth",2,"Color","k")
f_s=freq_forward_lin(1);
f_e=freq_forward_lin(end);
xlim([f_s,f_e])
legend
saveas(gcf,s+".fig"); 
saveas(gcf,s+".png");