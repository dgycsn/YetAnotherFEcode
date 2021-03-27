classdef Quad8Element < Element
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 2     % number of DOFs per node
        nNodes = 8          % number of nodes per element
        nDim = 2            % number of dimensions in local coordinates
    end
    
    properties
        thickness = 0	% element thickness, by default zero
        Material       	% Object of class Material
        quadrature    	% weights and points for gauss quadrature
        initialization 	% some 0-matrices to speedup numerical integration     
    end
    
    properties (Dependent)
        uniformBodyForce
        area                % area of the element
    end
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function self = Quad8Element(thickness, Material, Ngauss)
            % _____________________________________________________________
            %
            % SELF-FUNCTION
            % self = Quad8Element(Material,Ngauss)
            % defines element's properties
            %______________________________________________________________
            self.Material = Material;
            if nargin == 2
                Ngauss = 2;
            end
            [x, w] = lgwt(Ngauss,-1,1);
            self.quadrature.Ng = Ngauss;
            self.quadrature.X = x;      % gauss integration points
            self.quadrature.W = w;      % gauss integration weights
            self.thickness = thickness;
            % INIZIALIZATION of some matrices (this should speedup
            % numerical integration)
            C = self.Material.get_stress_strain_matrix_2D;
            H = [1 0 0 0; 
                0 0 0 1; 
                0 1 1 0];
            self.initialization.A = zeros(3,4); % nonlinear strain
            self.initialization.G = zeros(4,16);% shape function derivatives
            self.initialization.Z = zeros(8);   % zero-matrix
            self.initialization.K = zeros(16);  % stiffness-element matrix
            self.initialization.F = zeros(16,1);% internal forces (element)
            self.initialization.C = C;          % constitutive law matrix
            self.initialization.H = H;          % linear strain
        end
        
        function Mel = mass_matrix(self)
            % _____________________________________________________________
            %
            % Mel = mass_matrix_global(self,nodes,~);
            % Mel: element-level mass matrix (in global coordinates)
            %______________________________________________________________
            X = self.quadrature.X;
            W = self.quadrature.W;
            rho = self.Material.DENSITY;
            Mel = zeros(16);
           for ii = 1:self.quadrature.Ng
                for jj = 1:self.quadrature.Ng
                    g = X(ii);
                    h = X(jj);
                    % shape functions and detJ (ABAQUS ORDER)
                    N = self.shape_functions(g,h);
                    [~,detJ] = shape_function_derivatives(self,g,h);
                    NN = kron(N', eye(2));
                    % integration of K and M through GAUSS QUADRATURE
                    Mel = Mel + ( NN' * NN )*( W(ii) * W(jj) * detJ );
                end
            end
            Mel = sparse(rho*Mel);
        end
        
        function [K,F] = tangent_stiffness_and_force(self, x)
            displ = self.extract_element_data(x);
            X = self.quadrature.X;
            W = self.quadrature.W;
            K = self.initialization.K;
            F = self.initialization.F;
            C = self.initialization.C;
            H = self.initialization.H;
            ZZ = self.initialization.Z;
            for ii = 1:self.quadrature.Ng
                for jj = 1:self.quadrature.Ng
                    g = X(ii);
                    h = X(jj);
                    we = W(ii)*W(jj); % weights
                    [G,detJ,dH] = shape_function_derivatives(self,g,h);
                    th  = G*displ;
                    A =	[th(1)	0     th(3) 0;
                      	 0    	th(2) 0     th(4);
                         th(2)	th(1) th(4) th(3)];
                    % Green Strain tensor
                    E = (H + 1/2*A)*th;
                    % second Piola-Kirchhoff stress tensor
                    s = C*E; % s = [S11 S22 S12]
                    S = [s(1), s(3); s(3), s(2)];
                    Bnl = (H + A)*G;
                    % functions to integrate over volume
                    int_K1 = Bnl'*C*Bnl;
                    HSH = dH'*S*dH;
                    int_Ks = [HSH ZZ; ZZ HSH]; % (faster than blkdiag)
                    int_K = int_K1 + int_Ks;
                    int_F = Bnl'*s;
                    % integration of K and F through Gauss quadrature
                    K = K + int_K * (we * detJ);
                    F = F + int_F * (we * detJ);
                end
            end
        end
        
        function xe = extract_element_data(self, x)
            % x is a vector of full DOFs
            index = get_index(self.nodeIDs, self.nDOFPerNode);
            xe = x(index,:);
        end
        
        function f = uniform_body_force(self,~)
            % _____________________________________________________________
            %
            % F = uniform_body_force(self,direction)
            % This function computes a load along direction=2(Y) by
            % dividing the load on the 8 nodes according to the element
            % area (A/8) [it might not be the best way, but still...]
            %______________________________________________________________
            f = sparse(16,1);
            f(2:2:end) = self.area/8; % uniformly distributed pressure on the structure
        end
        
        % ANCILLARY FUNCTIONS _____________________________________________
        
        function A = get.area(self)
            % Integrate detJ (jacobian from isoparametric to physical
            % coordinates) over the area to get A
            detJ = 0;
            W = self.quadrature.W;
            X = self.quadrature.X;
            for ii = 1 : length( W )
                for jj = 1 : length( W )
                    g = X(ii);
                    h = X(ii);
                    [~, detJ_i] = shape_function_derivatives(self, g, h);
                    detJ = detJ + detJ_i * W(ii) * W(jj);
                end
            end
            A = detJ;
        end
        
        function [G,detJ,dH] = shape_function_derivatives(self,g,h)
            %______________________________________________________________
            %
            % [G,detJ,dH] = shape_function_derivatives(self,g,h,r)
            % G = shape function derivatives in physical coordinates, s.t.
            % th=G*p with th={ux uy vx vy}' (ux=du/dx...)
            % and p={u1,v1,...,u8,v8}' (nodal displacements).
            % detJ = det(J), J=jacobian
            %______________________________________________________________
            xy = self.nodes;
            % shape function derivatives in natural coordinates
            dHn = [ ...
                    -((h - 1)*(g + h + 1))/4 - ((g - 1)*(h - 1))/4, -((g - 1)*(g + h + 1))/4 - ((g - 1)*(h - 1))/4;
                    ((h - 1)*(h - g + 1))/4 - ((g + 1)*(h - 1))/4,   ((g + 1)*(h - g + 1))/4 + ((g + 1)*(h - 1))/4;
                    ((h + 1)*(g + h - 1))/4 + ((g + 1)*(h + 1))/4,   ((g + 1)*(g + h - 1))/4 + ((g + 1)*(h + 1))/4;
                    ((h + 1)*(g - h + 1))/4 + ((g - 1)*(h + 1))/4,   ((g - 1)*(g - h + 1))/4 - ((g - 1)*(h + 1))/4;
                                                        g*(h - 1),                                     g^2/2 - 1/2;
                                                      1/2 - h^2/2,                                      -h*(g + 1);
                                                       -g*(h + 1),                                     1/2 - g^2/2;
                                                      h^2/2 - 1/2,                                       h*(g - 1)]';
            J = dHn*xy;
            detJ = det(J);
            dH = J\dHn;	% derivatives in physical coordinates,
                      	% 2x16 matrix, [dNi_dx; dNi_dy]
                       	% with i = 1...10
            G = self.initialization.G;
            G(1:2,1:2:end) = dH;
            G(3:4,2:2:end) = dH;
        end
        
    end % methods
    
    methods (Static)
        
        function N = shape_functions(g,h)
            % N = shape_functions(g,h,r)
            % SHAPE FUNCTIONS FOR A 8-NODED QUADRILATERAL
            % - see Abaqus documentation:
            % Abaqus theory guide > Elements > Continuum elements > ...
            % ... > 3.2.4 Solid isoparametric quadrilaterals and hexahedra
            N = 1/4*[...
                    -(1-g)*(1-h)*(1+g+h); 
                    -(1+g)*(1-h)*(1-g+h);
                    -(1+g)*(1+h)*(1-g-h); 
                    -(1-g)*(1+h)*(1+g-h);
                    2*(1-g)*(1+g)*(1-h);  
                    2*(1-h)*(1+h)*(1+g);
                    2*(1-g)*(1+g)*(1+h);  
                    2*(1-h)*(1+h)*(1-g)];
        end
        
        function X = natural_coordinates
            X = [ ...
                -1  -1  % node 1 (corner)
                1   -1  % node 2 (corner)
                1   1	% node 3 (corner)
                -1  1   % node 4 (corner)
                0   -1  % node 5
                1   0   % node 6
                0   1   % node 7
                -1  0]; % node 8
        end
        
    end
        
end % classdef
    
