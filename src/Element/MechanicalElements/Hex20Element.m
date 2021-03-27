classdef Hex20Element < Element
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 3     % number of DOFs per node
        nNodes = 20         % number of nodes per element
        nDim = 3            % number of dimensions in local coordinates
    end
    
    properties
        quadrature              % weights and points for gauss quadrature
        Material                % Object of class Material
        initialization          % some 0-matrices to speedup numerical integration
    end
    
    properties (Dependent)
        uniformBodyForce
        vol                     % volume of the element
    end
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function self = Hex20Element(Material,Ngauss)
            % _____________________________________________________________
            %
            % SELF-FUNCTION
            % self = Hex20Element(Material,Ngauss)
            % defines element's properties
            %______________________________________________________________
            self.Material = Material;
            if nargin == 1
                Ngauss = 2;
            end
            [x,w]=lgwt(Ngauss,-1,1);
            self.quadrature.Ng = Ngauss;
            self.quadrature.X = x;	% gauss integration points
            self.quadrature.W = w;	% gauss integration weights
            % INIZIALIZATION of some matrices (this should speedup
            % numerical integration)
            C = self.Material.get_stress_strain_matrix_3D;
            H = [1 0 0 0 0 0 0 0 0;
                0 0 0 0 1 0 0 0 0;
                0 0 0 0 0 0 0 0 1;
                0 1 0 1 0 0 0 0 0;
                0 0 1 0 0 0 1 0 0;
                0 0 0 0 0 1 0 1 0];
            self.initialization.A = zeros(6,9); % nonlinear strain
            self.initialization.G = zeros(9,60);% shape function derivatives
            self.initialization.Z = zeros(20);  % zero-matrix
            self.initialization.K = zeros(60);  % stiffness-element matrix
            self.initialization.F = zeros(60,1);% internal forces (element)
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
            Mel = zeros(60);
            for ii = 1:self.quadrature.Ng
                for jj = 1:self.quadrature.Ng
                    for kk = 1:self.quadrature.Ng
                        g = X(ii);
                        h = X(jj);
                        r = X(kk);
                        N = self.shape_functions(g,h,r);
                        [~,detJ] = shape_function_derivatives(self,g,h,r);
                        NN(1,1:3:60) = N;
                        NN(2,2:3:60) = N;
                        NN(3,3:3:60) = N;
                        % integration of K and M through GAUSS QUADRATURE
                        Mel = Mel + W(ii)*W(jj)*W(kk)*(NN'*NN)*detJ;
                    end
                end
            end
            Mel = sparse(rho*Mel);
        end
        
        function [K,F] = tangent_stiffness_and_force(self,x)
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
                    for kk = 1:self.quadrature.Ng
                        g = X(ii);              % natural coordinates
                        h = X(jj);              % natural coordinates
                        r = X(kk);              % natural coordinates
                        we = W(ii)*W(jj)*W(kk); % weights
                        [G,detJ,dH] = shape_function_derivatives(self,g,h,r);
                        th  = G*displ;
                        A = self.initialization.A;
                        A(1,1)=th(1); A(4,1)=th(2); A(5,1)=th(3); A(2,2)=th(2); A(4,2)=th(1);
                        A(6,2)=th(3); A(3,3)=th(3); A(5,3)=th(1); A(6,3)=th(2); A(1,4)=th(4);
                        A(4,4)=th(5); A(5,4)=th(6); A(2,5)=th(5); A(4,5)=th(4); A(6,5)=th(6);
                        A(3,6)=th(6); A(5,6)=th(4); A(6,6)=th(5); A(1,7)=th(7); A(4,7)=th(8);
                        A(5,7)=th(9); A(2,8)=th(8); A(4,8)=th(7); A(6,8)=th(9); A(3,9)=th(9);
                        A(5,9)=th(7); A(6,9)=th(8);
                        % Green Strain tensor
                        E = (H + 1/2*A)*th;
                        % second Piola-Kirchhoff stress tensor
                        s = C*E; % s = [S11 S22 S33 S12 S13 S23]
                        S = [s(1) s(4) s(5); s(4) s(2) s(6); s(5) s(6) s(3)];
                        Bnl = (H + A)*G;
                        % functions to integrate over volume
                        int_K1 = Bnl'*C*Bnl;
                        HSH = dH'*S*dH;
                        int_Ks = [HSH ZZ ZZ; ZZ HSH ZZ; ZZ ZZ HSH]; % (faster than blkdiag)
                        int_K = (int_K1 + int_Ks)*detJ;
                        int_F = (Bnl'*s)*detJ;
                        % integration of K and F through Gauss quadrature
                        K = K + we*int_K;
                        F = F + we*int_F;
                    end
                end
            end
        end
        
        function xe = extract_element_data(self,x)
            % x is a vector of full DOFs
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            xe = x(index,:);
        end
        
        function  f = get.uniformBodyForce(self)
            % _____________________________________________________________
            %
            % F = uniform_body_force(self,direction)
            % This function computes a load along direction=3(Z) by
            % dividing the load on the 20 nodes according to the element
            % volume (V/20) [it might not be the best way, but still...]
            %______________________________________________________________
            f = sparse(60,1);
            f(3:3:end) = self.vol/20; % uniformly distributed pressure on the structure
        end
        
        % ANCILLARY FUNCTIONS _____________________________________________
        
        function V = get.vol(self)
            % volume is given by the integral of detJ (jacobian from 
            % isoparametric to physical space) over the volume of the
            % isoparametric element
            detJ = 0;
            W = self.quadrature.W;
            X = self.quadrature.X;
            for ii = 1 : length( w )
                g = X(1,ii);
                h = X(2,ii);
                r = X(3,ii);
                [~, detJ_i] = shape_function_derivatives(self,g,h,r);
                detJ = detJ + detJ_i * W(ii);
            end
            V = detJ;
        end
        
        function [G,detJ,dH] = shape_function_derivatives(self,g,h,r)
            %______________________________________________________________
            %
            % [G,detJ,dH] = G_HEX20(self,g,h,r)
            % G = shape function derivatives in physical coordinates, s.t.
            % th=G*p with th={ux uy uz vx vy vz wx wy wz}' (ux=du/dx...)
            % and p={u1,v1,w1,...,u20,v20,w20}' (nodal displacements).
            % detJ = det(J), J=jacobian
            %______________________________________________________________
            xyz = self.nodes;
            % shape functions derivatives (ABAQUS ORDER)
            dHn = [ ...
                (g/8 - 1/8)*(h - 1)*(r - 1) + ((h - 1)*(r - 1)*(g + h + r + 2))/8,   (g/8 + 1/8)*(h - 1)*(r - 1) - ((h - 1)*(r - 1)*(h - g + r + 2))/8, - (g/8 + 1/8)*(h + 1)*(r - 1) - ((h + 1)*(r - 1)*(g + h - r - 2))/8, - (g/8 - 1/8)*(h + 1)*(r - 1) - ((h + 1)*(r - 1)*(g - h + r + 2))/8, - (g/8 - 1/8)*(h - 1)*(r + 1) - ((h - 1)*(r + 1)*(g + h - r + 2))/8, - (g/8 + 1/8)*(h - 1)*(r + 1) - ((h - 1)*(r + 1)*(g - h + r - 2))/8, (g/8 + 1/8)*(h + 1)*(r + 1) + ((h + 1)*(r + 1)*(g + h + r - 2))/8, (g/8 - 1/8)*(h + 1)*(r + 1) + ((h + 1)*(r + 1)*(g - h - r + 2))/8, - (g/4 - 1/4)*(h - 1)*(r - 1) - ((g + 1)*(h - 1)*(r - 1))/4,                               (h/4 - 1/4)*(h + 1)*(r - 1), (g/4 - 1/4)*(h + 1)*(r - 1) + ((g + 1)*(h + 1)*(r - 1))/4,                                -(h/4 - 1/4)*(h + 1)*(r - 1), (g/4 - 1/4)*(h - 1)*(r + 1) + ((g + 1)*(h - 1)*(r + 1))/4,                                -(h/4 - 1/4)*(h + 1)*(r + 1), - (g/4 - 1/4)*(h + 1)*(r + 1) - ((g + 1)*(h + 1)*(r + 1))/4,                               (h/4 - 1/4)*(h + 1)*(r + 1),                                -(r/4 - 1/4)*(h - 1)*(r + 1),                               (r/4 - 1/4)*(h - 1)*(r + 1),                                -(r/4 - 1/4)*(h + 1)*(r + 1),                               (r/4 - 1/4)*(h + 1)*(r + 1);
                (g/8 - 1/8)*(h - 1)*(r - 1) + (g/8 - 1/8)*(r - 1)*(g + h + r + 2), - (g/8 + 1/8)*(h - 1)*(r - 1) - (g/8 + 1/8)*(r - 1)*(h - g + r + 2), - (g/8 + 1/8)*(h + 1)*(r - 1) - (g/8 + 1/8)*(r - 1)*(g + h - r - 2),   (g/8 - 1/8)*(h + 1)*(r - 1) - (g/8 - 1/8)*(r - 1)*(g - h + r + 2), - (g/8 - 1/8)*(h - 1)*(r + 1) - (g/8 - 1/8)*(r + 1)*(g + h - r + 2),   (g/8 + 1/8)*(h - 1)*(r + 1) - (g/8 + 1/8)*(r + 1)*(g - h + r - 2), (g/8 + 1/8)*(h + 1)*(r + 1) + (g/8 + 1/8)*(r + 1)*(g + h + r - 2), (g/8 - 1/8)*(r + 1)*(g - h - r + 2) - (g/8 - 1/8)*(h + 1)*(r + 1),                                -(g/4 - 1/4)*(g + 1)*(r - 1), (h/4 - 1/4)*(g + 1)*(r - 1) + ((g + 1)*(h + 1)*(r - 1))/4,                               (g/4 - 1/4)*(g + 1)*(r - 1), - (h/4 - 1/4)*(g - 1)*(r - 1) - ((g - 1)*(h + 1)*(r - 1))/4,                               (g/4 - 1/4)*(g + 1)*(r + 1), - (h/4 - 1/4)*(g + 1)*(r + 1) - ((g + 1)*(h + 1)*(r + 1))/4,                                -(g/4 - 1/4)*(g + 1)*(r + 1), (h/4 - 1/4)*(g - 1)*(r + 1) + ((g - 1)*(h + 1)*(r + 1))/4,                                -(r/4 - 1/4)*(g - 1)*(r + 1),                               (r/4 - 1/4)*(g + 1)*(r + 1),                                -(r/4 - 1/4)*(g + 1)*(r + 1),                               (r/4 - 1/4)*(g - 1)*(r + 1);
                (g/8 - 1/8)*(h - 1)*(r - 1) + (g/8 - 1/8)*(h - 1)*(g + h + r + 2), - (g/8 + 1/8)*(h - 1)*(r - 1) - (g/8 + 1/8)*(h - 1)*(h - g + r + 2),   (g/8 + 1/8)*(h + 1)*(r - 1) - (g/8 + 1/8)*(h + 1)*(g + h - r - 2), - (g/8 - 1/8)*(h + 1)*(r - 1) - (g/8 - 1/8)*(h + 1)*(g - h + r + 2),   (g/8 - 1/8)*(h - 1)*(r + 1) - (g/8 - 1/8)*(h - 1)*(g + h - r + 2), - (g/8 + 1/8)*(h - 1)*(r + 1) - (g/8 + 1/8)*(h - 1)*(g - h + r - 2), (g/8 + 1/8)*(h + 1)*(r + 1) + (g/8 + 1/8)*(h + 1)*(g + h + r - 2), (g/8 - 1/8)*(h + 1)*(g - h - r + 2) - (g/8 - 1/8)*(h + 1)*(r + 1),                                -(g/4 - 1/4)*(g + 1)*(h - 1),                               (h/4 - 1/4)*(g + 1)*(h + 1),                               (g/4 - 1/4)*(g + 1)*(h + 1),                                -(h/4 - 1/4)*(g - 1)*(h + 1),                               (g/4 - 1/4)*(g + 1)*(h - 1),                                -(h/4 - 1/4)*(g + 1)*(h + 1),                                -(g/4 - 1/4)*(g + 1)*(h + 1),                               (h/4 - 1/4)*(g - 1)*(h + 1), - (r/4 - 1/4)*(g - 1)*(h - 1) - ((g - 1)*(h - 1)*(r + 1))/4, (r/4 - 1/4)*(g + 1)*(h - 1) + ((g + 1)*(h - 1)*(r + 1))/4, - (r/4 - 1/4)*(g + 1)*(h + 1) - ((g + 1)*(h + 1)*(r + 1))/4, (r/4 - 1/4)*(g - 1)*(h + 1) + ((g - 1)*(h + 1)*(r + 1))/4];
            % jacobian
            J = dHn*xyz;
            J1 = [0 0 0; 0 0 0; 0 0 0];
            J1(1,1) = J(2,2)*J(3,3) - J(2,3)*J(3,2);
            J1(2,1) = J(2,3)*J(3,1) - J(2,1)*J(3,3);
            J1(3,1) = J(2,1)*J(3,2) - J(2,2)*J(3,1);
            J1(1,2) = J(1,3)*J(3,2) - J(1,2)*J(3,3);
            J1(2,2) = J(1,1)*J(3,3) - J(1,3)*J(3,1);
            J1(3,2) = J(1,2)*J(3,1) - J(1,1)*J(3,2);
            J1(1,3) = J(1,2)*J(2,3) - J(1,3)*J(2,2);
            J1(2,3) = J(1,3)*J(2,1) - J(1,1)*J(2,3);
            J1(3,3) = J(1,1)*J(2,2) - J(1,2)*J(2,1);
            detJ = J(1,1)*J1(1,1) + J(1,2)*J1(2,1) + J(1,3)*J1(3,1);
            J1 = J1/detJ;
            dH = J1*dHn;        % derivatives in physical coordinates,
            % 3x20 matrix, [dNi_dx; dNi_dy; dNi_dz]
            % with i = 1...20
            G = self.initialization.G;
            G(1:3,1:3:60) = dH;
            G(4:6,2:3:60) = dH;
            G(7:9,3:3:60) = dH;
        end
        
    end % methods
    
    methods (Static)
        
        function N = shape_functions(g,h,r)
            % N = shape_functions(g,h,r)
            % SHAPE FUNCTIONS FOR A 20-NODED HEXAHEDRON
            % - see Abaqus documentation:
            % Abaqus theory guide > Elements > Continuum elements > ...
            % ... > 3.2.4 Solid isoparametric quadrilaterals and hexahedra
            N = [ -1/8*(1-g)*(1-h)*(1-r)*(2+g+h+r); % n1
                -1/8*(1+g)*(1-h)*(1-r)*(2-g+h+r);   % n2
                -1/8*(1+g)*(1+h)*(1-r)*(2-g-h+r);   % n3
                -1/8*(1-g)*(1+h)*(1-r)*(2+g-h+r);   % n4
                -1/8*(1-g)*(1-h)*(1+r)*(2+g+h-r);   % n5
                -1/8*(1+g)*(1-h)*(1+r)*(2-g+h-r);   % n6
                -1/8*(1+g)*(1+h)*(1+r)*(2-g-h-r);   % n7
                -1/8*(1-g)*(1+h)*(1+r)*(2+g-h-r);   % n8
                +1/4*(1-g)*(1+g)*(1-h)*(1-r);       % n9
                +1/4*(1-h)*(1+h)*(1+g)*(1-r);       % n10
                +1/4*(1-g)*(1+g)*(1+h)*(1-r);       % n11
                +1/4*(1-h)*(1+h)*(1-g)*(1-r);       % n12
                +1/4*(1-g)*(1+g)*(1-h)*(1+r);       % n13
                +1/4*(1-h)*(1+h)*(1+g)*(1+r);       % n14
                +1/4*(1-g)*(1+g)*(1+h)*(1+r);       % n15
                +1/4*(1-h)*(1+h)*(1-g)*(1+r);       % n16
                +1/4*(1-r)*(1+r)*(1-g)*(1-h);       % n17
                +1/4*(1-r)*(1+r)*(1+g)*(1-h);       % n18
                +1/4*(1-r)*(1+r)*(1+g)*(1+h);       % n19
                +1/4*(1-r)*(1+r)*(1-g)*(1+h)];      % n20
        end
        
        function X = natural_coordinates
            X = [ ...
                -1 -1 -1
                +1 -1 -1
                +1 +1 -1
                -1 +1 -1
                -1 -1 +1
                +1 -1 +1
                +1 +1 +1
                -1 +1 +1
                 0 -1 -1
                +1  0 -1
                 0 +1 -1
                -1  0 -1
                 0 -1 +1
                +1  0 +1
                 0 +1 +1
                -1  0 +1
                -1 -1  0
                +1 -1  0
                +1 +1  0
                -1 +1  0];
        end
        
    end
    
end % classdef

