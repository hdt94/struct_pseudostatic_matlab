classdef Frame < handle & matlab.mixin.Heterogeneous
    % Frame element. 
    %
    % Obj = Frame(Material, Section, Loads)
    %    
    % Stiffness:
    %   Frame stiffness is compounded by following matrices which are defined as
    %   properties:    
    %       k = [kiij, kij; kji, kjji];
    %   Which, according to Kanderstuncer conventions, comply following
    %   equivalences:
    %       kiij = kjji
    %       kij = kji    
    %
    % Loads:
    %   Positive values of loads are supposed to be in a negative ortogonal
    %   direction. Example considering i < j:
    %
    %                P = 4      w = 0.5
    %                  |        ||||||
    %                  v        vvvvvv
    %          o-----------------------------o
    %       Node i                          Node j
    %
    
    properties (SetAccess = private)
        % stiff.Load derived array. Applied loads.
        Loads = stiff.Load.empty();
        %
        Material = stiff.Material.empty();
        
        Section = stiff.Section.empty();
    end
    
    properties (SetAccess = private, Transient)
        % Numeric matrix (3,3). Stiffness in node i.
        kiij = [];
        % Positive numeric scalar. Length of frame.
        Length = [];
        % Numeric matrix (3,3). Rotation from local to global in i-node.
        rij = [];
        % Numeric matrix (3,3). Rotation from local to global in j-node.
        rji = [];
        
    end
    
    properties (Dependent)
        % Numeric matrix (3,3). Stiffness contribution of node i in node j.
        kij
        % Numeric matrix (3,3). Stiffness contribution of node j in node i.
        kji
        % Numeric matrix (3,3). Stiffness in node j.
        kjji
        
        Rij
        
        Rji
    end
    
    methods
        function Obj = Frame(Material, Section, Loads)
            validateattributes(Material, {'stiff.Material'}, {'scalar'}, '', 'Material');
            validateattributes(Section, {'stiff.Section'}, {'scalar'}, '', 'Section');
            Obj.Material = Material;
            Obj.Section = Section;            
            if nargin > 2
                validateattributes(Loads, {'stiff.Load'}, {'row'}, '', 'Loads');
                Obj.Loads = Loads;
            end
        end
        
        function kij = get.kij(Obj)
            kij = Obj.kiij;
            kij(end) = 0.5 * kij(end);
        end
        
        function kji = get.kji(Obj)
            kji = Obj.kiij;
            kji(end) = 0.5 * kji(end);
        end
        
        function kjji = get.kjji(Obj)
            kjji = Obj.kiij;
        end
        
        function Rij = get.Rij(Obj)
            Rij = Obj.rij';
        end
        
        function Rji = get.Rji(Obj)
            Rji = Obj.rji';
        end
        
        function [Local] = calculateBasicReactions(Obj)
            % Calculate local reactions in basic system
            %
            % Output:
            %   Local: numeric array (3,NFrames,2). Local reactions in basic system.
            %       Notes:
            %           Local(1,:,:): axial forces.
            %           Local(2,:,:): shear forces.
            %           Local(3,:,:): moments.
            %           Local(:,:,1): reactions in i-nodes.
            %           Local(:,:,2): reactions in j-nodes.
            %
            % Appendix:
            %   NFrames: number of frames.
            
            NDoF = 3; % Number of DoF in each frame node
            NFrames = numel(Obj);
            Local = zeros(NDoF, NFrames, 2);
            for f = 1 : NFrames
                NLoads = numel(Obj(f).Loads);
                iLoads = zeros(NLoads, 2);
                jLoads = zeros(NLoads, 2);
                for lo = 1 : NLoads
                    [iLoads(lo,:), jLoads(lo,:)] = Obj(f).Loads(lo).calculateReactions(Obj(f).Length);
                end
                Local(:,f,1) = [0, sum(iLoads,1)];
                Local(:,f,2) = [0, sum(jLoads,1)];
            end
        end % function [...] = calculateBasicReactions(...)
        
        function defineLengthNRotation(Obj, Nodes, Connectivity)
            % Define frame length and rotation matrices
            %
            % Obj.defineLengthNRotation(Nodes, Connectivity);
            %
            % Inputs:
            %   Nodes: stiff.Node array. Nodes in which connectivity depends.
            %   Connectivity: positive integer matrix (NFrames,2). Node connectivity index of frames.
            %
            % Note: maximum index in Connectivity cannot be greater than NNodes.
            %
            % Appendix:
            %   NFrames: number of frames.
            %   NNodes: number of nodes.
            
            Location = [Nodes(:).Location];
            [Length, rij, rji] = Obj.calculateLengthNRotation(Location, Connectivity);
            for f = 1 : numel(Obj)
                Obj(f).Length = Length(f);
                Obj(f).rij = rij(:,:,f);
                Obj(f).rji = rji(:,:,f);
            end
        end
        
        function defineStiffness(Obj)
            % Define frame stiffness
            %
            % Obj.defineStiffness();
            %
            % Note: previous call to defineLengthNRotation() method is required.
            %
            % Example:
            %   Obj.defineLengthNRotation();
            %   Obj.defineStiffness();
            
            Sections = [Obj(:).Section];
            A = cat(3,Sections(:).Area) .* cat(3,Sections(:).Modifier_Area);
            I = cat(3,Sections(:).Inertia);
            Materials = [Obj(:).Material];
            E = cat(3,Materials(:).Young);
            Length = cat(3,Obj(:).Length);
            kiij = Obj.calculateStiffness(A, I, E, Length);
            for f = 1 : numel(Obj)
                Obj(f).kiij = kiij(:,:,f);
            end
        end % function defineStiffness(...)
        
        function [Masses] = getLoadMasses(Obj)
            %
            % All loads are supposed to be uniform distributed loads
            %
            %
            % Output:
            %   
            
            NFrames = numel(Obj);
            Weights = zeros(1, NFrames);
            for f = 1 : NFrames
                Weights(:,f) = Obj(f).Length * sum([Obj(f).Loads(:).Value]);                
            end
            Gravity = stiff.Definition.Constant.Gravity;
            Masses = Weights * (1/Gravity);
        end
        
        function [Masses] = getSelfMasses(Obj)            
            Lengths = [Obj(:).Length];
            Materials = [Obj(:).Material];
            Sections = [Obj(:).Section];
            Masses = Lengths .* [Sections(:).Area] .* [Materials(:).Density];
        end        
        
        function plotDiagrams(Obj, Ax, Connectivity, Nodes, Loads)
            %
            % Input:
            %   Ax:
            %
            %   Connectivity:
            %   Location:
            %   Loads: numeric array (NDoF,NFrames,2)
            %
            %

            for n = 1 : numel(Ax)
                colormap(Ax(n), [1,1,0;1,0,0]);
            end
            
            NFrames = numel(Obj);
            Location = [Nodes(:).Location];
            iLoads = Loads(:,:,1);
            jLoads = Loads(:,:,2);
            Max = max(abs([iLoads,jLoads]), [], 2);
            Factors = 0.1 * max([Obj(:).Length]) ./ Max;              
            
            % Parameters for normal-force diagram
            NormalFace = reshape(1:4*NFrames, 4, NFrames)';
            NormalColor = zeros(NFrames, 1);
            Logical = iLoads(1,:) >= 0;
            NormalColor(Logical) = 0; % Color for positive normal-forces
            Logical = not(Logical);
            NormalColor(Logical) = 1; % Color for negative normal-forces
            NormalVertex = zeros(4*NFrames, 3);
            Counter1 = 0;
            
            GProps = {'FaceAlpha',0.8};%'FaceColor','flat', 'FaceAlpha',1};
            
            GShears = gobjects(1, NFrames);
            GMoments = gobjects(1, NFrames);
            for f = 1 : NFrames
                iLoc = Location(:,Connectivity(f,1));
                jLoc = Location(:,Connectivity(f,2));
                
                % Normal-diagram calculation
                iRotation = [Obj(f).rij(1,2); 0; Obj(f).rij(2,2)];
                jRotation = [Obj(f).rji(1,2); 0; Obj(f).rji(2,2)];
                iDelta = iRotation * Factors(1) * -iLoads(1,f); % Considering shear local-convention
                jDelta = jRotation * Factors(1) * jLoads(1,f);
                Index = (1:4) + Counter1; % Vector index
                NormalVertex(Index,:) = [iLoc, jLoc, jLoc+jDelta, iLoc+iDelta]';
                Counter1 = Counter1 + 4;
                
                % Shear and moment diagrams calculation
                XVector = linspace(0, Obj(f).Length, 100);
                NPoints = numel(XVector);
                [VShear, VMoment] = Obj(f).Loads(:).getDiagrams(XVector, iLoads(:,f));
                
                % Defining colors
                ShColor = zeros(NPoints, 1); % Colors of Shear vertices
                MoColor = zeros(NPoints, 1); % Colors of Moment vertices
                Logical = VShear >= 0;
                ShColor(Logical,:) = 1;
                ShColor(not(Logical)) = 0;
                Logical = VMoment >= 0;
                MoColor(Logical) = 1;
                MoColor(not(Logical)) = 0;
                
                % Calculating data
                jRotation = [Obj(f).rji(1,:); zeros(1,3); Obj(f).rji(2,:)];
                MShear = [XVector; Factors(2)*VShear; zeros(1,NPoints)];
                MMoment = [XVector; Factors(3)*VMoment; zeros(1,NPoints)];
                ShData = zeros(NPoints,3);
                MoData = zeros(NPoints,3);
                for n = 1 : NPoints
                    ShData(n,:) = iLoc + jRotation*MShear(:,n);
                    MoData(n,:) = iLoc + jRotation*MMoment(:,n);
                end
                ShData = [iLoc'; ShData; jLoc'; iLoc'];
                MoData = [iLoc'; MoData; jLoc'; iLoc'];
                ShColor = [ShColor(1,:);ShColor;ShColor(end,:);ShColor(1,:)];
                MoColor = [MoColor(1,:);MoColor;MoColor(end,:);MoColor(1,:)];
                
                % Note:
                % fill3(..., 'Parent',Ax) is used instead of fill3(Ax, ...) to make
                % it compatible with previous versions to MATLAB 2016a
                GShears(f) = fill3(ShData(:,1), ShData(:,2), ShData(:,3), ShColor, GProps{:}, 'Parent',Ax(2));
                GMoments(f) = fill3(MoData(:,1), MoData(:,2), MoData(:,3), MoColor, GProps{:}, 'Parent',Ax(3));
            end
            %             GNormal = patch(Ax(1), 'Faces',NormalFace,'Vertices',NormalVertex,'FaceVertexCData',NormalColor, GProps{:});
        end        
    end
    
    methods (Static)
        function [Lengths, rij, rji] = calculateLengthNRotation(Location, Connectivity)
            % Calculate length and rotation matrices
            %
            % [Lengths, rij, rji] = stiff.Frame.calculateLengthNRotation(Location, Connectivity);
            %
            % Input:
            %   Location: numeric matrix (3,NNodes). Node global location considering X, Y, Z.
            %   Connectivity: positive integer matrix (NFrames,2). Node connectivity index of frames.
            %
            % Output:
            %   Lengths: numeric vector (1,NFrames).
            %   rij: numeric matrix (3,3,NFrames). Rotation from local to global in i-node.
            %   rji: numeric matrix (3,3,NFrames). Rotation from local to global in j-node.
            %
            % Appendix:
            %   NNodes: number of nodes.
            %   NFrames: number of frames.
            %
            % Example:
            %   Location = [0,0,0; 0,0,3; 3,0,3; 3,0,0]';
            %   Connectivity = [1,2; 2,3; 3,4; 1,3];
            %	[Lengths, rij, rji] = stiff.Frame.calculateLengthNRotation(Location, Connectivity);
            
            NFrames = size(Connectivity, 1);
            
            % Calculating lengths and angles
            iNode = Connectivity(:,1);
            jNode = Connectivity(:,2);
            XYZDiff = Location(:,jNode) - Location(:,iNode);
            Lengths = sum(XYZDiff.^2,1) .^ 0.5;
            Angle = acos(XYZDiff(1,:) ./ Lengths);
            Logical = XYZDiff(end,:) < 0;
            Angle(Logical) = 2*pi - Angle(Logical);
            Angle = reshape(Angle, 1, 1, NFrames);
            
            % Calculating components of rotation matrices (view from i-node)
            cosXx = - cos(Angle);
            cosZy = cosXx;
            cosXy = sin(Angle);
            cosZx = - cosXy;
            
            % Assembling rotation matrices
            Zeros = zeros(1, 1, NFrames);
            Ones = ones(1, 1, NFrames);
            rij = [
                cosXx	cosXy   Zeros
                cosZx   cosZy   Zeros
                Zeros   Zeros   -Ones
                ];
            rji = [
                -cosXx	-cosXy  Zeros
                -cosZx  -cosZy	Zeros
                Zeros   Zeros	-Ones
                ];
        end % function [...] = calculateLengthNRotation(...)
        
        function [kiij] = calculateStiffness(A, I, E, Length)
            % Calculate local stiffness matrices
            %
            % kiij = stiff.Frame.calculateStiffness(A, I, E, Length);
            %
            % Input:
            %   A: numeric array (1,1,NFrames). Frame areas.
            %   I: numeric array (1,1,NFrames). Frame intertias.
            %   E: numeric array (1,1,NFrames). Frame elasticity moduli.
            %   Length: numeric array (1,1,NFrames). Frame lengths.
            %
            % Output:
            %   kiij: numeric array (3,3,NFrames). Frame self-node stiffness.
            %
            % Appendix:
            %   NFrames: number of frames.
            %
            % Example:
            %   A = rand(1, 1, 4);
            %   I = rand(1, 1, 4);
            %   E = rand(1, 1, 4);
            %   Length = rand(1, 1, 4);
            %   kiij = stiff.Frame.calculateStiffness(A, I, E, Length);
            
            NFrames = size(E, 3);
            
            % Calculating components of stiffness matrix
            Zeros = zeros(1, 1, NFrames);
            EI = E .* I;
            k11 = E .* A ./ Length;
            k22 = 12 * EI ./ Length.^3;
            k23 = - 6 * EI ./ Length.^2;
            k33 = 4 * EI ./ Length;
            
            % Assembling stiffness matrix considering k32 = k23
            kiij = [
                k11     Zeros   Zeros
                Zeros   k22     k23
                Zeros	k23     k33];
        end % function [...] = calculateStiffness(...)
    end % methods (Static)
end % classdef