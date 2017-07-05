classdef Framed < stiff.Struct
    % Example:
    %   Obj = stiff.struct.Framed(Nodes, Location, Frames, Connectivity);
    
    properties (SetAccess = private)
        % Positive integer matrix (NFrames,2). Node connectivity of frames.
        Connectivity = [];
        % stiff.Frame array (1,NFrames). Frames of structure.
        Frames = stiff.Frame.empty;
    end
    
    properties (SetAccess = private, Transient)
        % Numeric vector (NDoF,1). Global displacements.
        Displacements = [];
        % Numeric vector (NDoF,1). Global loads.
        Loads = [];
        % Numeric matrix (NDoF,NDoF). Global stiffness.
        Stiffness = [];
    end
    
    methods
        function Obj = Framed(Nodes, Frames, Connectivity)
            Obj@stiff.Struct(Nodes);
            validateattributes(Frames, {'stiff.Frame'}, {'nonempty'}, '', 'Connectivity');
            validateattributes(Connectivity, {'numeric'}, {'integer','positive'}, '', 'Connectivity');
            if size(Connectivity,1) ~= numel(Frames)
                error('Framed:Framed:invalidConnectivityRows', 'Connectivity NRows must be equal to NFrames');
            elseif size(Connectivity,2) ~= 2
                error('Framed:Framed:invalidConnectivityCols', 'Connectivity NCols must be equal to 2');
            elseif numel(Obj.Nodes) < max(Connectivity(:))
                error('Framed:Framed:invalidConnectivityMax', 'Maximum value of Connectivity cannot be greater than NNodes');
            end
            Obj.Connectivity = sort(Connectivity, 2);
            Obj.Frames = Frames;
        end
        
        function defineStiffness(Obj)
            Obj.Frames.defineLengthNRotation(Obj.Nodes, Obj.Connectivity);
            Obj.Frames.defineStiffness();
            Obj.Stiffness = Obj.calculateStiffness(Obj.Frames, Obj.Connectivity);
        end
        
        function [FrameLoads] = solve(Obj)
            
            NDoF = size(Obj.Nodes(1).Restraints, 1);
            NNodes = numel(Obj.Nodes);
            
            
            Obj.defineStiffness(); % Stiffnesses
            
            % Loads
            %             [Loads] = stiff.struct.Framed.getLoadVector(Frames, Nodes, NodeIndex, Connectivity)
            
            LocalBasicReactions = Obj.Frames.calculateBasicReactions();
            BasicReactions = Obj.local2global(Obj.Frames, LocalBasicReactions);
            Obj.Loads = Obj.getNodeLoads(Obj.Nodes, 1:NNodes, BasicReactions, Obj.Connectivity);
            Obj.Loads = reshape(Obj.Loads, numel(Obj.Loads), 1);
            
            % Displacements
            Obj.Displacements = zeros(size(Obj.Loads));
            Obj.solveReducedStatic();
            
            % Resulting local loads in frames
            Disp = reshape(Obj.Displacements, NDoF, NNodes);
            FrameLoads = LocalBasicReactions + Obj.getFrameLocalLoads(Obj.Frames, Obj.Connectivity, Disp);
        end
        
        function solveReducedStatic(Obj)
            Restraints = [Obj.Nodes(:).Restraints];
            NodeIndex = find(not(all(Restraints,1)));
            Index = Obj.getSubIndex(NodeIndex, Restraints(:,NodeIndex));
            Obj.Displacements(Index) = Obj.Stiffness(Index,Index) \ Obj.Loads(Index);
        end
        
        function plot(Obj, Ax, FrameLoads)
            
            if nargin == 2
                hold(Ax, 'on');
                Obj.plotStruct(Ax, Obj.Nodes, Obj.Frames, Obj.Connectivity);
                view(Ax, 0, 0);
                hold(Ax, 'off');
            elseif numel(Ax) ~= 3
                error('Framed:plot:invalidAxNumber', 'Ax expected to be 3 numel');
            else
                for n = 1 : 3
                    hold(Ax(n), 'on');
                    view(Ax(n), 0, 0);
                end
                Obj.plotStruct(Ax, Obj.Nodes, Obj.Frames, Obj.Connectivity);
                Obj.Frames.plotDiagrams(Ax, Obj.Connectivity, Obj.Nodes, FrameLoads);
                for n = 1 : 3
                    hold(Ax(n), 'off');
                end
            end
        end % function plot(...)
    end % methods
    
    methods (Static)
        function [Loads] = getFrameLocalLoads(Frames, Connectivity, Disp)
            %
            %
            
            % Getting local displacements
            iNodeIndex = Connectivity(:,1);
            jNodeIndex = Connectivity(:,2);
            Global = cat(3, Disp(:,iNodeIndex), Disp(:,jNodeIndex));
            Local = stiff.struct.Framed.global2local(Frames, Global);
           
            % Calculating local frame loads
            NDoF = 3;
            NFrames = size(Connectivity, 1);
            Loads = zeros(NDoF, NFrames, 2);
            for f = 1 : NFrames
                kiij = Frames(f).kiij;
                kij = Frames(f).kij;
                kji = Frames(f).kji;
                kjji = Frames(f).kiij;                                
                iDisp = Local(:,f,1); % Local i-node displacements
                jDisp = Local(:,f,2); % Local j-node displacements
                Loads(:,f,1) = kiij*iDisp + kij*jDisp;
                Loads(:,f,2) = kji*iDisp + kjji*jDisp;
            end
        end
        
        
        function [Local] = global2local(Frames, Global)
            %
            % [Local] = stiff.struct.Framed.global2local(Frames, Global);
            
            Rij = permute(cat(3,Frames(:).rij), [2,1,3]); % Transformation global to local in i-node
            Rji = permute(cat(3,Frames(:).rji), [2,1,3]); % Transformation global to local in j-node
            NDoF = size(Rij, 1); % Number of degrees of freedom
            NFrames = numel(Frames);
            Local = zeros(NDoF, NFrames, 2);
            for f = 1 : NFrames
                Local(:,f,1) = Rij(:,:,f) * Global(:,f,1);
                Local(:,f,2) = Rji(:,:,f) * Global(:,f,2);
            end
        end
        
        function [Global] = local2global(Frames, Local)
            %
            
            rij = cat(3, Frames(:).rij); % Transformation local to global in i-node
            rji = cat(3, Frames(:).rji); % Transformation local to global in j-node
            NDoF = size(rij, 1); % Number of degrees of freedom
            NFrames = numel(Frames);
            Global = zeros(NDoF, NFrames, 2);
            for f = 1 : NFrames
                Global(:,f,1) = rij(:,:,f) * Local(:,f,1);
                Global(:,f,2) = rji(:,:,f) * Local(:,f,2);
            end
        end
        
        function [Loads] = getLoadVector(Frames, Nodes, NodeIndex, Connectivity)
            % [Loads] = stiff.struct.Framed.getLoadVector(Frames, Nodes, NodeIndex, Connectivity)
            LocalBasicReactions = Frames.calculateBasicReactions();
            BasicReactions = stiff.struct.Framed.local2global(Frames, LocalBasicReactions);
            Loads = stiff.struct.Framed.getNodeLoads(Nodes, NodeIndex, BasicReactions, Connectivity);
            Loads = reshape(Loads, numel(Loads), 1);
        end
        
        function [Loads] = getNodeLoads(Nodes, NodeIndex, Reactions, Connectivity)
            %
            %
            % Input:
            %   Nodes:
            %	NodeIndex:
            %   Reactions: numeric. Matrix (NDoF,NFrames,2) with frame reactions in
            %       basic system.
            %   Connectivity: numeric. Matrix (NFrames,2) with frame connectivity.
            %
            % Output:
            %   Loads: numeric. Matrix (NDoF,NNodes) with node loads.
            
            NDoF = size(Reactions, 1); % Number of degrees of freedom in each frame node
            NNodes = numel(NodeIndex);
            Loads = zeros(NDoF, NNodes);
            for n = 1 : NNodes
                Logical = Connectivity == NodeIndex(n);
                R1 = Reactions(:, Logical(:,1), 1);
                R2 = Reactions(:, Logical(:,2), 2);
                Loads(:,n) = Nodes(n).Loads - sum([R1,R2],2);
            end
        end
        
        function [Masses] = getNodeMasses(Nodes, NodeIndex, Frames, Connectivity)

            FrameMasses = 0.5 * (Frames.getLoadMasses() + Frames.getSelfMasses());
            NNodes = numel(NodeIndex);
            Masses = zeros(1, NNodes);
            for n = 1 : NNodes
                Logical = Connectivity == NodeIndex(n);
                M1 = FrameMasses(Logical(:,1));
                M2 = FrameMasses(Logical(:,2));
                Masses(n) = Nodes(n).Mass + sum([M1,M2],2);
            end
        end
        
        function [K] = calculateStiffness(Frames, Connectivity)
            % Calculate global stiffness matrix
            %
            % Example:
            %   K = stiff.struct.Framed.calculateStiffness(Frames, Connectivity);
            
            NNodes = max(Connectivity(:));
            NDoF = 3; % Number of global degrees of freedom in each node
            K = zeros(NDoF*NNodes);
            for n1 = 1 : NNodes
                IsConnected = Connectivity == n1;
                FrameIndex = find(any(IsConnected,2))'; % Row vector for looping
                IsConnected(FrameIndex,:) = not(IsConnected(FrameIndex,:)); % Indicates opposite nodes
                
                % Contributions of frames in which n1 corresponds to j-node                
                Index1 = (n1-1)*NDoF+1 : n1*NDoF;
                for f = find(IsConnected(:,1))'
                    rji = Frames(f).rji; % Transformation of local to global in j-node                    
                    Rji = Frames(f).Rji; % Transformation of global to local in j-node
                    kjji = Frames(f).kjji;
                    K(Index1,Index1) = K(Index1,Index1) + rji*kjji*Rji;                    
                end
                
                % Contributions of frames in which n1 corresponds to i-node
                for f = find(IsConnected(:,2))'
                    rij = Frames(f).rij; % Transformation of local to global in i-node
                    Rij = Frames(f).Rij; % Transformation of global to local in i-node
                    Rji = Frames(f).Rji; % Transformation of global to local in j-node
                    kiij = Frames(f).kiij;
                    kij = Frames(f).kij;
                    
                    
                    K(Index1,Index1) = K(Index1,Index1) + rij*kiij*Rij;
                    
                    n2 = Connectivity(f,IsConnected(f,:));
                    Index2 = (n2-1)*NDoF+1 : n2*NDoF;
                    K(Index1,Index2) = rij * kij * Rji;
                    K(Index2,Index1) = K(Index1,Index2)';
                end
            end
        end
        
        function [SubIndex] = getSubIndex(NodeIndex, Restraints)
            % [SubIndex] = stiff.struct.Framed.getSubIndex();
            %
            
            NNodes = numel(NodeIndex);
            NDoF = size(Restraints, 1);
            
            LoopIndex = 1 : NNodes; % Loop-index
            BasicIndex = 1 : NDoF; % Basic index for indexing and calculation
            AddIndex = NDoF * (LoopIndex-1); % Addition for loop-indexing
            AddSubIndex = NDoF * (NodeIndex-1); % Addition for loop-calculus
            
            SubIndex = zeros(1, NDoF*NNodes);
            for n = LoopIndex
                SubIndex(BasicIndex+AddIndex(n)) = BasicIndex + AddSubIndex(n);
            end
            if any(Restraints(:))
                SubIndex(Restraints(:)) = [];
            end
        end
        
        function plotStruct(Ax, Nodes, Frames, Connectivity)
            %
            % stiff.struct.Framed.plotStruct(Ax, Nodes, Frames, Connectivity)
            %
            % Input:
            %   Ax:
            %   Location: numeric. Matrix()
            
            
            % Plotting frames
            NFrames = numel(Frames);
            iLoc = [Nodes(Connectivity(:,1)).Location];
            jLoc = [Nodes(Connectivity(:,2)).Location];
            X = [iLoc(1,:); jLoc(1,:); nan(1,NFrames)];
            Y = [iLoc(2,:); jLoc(2,:); nan(1,NFrames)];
            Z = [iLoc(3,:); jLoc(3,:); nan(1,NFrames)];
            for n = 1 : numel(Ax)
                plot3(Ax(n), X(:), Y(:), Z(:));
            end
            
            % Plotting restraints
            Logical = any([Nodes(:).Restraints], 1);
            if any(Logical)
                Dist = 0.02 * max(abs([Frames(:).Length]));
                stiff.struct.Framed.plotSupports(Ax, Dist, Nodes(Logical));
            end
        end % function plotStruct(...)
        
        function plotSupports(Ax, Dist, Nodes)
            %
            % shmg7.stiff.struct.Framed.plotSupports(Ax, Dist, Restraints, Location)
                        
            FixedData =  [
                0       0       0
                -Dist   0       0
                -Dist   0       -Dist
                Dist    0       -Dist
                Dist    0       0
                0       0       0
                0       -Dist   0
                0       -Dist   -Dist
                0       Dist    -Dist
                0       Dist    0
                0       0       0
                ]';
            
            PinnedData = [
                0       0       0
                -Dist   0       -Dist
                Dist    0      -Dist
                0       0       0
                0       -Dist   -Dist
                0       Dist    -Dist
                0       0       0
                ]';
            
            
            Sum = sum([Nodes(:).Restraints], 1);
            IsFixed = Sum == 3;
            IsPinned = Sum == 2;
            FixedIndex = find(IsFixed);
            PinnedIndex = find(IsPinned);
            
            NCoordinates = 3;
            NFixed = numel(FixedIndex);
            NPinned = numel(PinnedIndex);
            
            NFixedPoints = size(FixedData, 2);
            NPinnedPoints = size(PinnedData, 2);
            NPoints = NFixed*NFixedPoints + NPinned*NPinnedPoints;
            Data = zeros(3, NPoints + NFixed + NPinned);
            ColCounter = 0;
            Location = [Nodes(:).Location];
            for n = 1 : NFixed
                NodeIndex = FixedIndex(n);
                DataIndex = (1:NFixedPoints) + ColCounter;
                for c = 1 : NCoordinates
                    Data(c,DataIndex) = FixedData(c,:) + Location(c,NodeIndex);
                end
                ColCounter = ColCounter + NFixedPoints + 1;
                Data(:,ColCounter) = nan;
            end
            
            for n = 1 : NPinned
                NodeIndex = PinnedIndex(n);
                DataIndex = (1:NPinnedPoints) + ColCounter;
                for c = 1 : NCoordinates
                    Data(c,DataIndex) = PinnedData(c,:) + Location(c,NodeIndex);
                end
                ColCounter = ColCounter + NPinnedPoints + 1;
                Data(:,ColCounter) = nan;
            end
            
            for n = 1 : numel(Ax)
                plot3(Ax(n), Data(1,:), Data(2,:), Data(3,:));
            end
        end % plotSupports(...)
    end % methods (Static)
end % classdef