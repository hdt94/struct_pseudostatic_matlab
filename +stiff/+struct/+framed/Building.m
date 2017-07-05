classdef Building < stiff.struct.Framed
    % stiff.struct.framed.Building
    
    methods
        function Obj = Building(Nodes, Frames, Connectivity)
            Obj@stiff.struct.Framed(Nodes, Frames, Connectivity);
        end
        
        function [K, M, ZLevels] = getReduced(Obj)
            NNodes = numel(Obj.Nodes);
            
            Restraints = [Obj.Nodes(:).Restraints];
            NodeIndex = find(not(all(Restraints,1)));
            SubIndex = Obj.getSubIndex(NodeIndex, Restraints(:,NodeIndex));
            
            Location = [Obj.Nodes(:).Location];
            ZLocation = Location(end,NodeIndex);
            ZLevels = unique(ZLocation);
            ZLevels(ZLevels <= 0) = [];                  
                        
            Obj.defineStiffness();
            Mass = Obj.getNodeMasses(Obj.Nodes, 1:NNodes, Obj.Frames, Obj.Connectivity);
            [K, M] = Obj.reduce(ZLevels, ZLocation, Obj.Stiffness(SubIndex,SubIndex), Mass(NodeIndex));
        end
        
        function [Displacement, Force, Shear] = solveResponseSpectrum(Obj, Spectrum)
            validateattributes(Spectrum, {'stiff.Spectrum'}, {'scalar'}, '', 'Spectrum');
            
            % Parameters
            [K, M, ZLevels] = Obj.getReduced();
            NNodes = size(K,1);            
            [Shapes, Values] = eig(K, M);
            Periods = 2*pi ./ sqrt(diag(Values));
            Sa = diag(Spectrum.getSa(Periods)); % Pseudo-acceleration
            L = (Shapes'*M*Shapes) \ diag(Shapes'*M*ones(NNodes,1)); % Modal participation factor
            
            % Modal results
            U = Shapes * L * (Values\Sa);                        
            F = M * Shapes * L * Sa;
            [~, Index] = sort(ZLevels, 'descend');
            V = cumsum(F(Index,:), 1);
            V = V(Index,:);
            
            % Modal combinations
            Force = struct('Absolute',sum(abs(F),2), 'SRSS',sqrt(sum(F.^2,2)));
            Displacement = struct('Absolute',sum(abs(U),2), 'SRSS',sqrt(sum(U.^2,2)));
            Shear = struct('Absolute',sum(abs(V),2), 'SRSS',sqrt(sum(V.^2,2)));
        end
    end % methods
    
    methods (Static)
        function [K, M] = reduce(ZLevels, ZLocation, Stiffness, Masses)
            NLevels = numel(ZLevels);
            K = Stiffness;
            M = zeros(NLevels);
            for n = 1 : NLevels
                Index = find(ZLocation == ZLevels(n));
                K(Index(1),:) = sum(K(Index,:),1);
                K(:,Index(1)) = sum(K(:,Index),2);
                M(n,n) = sum(Masses(Index));
                
                Index = Index(2:end);
                K(Index,:) = [];
                K(:,Index) = [];
                ZLocation(Index) = [];
                Masses(Index) = [];
            end
        end
        
        function Restraints = getNodeRestraints(ZLocation)
            NNodes = numel(ZLocation);
            Restraints = true(3, NNodes);
            IsUpperZero = ZLocation > 0;
            Restraints(1,IsUpperZero) = false; % Release X-global defree of freedom
        end
    end % methods (Static)
end % classdef