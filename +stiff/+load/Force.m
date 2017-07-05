classdef Force < stiff.Load
    % stiff.load.Force()
    
    properties (SetAccess = private)
        % Numeric positive scalar. Relative distance (0-1 fraction) from i-node.
        Distance = 0.5;
    end
    
    methods
        function Obj = Force(Value, Distance)
            Obj@stiff.Load(Value);
            validateattributes(Distance, {'numeric'}, {'scalar','>=',0,'<=',1}, '', 'Distance');
            Obj.Distance = Distance;
        end
        
        function [Shear, Moment] = calculateDiagrams(Obj, XVector)
            Length = XVector(end);
            NPoints = numel(XVector);
            Shear = zeros(1, NPoints);
            Moment = zeros(1, NPoints);
            
            XPoint = Obj.Distance * Length; % Point where load is applied
            Logical = XVector >= XPoint;
            Shear(Logical) = Obj.Value;
            Moment(Logical) = Obj.Value * (XVector(Logical)-XPoint);
        end
        
        function [iNode, jNode] = calculateReactions(Obj, Length)
            Shear = Obj.Value * 0.5;
            
            D1 = Obj.Distance * Length;
            D2 = (1-Obj.Distance) * Length;
            L2 = Length ^ 2;
            iMoment = Obj.Value * D1^2 * D2 / L2;
            jMoment = - Obj.Value * D1 * D2^2 / L2;
            
            iNode = [-Shear, iMoment];
            jNode = [Shear, jMoment];
        end
    end % methods
end % classdef