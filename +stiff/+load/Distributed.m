classdef Distributed < stiff.Load
    % stiff.load.Distributed(Value)
    
    methods
        function Obj = Distributed(Value)
            Obj@stiff.Load(Value);
        end
        
        function [iNode, jNode] = calculateReactions(Obj, Length)
            Shear = Obj.Value * Length * 0.5;
            Moment = Obj.Value * Length^2  / 12;
            iNode = [-Shear, Moment];
            jNode = [Shear, -Moment];
        end
        
        function [Shear, Moment] = calculateDiagrams(Obj, XVector)
            Shear = Obj.Value * XVector;
            Moment = 0.5 * Obj.Value * (XVector.^2);
        end
    end % methods
end % classdef