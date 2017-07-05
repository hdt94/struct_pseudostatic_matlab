classdef Load < matlab.mixin.Heterogeneous
    
    properties (SetAccess = private)
        % Numeric scalar (1,1). Load value.
        Value = nan;
    end
    
    methods
        function Obj = Load(Value)
            validateattributes(Value, {'numeric'}, {'scalar'}, '', 'Value');
            Obj.Value = Value;
        end
    end
    
    methods (Abstract)
        [Shear1, Moment1, Shear2, Moment2] = calculateReactions(Obj);
    end
    
    methods (Sealed)
        function [Shear, Moment] = getDiagrams(Obj, XVector, iLoads)
            NLoads = numel(Obj);
            if NLoads == 1
                [Shear, Moment] = Obj.calculateDiagrams(XVector);
            else
                NPoints = numel(XVector);
                Sh = zeros(NLoads, NPoints);
                Mo = zeros(NLoads, NPoints);
                for n = 1 : NLoads
                    [Sh(n,:), Mo(n,:)] = Obj(n).calculateDiagrams(XVector);
                end
                Shear = sum(Sh, 1);
                Moment = sum(Mo, 1);
            end
            Shear = Shear + iLoads(2);
            Moment = Moment + iLoads(2)*XVector + iLoads(3);
        end
    end % methods (Sealed)
end % classdef