classdef Constant < handle
    % Class to storage values related to numeric constant.
    %
    % Example:
    %   Constant = stiff.Definition.Constant;
    %   Constant.setGravity(9.80665);
    
    properties (SetAccess = private)
        Gravity = 0;
    end
    
    methods
        function setGravity(Obj, Gravity)
            validateattributes(Gravity, {'numeric'}, {'scalar','nonnegative'}, '', 'Gravity');
            Obj.Gravity = Gravity;
        end
    end
end