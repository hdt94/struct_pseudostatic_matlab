classdef Section < handle
    
    properties (SetAccess = private)
        % Positive numeric scalar (1,1). Area.
        Area = nan;
        % Positive numeric scalar (1,1). Inertia.
        Inertia = nan;
        % Positive numeric scalar (1,1). Mutiplier to modify area.
        Modifier_Area = 1;
    end
    
    methods
        function Obj = Section(Area, Inertia, Modifier_Area)
            validateattributes(Area, {'numeric'}, {'scalar','positive'}, '', 'Area');
            validateattributes(Inertia, {'numeric'}, {'scalar','positive'}, '', 'Inertia');
            Obj.Area = Area;
            Obj.Inertia = Inertia;
            if nargin > 2
                validateattributes(Modifier_Area, {'numeric'}, {'scalar','nonnegative'}, '', 'Modifier_Area');
                Obj.Modifier_Area = Modifier_Area;
            end
        end % function Obj = Section(...)
    end % methods
end % classdef