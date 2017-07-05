classdef Rectangular < stiff.Section
    % Rectangular section.
    %
    % Example:
    %   Obj = stiff.section.Rectangular(Height, Width, Modifier_Area);
    
    properties (SetAccess = private)
        % Numeric positive scalar (1,1).
        Height = nan;
        % Numeric positive scalar (1,1).
        Width = nan;
    end
    
    methods
        function Obj = Rectangular(Height, Width, varargin)
            validateattributes(Height, {'numeric'}, {'scalar','positive'}, '', 'Height');
            validateattributes(Width, {'numeric'}, {'scalar','positive'}, '', 'Width');
            Area = Height * Width;
            Inertia = Height^3 * Width / 12;
            Obj@stiff.Section(Area, Inertia, varargin{:});
            Obj.Height = Height;
            Obj.Width = Width;
        end % function Obj = Rectangular(...)
    end % methods
end % classdef