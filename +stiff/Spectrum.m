classdef Spectrum < handle
    
    properties (SetAccess = private)
        % Numeric scalar (1,1). Multiplier of spectrum.
        ScaleFactor = nan;
    end
    
    methods
        function Obj = Spectrum(ScaleFactor)
            validateattributes(ScaleFactor, {'numeric'}, {'scalar'}, '', 'ScaleFactor');
            Obj.ScaleFactor = ScaleFactor;
        end
    end
    
    methods (Abstract)
        [Sa] = getSa(Obj, Period);
    end
end