classdef NSR_10 < stiff.Spectrum   
    % Example:
    %   ScaleFactor = 1;
    %   Obj = stiff.spectrum.NSR_10(ScaleFactor);
    %   Obj.setAa(0.25);
    %   Obj.setAv(0.25);
    %   Obj.setFa(1.15);
    %   Obj.setFv(1.55);
    %   Obj.setI(1.5);
    %   Period = [0:0.1:1, 1.2, 1.5, 1.7, 2, 2.5, 3, 3.5, 4, 5, 8, 11, 15];
    %   Sa = Obj.getSa(Period);
    %   display([Period', Sa']);
    
    properties (SetAccess = private)
        % Numeric scalar (1,1).
        Aa = nan;
        % Numeric scalar (1,1).
        Av = nan;
        % Numeric scalar (1,1).
        Fa = nan;
        % Numeric scalar (1,1).
        Fv = nan;
        % Numeric scalar (1,1).
        I = nan;
    end
    
    properties (Dependent)
        % Numeric scalar (1,1).
        Tc
        % Numeric scalar (1,1).
        Tl
    end
    
    methods
        function Obj = NSR_10(ScaleFactor)
            Obj@stiff.Spectrum(ScaleFactor);
        end
        
        function [Sa] = getSa(Obj, Period)
            validateattributes(Period, {'numeric'}, {'vector','nonnegative'}, '', 'Period');
            Sa = Obj.ScaleFactor * Obj.I * Obj.calculateSa(Period, Obj.Tc, Obj.Tl, Obj.Aa, Obj.Fa, Obj.Av, Obj.Fv);
        end
        
        function Tc = get.Tc(Obj)
            Tc = 0.48 * Obj.Av * Obj.Fv / (Obj.Aa * Obj.Fa);
        end
        
        function Tl = get.Tl(Obj)
            Tl = 2.4 * Obj.Fv;
        end
        
        function setAa(Obj, Aa)
            validateattributes(Aa, {'numeric'}, {'scalar','positive'}, '', 'Aa');
            Obj.Aa = Aa;
        end
        
        function setAv(Obj, Av)
            validateattributes(Av, {'numeric'}, {'scalar','positive'}, '', 'Av');
            Obj.Av = Av;
        end
        
        function setFa(Obj, Fa)
            validateattributes(Fa, {'numeric'}, {'scalar','positive'}, '', 'Fa');
            Obj.Fa = Fa;
        end
        
        function setFv(Obj, Fv)
            validateattributes(Fv, {'numeric'}, {'scalar','positive'}, '', 'Fv');
            Obj.Fv = Fv;
        end
        
        function setI(Obj, I)
            validateattributes(I, {'numeric'}, {'scalar','positive'}, '', 'I');
            Obj.I = I;
        end
    end % methods
    
    methods (Static)
        function [Sa] = calculateSa(Period, Tc, Tl, aa, fa, av, fv)
            C1 = 2.5 * aa * fa;
            C2 = 1.2 * av * fv;
            C3 = 1.2 * av * fv * Tl;
            
            Sa = zeros(1, numel(Period));
            Logical = Period <= Tc;
            Sa(Logical) = C1;
            Logical = (Period > Tc) & (Period <= Tl);
            Sa(Logical) = C2 ./ Period(Logical);
            Logical = Period > Tl;
            Sa(Logical) = C3 ./ (Period(Logical).^2 + C2/Tl);
        end
    end % methods (Static)
end % classdef