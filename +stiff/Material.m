classdef Material < handle
        
    properties (SetAccess = private)
        % Nonnegative numeric scalar (1,1). Material's density.
        Density = nan;        
        % Nonnegative numeric scalar (1,1). Elasticity modulus.
        Young = nan;
    end
    
    methods        
        function setDensity(Obj, Density)
            validateattributes(Density, {'numeric'}, {'vector','nonnegative'}, 'Density');
            NMaterials = numel(Obj);
            if numel(Density) ~= NMaterials
                error('Material:setDensity:invalidNumEl', 'Numel of Density must be equal to NMaterials');
            end
            for n = 1 : NMaterials
                Obj(n).Density = Density(n);
            end
        end
        
        function setYoung(Obj, Young)
            validateattributes(Young, {'numeric'}, {'vector','nonnegative'}, 'Young');
            NMaterials = numel(Obj);
            if numel(Young) ~= NMaterials
                error('Material:setYoung:invalidNumEl', 'Numel of Young must be equal to NMaterials');
            end
            for n = 1 : NMaterials
                Obj(n).Young = Young(n);
            end
        end
    end % methods
end % classdef