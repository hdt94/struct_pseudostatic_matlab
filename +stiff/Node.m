classdef Node < handle
    
    properties
        % Numeric vector (3,1). Loads considering global loads Px, Pz, My.
        Loads = zeros(3,1);
        % Numeric vector (3,1). Location considering global location X, Y, Z.
        Location = zeros(3,1);
        % Numeric scalar (1,1). Mass assigned to node.
        Mass = 0;
        % Logical vector (3,1). Restraints considering global degrees of freedom DX, DZ, GY.
        Restraints = false(3,1);
    end
    
    methods (Static)
        function Obj = create(NNodes, Location, Restraints, Loads)
            Obj(NNodes) = stiff.Node;
            Obj(:).setLocation(Location);
            Obj(:).setRestraints(Restraints);
            Obj(:).setLoads(Loads);
        end
    end
    
    methods
        function setLoads(Obj, Loads)
            validateattributes(Loads, {'numeric'}, {'nonnan'}, '', 'Loads');
            NNodes = numel(Obj);
            [Rows, Cols] =  size(Loads);
            if Rows ~= 3
                error('Node:setLoads:invalidRows', 'Loads must be NRows equal to 3');
            elseif Cols ~= NNodes
                error('Node:setLoads:invalidCols', 'Loads must be NCols equal to NNodes');
            end
            for n = 1 : NNodes
                Obj(n).Loads = Loads(:,n);
            end
        end
        
        function setLocation(Obj, Location)
            validateattributes(Location, {'numeric'}, {'nonnan'}, '', 'Location');
            NNodes = numel(Obj);
            [Rows, Cols] =  size(Location);
            if Cols ~= NNodes
                error('Node:setLocation:invalidCols', 'Location must be NCols equal to NNodes');
            elseif Rows == 2
                Location = [Location(1,:); zeros(1,NNodes); Location(2,:)];
            elseif Rows ~= 3
                error('Node:setLocation:invalidRows', 'Location must be NRows equal to 3');
            else
            end
            for n = 1 : NNodes
                Obj(n).Location = Location(:,n);
            end
        end
        
        function setMass(Obj, Mass)
            validateattributes(Mass, {'numeric'}, {'vector','>=',0}, '', 'Mass');
            NNodes = numel(Obj);
            if numel(Mass) ~= NNodes
                error('Node:setMass:invalidNumEl', 'Mass must be length equal to NNodes');
            end
            for n = 1 : NNodes
                Obj(n).Mass = Mass(n);
            end
        end
        
        function setRestraints(Obj, Restraints)
            validateattributes(Restraints, {'logical'}, {}, '', 'Restraints');
            NNodes = numel(Obj);
            [Rows, Cols] =  size(Restraints);
            if Rows ~= 3
                error('Node:setRestraints:invalidRows', 'Restraints must be NRows equal to 3');
            elseif Cols ~= NNodes
                error('Node:setRestraints:invalidCols', 'Restraints must be NCols equal to NNodes');
            end
            for n = 1 : NNodes
                Obj(n).Restraints = Restraints(:,n);
            end
        end
    end % methods
end % classdef