classdef Struct < handle
    
    properties (SetAccess = private)
        % stiff.Node vector (1,NNodes). Structure nodes.
        Nodes = stiff.Node.empty;
    end
    
    methods
        function Obj = Struct(Nodes)
            validateattributes(Nodes, {'stiff.Node'}, {'row'}, '', 'Nodes');
            Obj.Nodes = Nodes;
        end
    end
end