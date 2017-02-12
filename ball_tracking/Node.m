classdef Node
    properties
        model;
        frame;
    end
    
    methods
        
        function obj = Node(varargin)
            nb_arg = length(varargin);
            if nb_arg == 0
                return;
            end
            obj.model = varargin{1};
            obj.frame = varargin{2};
        end
    end
    
end