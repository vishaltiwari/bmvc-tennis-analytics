classdef PathClass
    properties
        node_arr
    end
    
    methods
        function obj = PathClass(varargin)
            nb_arg = length(varargin);
            if nb_arg == 0
                return;
            end
            obj.node_arr = varargin{1};
        end
    end
end