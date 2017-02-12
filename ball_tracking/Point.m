classdef Point
    properties
        x,
        y,
        k
    end
    methods
        function obj = Point(varargin)
            nVarargs = length(varargin);
            if nVarargs == 0
                obj.x=0;
                obj.y=0;
                obj.k=0;
                return;
            end
            obj.x = varargin{1};
            obj.y = varargin{2};
            obj.k = varargin{3};
        end
        
        function tf = eq(obj1,obj2)
            if obj1.x == obj2.x && obj1.y == obj2.y && obj1.k == obj2.k
                tf = true;
            else
                tf = false;
            end
        end
        
    end
end