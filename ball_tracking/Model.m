classdef Model
    properties
        p1, % For each frame in point p,  p1.k, They are in the folo k1 < k2 < k3.
        p2,
        p3,
        support_set,
        len
    end
    
    methods
        
        function obj = Model(varargin)
            nb_arg = length(varargin);
            if nb_arg == 0
                return;
            end
            obj.p1 = varargin{1};
            obj.p2 = varargin{2};
            obj.p3 = varargin{3};
        end
        
        function [] = updateModel(obj,p1,p2,p3)
            obj.p1 = p1;
            obj.p2 = p2;
            obj.p3 = p3;
        end
        
        % v1::(vx,vy) , p1::(px,py,pk) , a::(ax,ay)
        function [v1] = compute_v1(obj,a)
            d_k21 = obj.p2.k - obj.p1.k;
            v1 = [((obj.p2.x - obj.p1.x) / d_k21) - (a(1)*d_k21*0.5) ... 
                  ((obj.p2.y - obj.p1.y) / d_k21) - (a(2)*d_k21*0.5)];
        end
        
        function [a] = compute_a(obj)
            d_k21 = obj.p2.k - obj.p1.k;
            d_k32 = obj.p3.k - obj.p2.k;
            a = [2* (d_k21*(obj.p3.x - obj.p2.x) - d_k32*(obj.p2.x-obj.p1.x)) / ((d_k21 + d_k32) * d_k32 * d_k21) ...
                 2* (d_k21*(obj.p3.y - obj.p2.y) - d_k32*(obj.p2.y-obj.p1.y)) / ((d_k21 + d_k32) * d_k32 * d_k21)];
        end
        
        % Hypothesis model which is re-computed, after every iteration.
        % For the fitted model, get the
        % position of the ball at frame k and i-N < k < i+N;

        function [p] = compute_k_prime(obj,v1,a,k)
            d_k = k - obj.p1.k;
            p = [obj.p1.x+ d_k *v1(1) + d_k * d_k *a(1)*0.5 ...
                 obj.p1.y+ d_k *v1(2) + d_k * d_k *a(2)*0.5];
        end
        
    end
end