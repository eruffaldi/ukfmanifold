classdef AManifoldWrapLie < AManifoldWrap
    %AMANIFOLDWRAP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = AManifoldWrapLie(fx)
            obj@AManifoldWrap(fx);
        end
        
        function b = islie(obj)
            b = true;
        end
        
        function v = log(obj,x)
        end
        
        function x = exp(obj,v)
        end
        
        function z = prod(obj,x,y)
        end
        
        function y = inv(obj,x)
        end
        
    end
end

