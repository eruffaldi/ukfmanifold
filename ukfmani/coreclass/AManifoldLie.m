classdef AManifoldLie < AManifold
    %AMANIFOLD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function Y = step(obj,X,v)
                Y = X*obj.exp(v);
        end
    
        function Y = delta(obj,X1,X2)
            Y = obj.log(X1*obj.inv(X2));
        end
        
        function b = islie(obj)
            b = true;
        end
    end
    
    methods(Abstract)
        v = log(obj,x);
        x = exp(obj,v);        
        z = prod(obj,x,y);
        y = inv(obj,x);        
    end
        
end