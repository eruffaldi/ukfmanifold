classdef AManifoldWrap < AManifold
    %AMANIFOLDWRAP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
       fx
    end
    
    methods
        function obj = AManifoldWrap(fx)
           obj.fx = fx;
           obj.G = fx.group;
           obj.a = fx.alg;
        end
        
        function b = islie(obj)
        end
        
        function y = step(obj,x,v)
        end
        function v = delta(obj,x,y)
        end
         
        % given [n,C] cell array emits [n,G] array
        function xvals = pack(obj,xcells)
        end
        % given [n,G] array emits [n,C] cell
        function xcells = unpack(obj,xvals)
        end        
    end
    
    methods(Static)
        function e = wrap(model)
            if hasattr(model,'log')
                e = AManifoldWrapLie(model);
            else
                e = AManifoldWrap(model);
            end
        end
    end
    
end

