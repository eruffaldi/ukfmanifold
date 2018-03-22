classdef AManifoldProject < AManifold
    %AMANIFOLDPROJECT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n
        other
    end
    
    methods
        function obj = AManifoldProject(otherManifold,T,k)
            obj.other = otherManifold;
            obj.n = n
            obj.k = k;
            obj.T = T;
        end
        
        function outputArg = method1(obj,inputArg)
        end
    end
end

