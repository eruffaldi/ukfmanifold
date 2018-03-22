classdef AManifold
    %AMANIFOLD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        G
        a
        count
    end
    
    methods (Abstract)        
        y = step(obj,x,v);
        v = delta(obj,x,y);        
        % given [n,C] cell array emits [n,G] array
        xvals = pack(obj,xcells);
        % given [n,G] array emits [n,C] cell
        xcells = unpack(obj,xvals);
        
        b = islie(obj);
    end
end

