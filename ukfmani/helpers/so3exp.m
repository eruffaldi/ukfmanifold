% Emanuele Ruffaldi 2017 @ SSSA
function R = so3exp(omega)

theta = norm(omega);
if theta < 1e-12
    % Original
    %     A = 1;
    %     B = 0.5;
    %     C = 1/6;
    %     S = zeros(3);
    %     R = eye(3) + A*S + B*S^2;
    %     V = eye(3) + B*S + C*S^2;
    
        N = 10;
        R = eye(3);
        xM = eye(3);
        cmPhi = skew(omega);
        for n = 1:N
            xM = xM * (cmPhi / n);
            R = R + xM;
        end
        
        % Project the resulting rotation matrix back onto SO(3)
        R = R / sqrtm( R'*R ) ;
    
else
    %Original
    A = sin(theta)/theta;
    B = (1-cos(theta))/(theta^2);
    C = (theta-sin(theta))/(theta^3);
    S = skew(omega);
    R = eye(3) + A*S + B*S^2;
    %V = eye(3) + B*S + C*S^2;
    
    %Barfoot
    if 0==1    
        axis = omega/theta;
        cp = cos(theta);
        sp = sin(theta);
        sa = skew(axis);
        
        R = cp*eye(3) + (1-cp)*axis*(axis') + sp*sa;
    end
        assert(abs(det(R)-1) < 1e-5,'unitary');
end
