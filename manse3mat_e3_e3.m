
		function m = customManifold()
			m = [];
			m.type = 'special';
			m.inv = @minvert;
			m.prod = @mproduct;
			m.log = @mlog;
				m.exp = @mexp;
			
			m.delta = @mdelta;
			m.step = @mstep;
			m.group = 22;
			m.alg = 12;
			m.count = 3;
			m.pack = @mpack;
			m.unpack = @munpack;

		end
		
		function z = minv(x1)
			z = zeros(22,1);
			z(1:16) = se3mat_inv(x1(1:16));
			z(17:19) = -x1(17:19);
			z(20:22) = -x1(20:22);
		end

		function z = mproduct(x1,x2)
			z = zeros(22,1);
			z(1:16) = reshape(reshape(x1(1:16),4,4)*reshape(x2(1:16),4,4),1,[]);
			z(17:19) = x1(17:19)+x2(17:19);
			z(20:22) = x1(20:22)+x2(20:22);

		end

			function z = mlog(x1)
				z = zeros(12,1);
				z(1:6) = se3mat_log(x1(1:16));
			z(7:9) = x1(17:19);
			z(10:12) = x1(20:22);

			end
	
			function z = mexp(x1)
				z = zeros(22,1);
				z(1:16) = se3mat_exp(x1(1:6));
			z(17:19) = x1(7:9);
			z(20:22) = x1(10:12);
		end
	
		function z = mdelta(x1,x2)
			z = zeros(12,1);
			z(1:6) = se3mat_log(reshape(x1(1:16),4,4)*reshape(se3mat_inv(x2(1:16)),4,4));
			z(7:9) = x1(17:19)-x2(17:19);
			z(10:12) = x1(20:22)-x2(20:22);
		end

		function z = mstep(x1,x2)
			z = zeros(22,1);
			z(1:16) = reshape(reshape(se3mat_exp(x2(1:6)),4,4)*reshape(x1(1:16),4,4),1,[]);
			z(17:19) = x1(17:19)+x2(7:9);
			z(20:22) = x1(20:22)+x2(10:12);
		end

		function z = mpack(x1)
			z = zeros(22,size(x1,2));
			for j=1:size(x1,2)
				z(1:16,j) = reshape(x1{1,j},1,[]);
				z(17:19,j) = x1{2,j};
				z(20:22,j) = x1{3,j};
			end
		end

		function z = munpack(x1)
			z = cell(3,size(x1,2));
			for j=1:size(x1,2)
				z{1,j} = reshape(x1(1:16,j),4,4);
				z{2,j} = x1(17:19,j);
				z{3,j} = x1(20:22,j);
			end
		end


function y = se3mat_inv(x)

	M = reshape(x,4,4);
	R = M(1:3,1:3);
	y = eye(4);
	y(1:3,1:3) = R';
	y(1:3,4) = -y(1:3,1:3)*M(1:3,4);
	y = reshape(y,1,[]);
end


function y = se3mat_log(x)
	
	M = reshape(x,4,4);
	t = M(1:3,4);
	R = M(1:3,1:3);

    theta = acos((trace(R)-1)/2); %acos(max(-1,min((trace(R)-1)/2,1)));
    if isreal(theta) == 0
        R = R/abs(det(R));
        theta =  acos((trace(R)-1)/2);
    end

    if abs(theta) < 1e-10
        B = 0.5;
        SO = (1/(2))*(R-R');  % =skew(omega)
        iV = eye(3); % use directly skew of omega
    else
        A = sin(theta)/theta;
        B = (1-cos(theta))/(theta*theta);
        SO = (1/(2*A))*(R-R');  % =skew(omega)
        %??
        % syms x real
        % A = sin(x)/x
        % B= (1-cos(x))/(x*x)
        % Q=1/(x*x)*(1 - A/2/B)
        % taylor(Q,x,0)
        %       x^4/30240 + x^2/720 + 1/12
        Q= 1/(theta^2)*(1 - A/2/B);
        iV = eye(3) - 1/2*SO + Q*SO^2; % use directly skew of omega
    end
    omega = [SO(3,2) SO(1,3) SO(2,1)];

	%y = [m1.log(getrot(x)), getpos(x)]; % not exact
	u = iV*t;
	y = [omega(:);u(:)]';
end

function y = se3mat_exp(x)

	omega = x(1:3);
	u = x(4:6); 

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
	        V = eye(3);
	        pxn = eye(3);
	        for n = 1:N
	            xM = xM * (cmPhi / n);
	            pxn = pxn*cmPhi/(n + 1);
	            R = R + xM;
	            V = V + pxn;
	        end
	        
	        % Project the resulting rotation matrix back onto SO(3)
	        R = R / sqrtm( R'*R ) ;
	    
	else
	    %Original
	    if 1==1
	        A = sin(theta)/theta;
	        B = (1-cos(theta))/(theta^2);
	        C = (theta-sin(theta))/(theta^3);
	        S = skew(omega);
	        R = eye(3) + A*S + B*S^2;
	        V = eye(3) + B*S + C*S^2;
	    else
	    %Barfoot
	        
	        axis = omega/theta;
	        cp = cos(theta);
	        sp = sin(theta);
	        cph = (1 - cos(theta))/theta;
	        sph = sin(theta)/theta;
	        sa = skew(axis);
	        
	        R = cp*eye(3) + (1-cp)*axis*(axis') + sp*sa;
	        V = sph * eye(3) + (1 - sph) * axis * (axis') + cph * sa;
	    end
	    
	end

	y = eye(4);
	y(1:3,1:3) = R;
	y(1:3,4) = V*u(:);
	y = reshape(y,1,[]);

	end

