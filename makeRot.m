% SO3 as matrix
function m = makeRot()


m = [];
m.inv = @(X) flatten(unflatten(X))';
m.prod = @(X,Y) flatten(unflatten(X)*unflatten(Y));
m.log = @(X) so3log(unflatten(X));
m.exp = @(x) flatten(so3exp(x));
m.delta = @(X,Y) so3log(unflatten(X)*unflatten(Y)'); 
m.step = @(X,y) flatten(unflatten(X)*so3exp(y));
m.meancov = @manimeancov;
m.count = 1;
m.group = 9; % as matrix
m.alg = 3;
m.pack = @(x) x(:);
m.unpack = @(x) reshape(x,3,3);

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
    %A = sin(theta)/theta;
    %B = (1-cos(theta))/(theta^2);
    %C = (theta-sin(theta))/(theta^3);
    %S = skew(omega);
    %R = eye(3) + A*S + B*S^2;
    %V = eye(3) + B*S + C*S^2;
    
    %Barfoot
        
        axis = omega/theta;
        cp = cos(theta);
        sp = sin(theta);
        sa = skew(axis);
        
        R = cp*eye(3) + (1-cp)*axis*(axis') + sp*sa;
       
    
end

y = eye(4);
y(1:3,1:3) = R;

function S = skew(v)
S = [  0   -v(3)  v(2)
    v(3)  0    -v(1)
    -v(2) v(1)   0];

function omega = so3log(R)

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

    
function y = flatten(x)
y = x(:);

function y = unflatten(x)
y = reshape(x,3,3);
