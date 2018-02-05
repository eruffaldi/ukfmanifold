% SE3 in 4x4 matrix
%
% Emanuele Ruffaldi 2017 @ SSSA
function m = makeSE3Mat()

m1= makeRot();
m = [];
m.type = {'SE3Mat'};
m.inv = @se3inv;
m.prod = @(x,y) mflat(munflat(x)*munflat(y));
m.m1 = makeRot();
m.log = @se3log;
m.exp = @se3exp;
m.delta = @(x,y) se3log(munflat(x)*munflat(se3inv(y)));
m.pack = @(x) reshape(x{1},[],1);
m.unpack = @(x) {reshape(x,4,4)};
m.flat = @(x) reshape(x,1,[]);
m.unflat = @(x) reshape(x,4,4);

% Jacobian(y) * X
m.step = @(X,y) mflat(munflat(m.exp(y))*munflat(X));
m.group = 16;
m.alg = 6;
m.count = 1;

function y = se3log(x)

t = getpos(x);
R = getrot(x);

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
u = iV*t';
y = [omega(:);u(:)]';

function y = se3exp(x)

omega = getomega(x);
u = getvel(x);

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
y = mflat(y);



function u = mflat(x)

u = x(:);

function u = munflat(x)

u = reshape(x,4,4);

function y = se3inv(x)

R = getrot(x);
y = eye(4);
y(1:3,1:3) = R';
y(1:3,4) = -y(1:3,1:3)*getpos(x)';
y = mflat(y);

function u = ubuild2(rot,pos)

u = eye(4);
u(1:3,1:3) = reshape(rot,3,3);
u(1:3,4) = pos;

function R = getrot(x)

M = munflat(x);
R = M(1:3,1:3);

function p = getpos(x)

p = x(13:15);
p = p(:)';

function u = getomega(x)

u = x(1:3);

function u = getvel(x)

u = x(4:6);