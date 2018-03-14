% SO3 as matrix
%
% Emanuele Ruffaldi 2017 @ SSSA
function m = makeRot()


m = [];
m.type = {'Rot'};
m.inv = @(X) flatten(unflatten(X)');
m.prod = @(X,Y) flatten(unflatten(X)*unflatten(Y));
m.log = @(X) so3log(unflatten(X));
m.exp = @(x) flatten(so3exp(x));
m.delta = @(X,Y) so3log(unflatten(X)*unflatten(Y)'); 
m.step = @(X,y) flatten(so3exp(y)*unflatten(X));
m.meancov = @manimeancov;
m.count = 1;
m.transport = @(X,t,Y) t;
m.group = 9; % as matrix
m.alg = 3;
m.pack = @(x) x(:)';
m.unpack = @(x) reshape(x,3,3);
m.islie = 1;
m.s = int_manisetup([],[],m);

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
y = x(:)';

function y = unflatten(x)
y = reshape(x,3,3);
