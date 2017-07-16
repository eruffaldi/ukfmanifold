% TODO: explain better the fact that angular velocity of the state is
% GLOBAL rather than LOCAL
%
%http://it.mathworks.com/matlabcentral/fileexchange/1176-quaternion-toolbox
mysetup('quaternions');
% build input and output manifolds
mx = manisetup(makeCom(makeQuat(),makeRn(3))); % quat and vel 
mz = manisetup(makeQuat());

% initial state and noise definition
x0 = mx.step([0,0,0,1, 0,0,0],[pi/2,0.2,0, 0,0,0]);
P0 = 0.5*eye(mx.alg);
Q = 0.01*eye(mx.alg); % process noise
R = 1e-3*eye(mz.alg); % measure noise
zobs = qomega2q([pi/2,0,0]);

wsigmax = ut_mweights2(mx.group,mx.alg,0.5);

% observation is identity
% process is the integral
dt = 0.1;

% Process is the integral of the omega by time
% Note: in the paper Kraft: they state qk q_noise q_vel
f_fx = @(qk,ok) deal(qmult(qomega2q(dt*ok),qk),ok);
h_fx = @(qk,ok) qk;

% loop
deltas = zeros(200,mz.alg);
states = zeros(size(deltas,1),mx.group);
for L=1:size(deltas,1)
    states(L,:) = x0;
    
    [xp,Pp] = manistatestep(mx,x0,P0,f_fx,Q,wsigmax);
    [zm,Czz,Cxz] = manievalh(mx,mz,xp,Pp,h_fx,wsigmax);
    
    % Kalman update with observation noise (additive)    
    Pvv = Czz + R;
    K = Cxz/Pvv;
    P0 = (eye(size(P0)) - K * Pvv * K') * P0;
    delta = mz.delta(zobs,zm);
    x0 = mx.step(xp,(K*delta')');
    deltas(L,:) = delta;
end
figure(1)
plot(deltas(10:end,:))
figure(2)
plot(states(10:end,:))
