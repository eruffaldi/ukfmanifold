% TODO: explain better the fact that angular velocity of the state is
% GLOBAL rather than LOCAL

% build input and output manifolds
mx = manisetup(makeProduct(makeSE3Mat(),makeRn(3),makeRn(3))); 
mxt = makeSE3Mat(); % helper
mz = manisetup(makeSE3Mat());

% initial state and noise definition
x0 = mx.step(mx.exp([0,0,0,  0,
    1,0,   0,0,0,   0,0,0]),[pi/2,0.2,0,  0,0,0,   0,0,0,   0,0,0]);
P0 = 0.5*eye(mx.alg);
Q = 0.01*eye(mx.alg); % process noi!se
R = 1e-3*eye(mz.alg); % measure noise
zobs = mz.exp([pi/2,0,0, 0,0,1]);

wsigmax = ut_mweights2(mx.group,mx.alg,0.5);
wsigmax.sqrt = @svdsqrt; 

% observation is identity
% process is the integral
dt = 0.1;

% integrate 
f_fx = @(Tk,wk,vk) deal(mxt.step(Tk,[wk,vk]),wk,vk);
h_fx = @(Tk,wk,vk) Tk;

tic
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
toc
figure(1)
plot(deltas(10:end,:))
figure(2)
plot(states(10:end,:))
