% TODO: explain better the fact that angular velocity of the state is
% GLOBAL rather than LOCAL

% build input and output manifolds
mx = manisetup(makeProduct(makeSE3Mat(),makeRn(3),makeRn(3))); 
mz = manisetup(makeSE3Mat());

mxt = makeSE3Mat(); % helper without the need of setup

% initial state and noise definition
x0 = mx.step(mx.exp([0,0,0,  0,1,0,   0,0,0,   0,0,0]),[pi/2,0.2,0,  0,0,0,   0,0,0,   0,0,0]);
P0 = 0.5*eye(mx.alg);
Q = 0.01*eye(mx.alg); % process noi!se
R = 1e-3*eye(mz.alg); % measure noise

z0 = mz.exp([pi/2,0,0, 0,0.1,1]);
zobsval = zeros(200,16);

% costant se3 velocity
o0 = zeros(6,1);
o0(4) = 0.1;
o0(1) = 0.2;
zobsval(1,:) = z0;
lzobsval = zeros(size(zobsval,1),mz.alg);
lzobsval(1,:) = mz.log(z0);
for I=2:size(zobsval,1)
    zobsval(I,:) = mz.step(zobsval(I-1,:),o0);
    lzobsval(I,:) = mz.log(zobsval(I,:));
end
zobs = @(t) zobsval(t,:);

wsigmax = ut_mweights2(mx.group,mx.alg,0.5);
wsigmax.sqrt = @svdsqrt; 

% observation is identity
% process is the integral
dt = 0.1;


f_fx = @(Tk,wk,vk) deal(mxt.step(Tk,[wk,vk]),wk,vk); % Xk = (Tk,wk,vk)
h_fx = @(Tk,wk,vk) Tk;

tic
% loop
deltas = zeros(200,mz.alg);
states = zeros(size(deltas,1),mx.group);
lstates = zeros(size(deltas,1),mx.alg);
for L=1:size(deltas,1)
    states(L,:) = x0;
    lstates(L,:) = mx.log(x0);
    
    [xp,Pp] = manistatestep(mx,x0,P0,f_fx,Q,wsigmax);
    [zm,Czz,Cxz] = manievalh(mx,mz,xp,Pp,h_fx,wsigmax);
    
    % Kalman update with observation noise (additive)    
    Pvv = Czz + R;
    K = Cxz/Pvv;
    P0 = (eye(size(P0)) - K * Pvv * K') * P0;
    delta = mz.delta(zobs(L),zm);
    x0 = mx.step(xp,(K*delta')');
    deltas(L,:) = delta;
end
%%
toc
figure(1)
subplot(3,1,1);
plot(deltas(:,1:3))
title('Deltas observation-prediction');
xlabel('sample');
subplot(3,1,2);
plot(deltas(:,4:6))
title('Deltas observation-prediction');
xlabel('sample');
subplot(3,1,3);
plot(sum(deltas.^2,2));
title('Norm of error');
xlabel('sample');
figure(2)
plot(states(10:end,:))
title('All states as matrix');

%% Latency estimation
dd = zeros(6,1);
for J=1:6
    figure(2+J);
    dd(J) = finddelay(lstates(:,J),lzobsval(:,J));

    plot([lstates(:,J),lzobsval(:,J)]);
end
disp('delay')
dd
