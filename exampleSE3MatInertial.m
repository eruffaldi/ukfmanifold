

mx = manisetup(makeProduct(makeSE3Mat(),makeRn(3),makeRn(3))); 
mz = manisetup(makeRn(3),makeRn(3));

mxt = makeSE3Mat(); % helper without the need of setup

% initial state and noise definition
x0 = mx.step(mx.exp([0,0,0,  0,1,0,   0,0,0,   0,0,0]),[pi/2,0.2,0,  0,0,0,   0,0,0,   0,0,0]);
P0 = 0.5*eye(mx.alg);
Q = 0.01*eye(mx.alg); % process noi!se
R = 1e-3*eye(mz.alg); % measure noise

% BUILD OBSERVATION
z0 = mz.exp([pi/2,0,0, 0,0.1,1]);
zobsval = zeros(200,16);
v0 = zeros(6,1);
v0(4) = 0.1;
v0(1) = 0.2;
zobsval(1,:) = z0;
lzobsval = zeros(size(zobsval,1),mz.alg);
lzobsval(1,:) = mz.log(z0);
for I=2:size(zobsval,1)
    v0(2) = sin(I/100);
    zobsval(I,:) = mz.step(zobsval(I-1,:),v0);
    lzobsval(I,:) = mz.log(zobsval(I,:));
end
zobs = @(t) zobsval(t,:);




wsigmax = ut_mweights2(mx.group,mx.alg,0.5);
wsigmax.sqrt = @svdsqrt; 

dt = 0.1;

II = eye(3); 
% dot omega = - inv(I) ( omega cross I omega)
% I that is local
f_fx = @(Tk,wk,vk) se3step_inertial_no_acc(Tk,wk,vk,II,dt); 
h_fx = []; % @(Tk,wk,vk) deal(wk,vk);

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
