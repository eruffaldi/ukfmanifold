clear all;
clc;
close all;
%% load and divide
load 'head_pose';
 
%addpath 'C:\Users\bigDaddy\Dropbox\project_ale\pose_estimator\ukfmanifold-master'
%addpath '/home/cattaneo/Dropbox/project_ale/pose_estimator/ukfmanifold-master'
addpath C:\Users\ing_Cattaneo\Dropbox\project_ale\pose_estimator\ukfmanifold-master;
addpath C:\Users\ing_Cattaneo\Documents\MATLAB\lib\quaternions-1.3\quaternions;
Tx=headpose(:,10);
Ty=headpose(:,11);
Tz=headpose(:,12);
Px=headpose(:,13);
Py=headpose(:,14);
Pz=headpose(:,15);
 
% build input and output manifolds
 
 
 
mx = manisetup(makeProduct(makeSE3Mat(),makeRn(3),makeRn(3))); 
mz = manisetup(makeProduct(makeSE3Mat()));
%% data
start=20;
stop=1540;
 
Rr= headpose(start:stop,10:12);
Tr= headpose(start:stop,13:15);
zobs=nan(size(Tr,1),(4+3));
missingdata = sum(isnan(headpose(start:stop,10:15)),2) > 0;
 
 
 
% plot data before ukf
figure,
subplot(3,1,1); plot(Tr(:,1)); grid on; legend('Tx'); title('x_{traslation}'); xlabel('sec'); ylabel('pixel');
subplot(3,1,2); plot(Tr(:,2)); grid on; legend('Ty'); title('y_{traslation}'); xlabel('sec'); ylabel('pixel');
subplot(3,1,3); plot(Tr(:,3)); grid on; legend('Tz'); title('z_{traslation}'); xlabel('sec'); ylabel('pixel');
 
figure,
subplot(3,1,1); plot(Rr(:,1)); grid on; legend('Rx'); title('x_{rot}'); xlabel('sec'); ylabel('pixel');
subplot(3,1,2); plot(Rr(:,2)); grid on; legend('Ry'); title('y_{rot}'); xlabel('sec'); ylabel('pixel');
subplot(3,1,3); plot(Rr(:,3)); grid on; legend('Rz'); title('z_{rot}'); xlabel('sec'); ylabel('pixel');
close all;
%% 
for i=1:stop-start
    zobs(i,:)=[eul2quat(Rr(i,:)),Tr(i,:)];
end
 
%% UKF pose
 
N=stop;
 
% initial state and noise definition rotation and lineaqr traslation
initialguess=eul2rotm(Rr(start,:));
x0 = mx.step(mx.exp([initialguess(1,:), initialguess(2,:),initialguess(3,:),   0,0,0]),[Rr(start,:), Tr(start,:),  0,0,0,   0,0,0]);
 
P0 = 0.5*eye(mx.alg);
Q = 0.01*eye(mx.alg); % process noi!se
R = 1e-3*eye(mz.alg); % measure noise
 
wsigmax = ut_mweights2(mx.group,mx.alg,0.5);
wsigmax.sqrt = @svdsqrt; 
 
% observation is identity
% process is the integral
dt = 0.1;
 
 
f_fx = @(Tk,wk,vk) deal(mxt.step(Tk,[wk,vk]),wk,vk); % Xk = (Tk,wk,vk)
h_fx = @(Tk,wk,vk) Tk;
 
tic
% loop
deltas = zeros(N,mz.alg);
states = zeros(size(deltas,1),mx.group);
lstates = zeros(size(deltas,1),mx.alg);
for L=1:size(deltas,1)
    states(L,:) = x0;
    lstates(L,:) = mx.log(x0);
     
    
    [xp,Pp] = manistatestep(mx,x0,P0,f_fx,Q,wsigmax);
    [zm,Czz,Cxz] = manievalh(mx,mz,xp,Pp,h_fx,wsigmax);
     
    if missingdata(L) == 1
        % Kalman update with observation noise (additive)    
        Pvv = Czz + R;
        K = Cxz/Pvv;
        P0 = (eye(size(P0)) - K * Pvv * K') * P0;
        delta = mz.delta(zobs(L),zm);
        x0 = mx.step(xp,(K*delta')');
        deltas(L,:) = delta;
    end
     
     
     
end
%% plotting graph and data converter
 
Rr1=zeros(size(states,1),3);
for i=1:size(states,1)
    Rr1(i,:)=quat2eul(states(i,1:4));
end
 
figure('name','measures');
subplot(6,1,1); plot(deltas(:,1)); grid on;legend('Rx'); xlabel('sec'); ylabel('pixel');
subplot(6,1,2); plot(deltas(:,2)); grid on;legend('Ry'); xlabel('sec'); ylabel('pixel');
subplot(6,1,3); plot(deltas(:,3)); grid on;legend('Rz'); xlabel('sec'); ylabel('pixel');
subplot(6,1,4); plot(deltas(:,4)); grid on;legend('x'); xlabel('sec'); ylabel('pixel');
subplot(6,1,5); plot(deltas(:,5)); grid on;legend('y'); xlabel('sec'); ylabel('pixel');
subplot(6,1,6); plot(deltas(:,6)); grid on;legend('z'); xlabel('sec'); ylabel('pixel');
 
 
figure('name','stima_rotation');
subplot(3,1,1); plot(Rr1(:,1)); hold on; plot(Rr(start:stop,1)); grid on; xlabel('sec'); ylabel('pixel'); legend('sR_x','R_x');
subplot(3,1,2); plot(Rr1(:,2)); hold on; plot(Rr(start:stop,2)); grid on; xlabel('sec'); ylabel('pixel'); legend('sR_y','R_y');
subplot(3,1,3); plot(Rr1(:,3)); hold on; plot(Rr(start:stop,3)); grid on; xlabel('sec'); ylabel('pixel'); legend('sR_z','R_z');
 
figure('name','stima_traslation');
subplot(3,1,1); plot(Tr(start:stop,1),'-'); hold on; plot(states(:,8),'--'); grid on; xlabel('sec'); ylabel('pixel'); legend('x','s_x');
subplot(3,1,2); plot(Tr(start:stop,2),'-'); hold on; plot(states(:,9),'--'); grid on; xlabel('sec'); ylabel('pixel'); legend('y','s_y');
subplot(3,1,3); plot(Tr(start:stop,3),'-'); hold on; plot(states(:,10),'--'); grid on; xlabel('sec'); ylabel('pixel'); legend('z','s_z');