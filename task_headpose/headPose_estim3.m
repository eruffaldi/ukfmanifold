clear all;
close all;
%% load and divide
load headPose;
 
%addpath 'C:\Users\bigDaddy\Dropbox\project_ale\pose_estimator\ukfmanifold-master'
%addpath '/home/cattaneo/Dropbox/project_ale/pose_estimator/ukfmanifold-master'
% addpath C:\Users\ing_Cattaneo\Dropbox\project_ale\pose_estimator\ukfmanifold-master;
% addpath C:\Users\ing_Cattaneo\Documents\MATLAB\lib\quaternions-1.3\quaternions;
Tx=pose_Tx1;
Ty=pose_Ty1;
Tz=pose_Tz1;
Px=pose_Rx1;
Py=pose_Ry1;
Pz=pose_Rz1;
 
% build input and output manifolds
 
 
 
mx = manisetup(makeProduct(makeSE3Mat(),makeRn(3),makeRn(3))); 
mz = manisetup(makeProduct(makeSE3Mat()));
 
mxt = makeSE3Mat(); % helper without the need of setup
%% data
start=20;
stop=1540;
 
Rr= [Px(start:stop,:),Py(start:stop,:),Pz(start:stop,:)];
Tr= [Tx(start:stop,:),Ty(start:stop,:),Tz(start:stop,:)];
t = (1:length(Rr))*(1/30);

%ZYX
 
figure,
subplot(3,1,1); plot(t,Tr(:,1)); grid on;  title('x_{traslation}'); xlabel('sec'); ylabel('mm');
subplot(3,1,2); plot(t,Tr(:,2)); grid on;  title('y_{traslation}'); xlabel('sec'); ylabel('mm');
subplot(3,1,3); plot(t,Tr(:,3)); grid on; title('z_{traslation}'); xlabel('sec'); ylabel('mm');
  
figure,
subplot(3,1,1); plot(t,Rr(:,1)); grid on; title('Z_{rot}'); xlabel('sec'); ylabel('angles (rad)');
subplot(3,1,2); plot(t,Rr(:,2)); grid on;  title('Y_{rot}'); xlabel('sec'); ylabel('angles (rad)');
subplot(3,1,3); plot(t,Rr(:,3)); grid on; title('X_{rot}'); xlabel('sec'); ylabel('angles (rad)');
  
zobs=zeros(size(Tr,1),16);
val=zeros(4,4);
%% 
% ZYX = roll yaw pitch terna sinistrosra
for i=1:stop-start
    val=[eul2rotm(Rr(i,:)),Tr(i,:)';[0,0,0,0]]; 
    zobs(i,:)=mz.exp([val(1,:),val(2,:),val(3,:),val(4,:)]);
end
 
%% UKF pose
 
N=stop-start;
 
% initial state and noise definition rotation and lineaqr traslation
initialguess=eul2rotm(Rr(start,:));
x0 = mx.step(mx.exp([1,0,0, 0,1,0,0,0,1,  0,0,0]),[Rr(start,:), Tr(start,:),  0,0,0,   0,0,0]);
 
% P0 = blkdiag(1e-6*eye(mx.alg-3),1e-6*eye(3));
% Q = blkdiag(1e-2*eye(mx.alg-3),1e-4,1e-4,1e-1); % process noise
% R = blkdiag(1e-3*eye(mz.alg-3),1e1,1e-2,1e-1); % measure noise
 
P0 = blkdiag(0.5*eye(mx.alg-3),0.5*eye(3));
Q = blkdiag(0.01*eye(mx.alg-3),1e-2,1e-2,1e-2); % process noise
R = blkdiag(1e-4*eye(mz.alg-3),1e-1,1e-1,1e-1); % measure noise
 
wsigmax = ut_mweights2(mx.group,mx.alg,0.5);
wsigmax.sqrt = @svdsqrt; 
 
% observation is identity
% process is the integral
dt = 1/30*10;
 
 
f_fx = @(Tk,wk,vk) deal(mxt.step(Tk,[wk,vk]),wk,vk); % Xk = (Tk,wk,vk)
h_fx = @(Tk,wk,vk) Tk;
 
tic
% loop
deltas = zeros(N,mz.alg);
states = zeros(size(deltas,1),mx.group);
lstates = zeros(size(deltas,1),mx.alg);
for L=1:N
    states(L,:) = x0;
    lstates(L,:) = mx.log(x0);
     
    
    [xp,Pp] = manistatestep(mx,x0,P0,f_fx,Q,wsigmax);
%     if  sum(isnan(xp),2) >0
%         xp=zeros(size(xp,1),size(xp,2));
%     end
%     
   [zm,Czz,Cxz] = manievalh(mx,mz,xp,Pp,h_fx,wsigmax);
     
    %if missingdata(L) == 1 
        % Kalman update with observation noise (additive)    
        Pvv = Czz + R;
        K = Cxz/Pvv;
        P0 = (eye(size(P0)) - K * Pvv * K') * P0;
        delta = mz.delta(zobs(L,:),zm);
        x0 = mx.step(xp,(K*delta')');
        deltas(L,:) = delta;
    %end
     
     
     
end
%% plotting graph and data converter from 4x4 to Rotation and traslation
toc
Rr1=zeros(size(states,1),3);
Tr1=zeros(size(states,1),3);
for i=1:size(states,1)
    T=[states(i,1:3);states(i,4:6);states(i,7:9)];
    Rr1(i,:)=quat2eul(rotm2quat(T));
    Tr1(i,:)=[states(i,4),states(i,8),states(i,12)];
end
%% 
figure('name','measures');
subplot(6,1,1); plot(t(1,1:end-1),deltas(:,1)); grid on; xlabel('sec'); ylabel('rad');
subplot(6,1,2); plot(t(1,1:end-1),deltas(:,2)); grid on; xlabel('sec'); ylabel('rad');
subplot(6,1,3); plot(t(1,1:end-1),deltas(:,3)); grid on; xlabel('sec'); ylabel('rad');
subplot(6,1,4); plot(t(1,1:end-1),deltas(:,4)); grid on; xlabel('sec'); ylabel('mm');
subplot(6,1,5); plot(t(1,1:end-1),deltas(:,5)); grid on; xlabel('sec'); ylabel('mm');
subplot(6,1,6); plot(t(1,1:end-1),deltas(:,6)); grid on; xlabel('sec'); ylabel('mm');
 
 
figure('name','stima_rotation');
subplot(3,1,1); plot(t(1,1:end-1),Rr1(:,1)); hold on; plot(t(1,1:end-1),Rr(1:end-1,1)); grid on; xlabel('sec'); ylabel('angles (rad)'); legend('sR_z','R_z');
subplot(3,1,2); plot(t(1,1:end-1),Rr1(:,2)); hold on; plot(t(1,1:end-1),Rr(1:end-1,2)); grid on; xlabel('sec'); ylabel('angles (rad)'); legend('sR_y','R_y');
subplot(3,1,3); plot(t(1,1:end-1),Rr1(:,3)); hold on; plot(t(1,1:end-1),Rr(1:end-1,3)); grid on; xlabel('sec');  ylabel('angles (rad)'); legend('sR_x','R_x');
 
figure('name','stima_traslation');
subplot(3,1,1); plot(t(1,1:end-1),Tr(1:end-1,1),'-'); hold on; plot(t(1,1:end-1),Tr1(:,1),'--'); grid on; xlabel('sec'); ylabel('mm'); legend('x','s_x');
subplot(3,1,2); plot(t(1,1:end-1),Tr(1:end-1,2),'-'); hold on; plot(t(1,1:end-1),Tr1(:,2),'--'); grid on; xlabel('sec'); ylabel('mm'); legend('y','s_y');
subplot(3,1,3); plot(t(1,1:end-1),Tr(1:end-1,3),'-'); hold on; plot(t(1,1:end-1),Tr1(:,3),'--'); grid on; xlabel('sec'); ylabel('mm'); legend('z','s_z');