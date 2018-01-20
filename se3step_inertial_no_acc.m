function [Tk,wk,vk] = se3step_inertial_no_acc(Tk0,wk0,vk0,II,dt)

IIL = inv(Tk0)*II;
vk = vk0;
wk = wk0 - inv(IIL)*(cross(wk0,IIL*wk0));
Tk = se3step(Tk0,[wk,vk]);