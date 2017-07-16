% reverse of qdecomp
% Emanuele Ruffaldi 2017 @ SSSA
function q = qomega2q(omega)

momega = norm(omega);
if momega == 0
    q = [0,0,0,1];
else
    domega = omega/momega;

    q = [domega*sin(momega/2),cos(momega/2)];
end