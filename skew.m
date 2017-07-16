% Emanuele Ruffaldi 2017 @ SSSA
function S = skew(v)
S = [  0   -v(3)  v(2)
    v(3)  0    -v(1)
    -v(2) v(1)   0];