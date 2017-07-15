
# UKF over General Manifolds

We express first the UKF over a Manifold (SO3, SE3 or their combination) extending the existing works on the topic
(e.g. the original Quaternion UKF). In particular we acknowledge that

- state is in a Lie Group (or combination of there of): dimension n
- noise is expressed in the Algebra of the Lie Group: dimension m != n
- we express a multivariate Gaussian as N(mu, Sigma) as usual

We define the operator ** as the combination of the group variable with its algebra,
e.g. quaternion combined with an angular velocity. We also define the inverse in the group, and the log as the logarithm of a value.

First we define the sigma point construction and synthesis:

Chi_i = mu ** [ 0_alg +v1 ... +vm -v1...-vm ]  
mu = mean(Chi_i)
Chi_not_mu_i = log(Chi_i ** inv(mu))
C  = 1/(2m+) Sum (Chi_not_mu_i' Chi_not_mu_i)

The mean(Chi_i) is an iterative algorithm:

    mu_0 = x0
    v_{i,k} = log(x_i ** inv(mu_k))
    mu_{k+1} = exp(1/N Sum v_{i,k}) ** mu_k

Then we deal with the Kalman aspect in particular the correction:

x_i =  log(K log(z ** inv(z_{obs})) ** x_est_{i-1}


Case of SO3
- group: Rotation 3x3
- algenra: vector R3
- product: R1 R2
- inverse: R'
- logarithm/exponential: rodriguez and its inverse

Case of quaternion
- group: quaternion 4
- algenra: vector R3
- product: q1 q2
- inverse: q*
- logarithm/exponential: rodriguez and its inverse

# General Manifold
Instead of log/exp we define: 
    delta: group,group -> alg
    int: group,alg -> alg

Chi_i = mu ** [ 0_alg +v1 ... +vm -v1...-vm ]  
mu = mean(Chi_i)
Chi_not_mu_i = delta(Chi_i,mui);
C  = 1/(2m+) Sum (Chi_not_mu_i' Chi_not_mu_i)

The mean(Chi_i) is an iterative algorithm:

    mu_0 = x0
    v_{i,k} = delta(x_i,mu_k)
    mu_{k+1} = int(mu_k,1/N Sum v_{i,k})

Then we deal with the Kalman aspect in particular the correction:

x_i =  int(x_est_{i-1}, K delta(z, z_{obs}))

# General


# Additional

The reference paper on UKF and Quaternion is not using the algebra but directly quaternions

Paper: http://kodlab.seas.upenn.edu/uploads/Arun/UKFpaper.pdf (http://ieeexplore.ieee.org/document/1257247/)
Our Paper: https://dl.dropboxusercontent.com/u/146279/papers/2017_C_ETFADiStefano.pdf
Bonus for large dimensions: http://www.cs.unc.edu/~welch/kalman/media/pdf/ACC02-IEEE1357.PDF
Averaging: http://www.acsu.buffalo.edu/~johnc/uf_att.pdf
Italian: https://re.public.polimi.it/retrieve/handle/11311/961634/40508/1_Manuscript.pdf

