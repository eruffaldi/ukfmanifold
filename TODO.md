

- test fusion
- test projective
- manifold of positive
- manifold for estimating covariances
- manifold SE3 right
- HARD Epipolar constraint between projected points in two perspective views, see Roberto Tron's page
- MAYVE Symmetric, positive definite matrices

    A point X on the manifold is represented as a symmetric positive definite
     matrix X (nxn). Tangent vectors are symmetric matrices of the same size
     (but not necessarily definite).

     The Riemannian metric is the bi-invariant metric, described notably in
     Chapter 6 of the 2007 book "Positive definite matrices"
     by Rajendra Bhatia, Princeton University Press.

- example spd vs wishart
- examples
- smoothing
    http://becs.aalto.fi/en/research/bayes/ekfukf/documentation.pdf
    eq. 3.61
    
    input: all the states
    output:
        use X(k, kalman) and f to computer X(k+1,pred)
        D = C(k+1 pred, k) / P(k+1,pred)
        P(k,smoothed) = P(k,kalman) D (P(k+1,smoothed)-P(k+1,pred)) D'
        mu(k,smoothed) = mu(k,smoothed) boxplus D (mu(k+1,smoothed) boxminus mu(k+1,pred))
    
    1) what about other variables/parameters? e.g. a common set of parameters? See slides about total variance and uncoditioning
        we aim at: (Xt|t|Xt+1|t = xt+1)
        but we don't know xt+1 so we use law of total exèectatopm and variance:
            E(X) = EZ( E(X|Y = Z) )
            Var(X) = EZ( Var(X|Y = Z) ) + VarZ( E(X|Y = Z) )
        so we obtain
            X(t|T)
            noting that Xt|t | Xt+1|t=Xt+1|T
        
    2) if we save both the X(k,kalman) and the X(k,pred) then we can skip the f
  
- svdsqrt reduction if too small