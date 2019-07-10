#References:
#[1] 'Finding Prinicpal Paths in Data Space', M.J.Ferrarotti, W.Rocchia, S.Decherchi
#[2] 'Design and HPC Implementation of Unsupervised Kernel Methods in the Context of Molecular Dynamics', M.J.Ferrarotti, PhD Thesis.
#[3] https://github.com/mjf-89/PrincipalPath/blob/master/principalpath.py

import numpy
from scipy.spatial import distance

def initMedoidsPY(X, n, init_type, exclude_ids=[]):
    """
    Initialize NC medoids with init_type rational.

    Args:
        [ndarray float] X: data matrix

        [int] n: number of medoids to be selected

        [string] init_type: rational to be used
            'uniform': randomly selected with uniform distribution
            'kpp': k-means++ algorithm

        [ndarray int] exclude_ids: blacklisted ids that shouldn't be selected

    Returns:
        [ndarray int] med_ids: indices of the medoids selected
    """

    N=X.shape[0]
    D=X.shape[1]
    n=int(n)
    exclude_ids = numpy.array(exclude_ids, dtype=int)

    med_ids=-1*numpy.ones(n,int)
    numpy.random.seed(123)

    if(init_type=='uniform'):
        while(n>0):
            med_id = numpy.random.randint(0,N)
            if(numpy.count_nonzero(med_ids==med_id)==0 and numpy.count_nonzero(exclude_ids==med_id)==0):
                med_ids[n-1]=med_id
                n = n-1

    elif(init_type=='kpp'):
        accepted = False
        while(not accepted):
            med_id = numpy.random.randint(0,N)
            if(numpy.count_nonzero(exclude_ids==med_id)==0):
                accepted = True
        med_ids[0]=med_id

        for i in range(1,n):
            Xmed_dst = distance.cdist(X,numpy.vstack([X[med_ids[0:i],:],X[exclude_ids,:]]),'sqeuclidean')
            D2 = Xmed_dst.min(1)
            D2_n = 1.0/numpy.sum(D2)
            accepted = False
            while(not accepted):
                med_id = numpy.random.randint(0,N)
                if(numpy.random.rand()<D2[med_id]*D2_n):
                    accepted = True
            med_ids[i]=med_id
    else:
        raise ValueError('init_type not recognized.')

    return(med_ids)

def rkmPY(X, init_W, s):
    """
    Regularized K-means for principal path, MINIMIZER.

    Args:
        [ndarray float] X: data matrix

        [ndarray float] init_W: initial waypoints matrix

        [float] s: regularization parameter

        [matplotlib.axis.Axes] plot_ax: Axes for the 2D plot (first 2 dim of X), None to avoid plotting

    Returns:
        [ndarray float] W: final waypoints matrix

        [ndarray int] labels: final

    References:
        [1] 'Finding Prinicpal Paths in Data Space', M.J.Ferrarotti, W.Rocchia, S.Decherchi, [submitted]
        [2] 'Design and HPC Implementation of Unsupervised Kernel Methods in the Context of Molecular Dynamics', M.J.Ferrarotti, PhD Thesis.
    """

    #extract useful info from args
    N = X.shape[0]
    d = X.shape[1]
    NC = init_W.shape[0]-2

    #construct boundary matrix
    boundary = init_W[[0,NC+1],:]
    B=numpy.zeros([NC,d],float)
    B[[0,NC-1],:]=boundary

    #construct regularizer hessian
    AW = numpy.diag(numpy.ones(NC))+numpy.diag(-0.5*numpy.ones(NC-1),1)+numpy.diag(-0.5*numpy.ones(NC-1),-1)

    #compute initial labels
    XW_dst = distance.cdist(X,init_W,'sqeuclidean')
    u = XW_dst.argmin(1)

    #iterate the minimizer
    converged = False
    it = 0
    while(not converged):
        it = it+1
        #print('iteration '+repr(it))

        #compute cardinality
        W_card=numpy.zeros(NC+2,int)
        for i in range(NC+2):
            W_card[i] = numpy.sum(u==i)

        #compute centroid matrix
        C = numpy.ndarray([NC,d],float)
        for i in range(NC):
            C[i,:] = numpy.sum(X[u==i+1,:],0)

        #construct k-means hessian
        AX = numpy.diag(W_card[1:NC+1])

        #update waypoints
        W = numpy.matmul(numpy.linalg.pinv(AX+s*AW),C+0.5*s*B)
        W = numpy.vstack([boundary[0,:],W,boundary[1,:]])

        #compute new labels
        XW_dst = distance.cdist(X,W,'sqeuclidean')
        u_new = XW_dst.argmin(1)

        #check for convergence
        converged = not numpy.sum(u_new!=u)
        u=u_new

    return W
