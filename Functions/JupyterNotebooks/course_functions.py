import numpy as np
import scipy.linalg
import control
control.use_numpy_matrix(False)


def dlqr(A, B, Q, R):
    """Solve the discrete time lqr controller.

    x[k+1] = A x[k] + B u[k]

    cost = sum x[k].T*Q*x[k] + u[k].T*R*u[k]
    """

    # ref Bertsekas, p.151

    # first, try to solve the ricatti equation
    X = scipy.linalg.solve_discrete_are(A, B, Q, R)

    # compute the LQR gain
    K = scipy.linalg.inv(B.T @ X @ B + R) @ (B.T @ X @ A)

    eigVals = scipy.linalg.eig(A - B @ K)[0]

    return K, X, eigVals


def kalman(A, G, C, QN, RN, NN=None, continuous=True):
    
    A, G, C = np.array(A, ndmin=2), np.array(G, ndmin=2), np.array(C, ndmin=2)
    QN, RN =  np.array(QN, ndmin=2), np.array(RN, ndmin=2)
    if continuous:
        P, E, LT = control.care(A.T, C.T, np.dot(np.dot(G, QN), G.T), RN)
        MX = None
    else:
        P, E, LT = control.dare(A.T, C.T, np.dot(np.dot(G, QN), G.T), RN)
        MX = P@C.T/(C@P@C.T+RN)

    return control.statesp._ssmatrix(LT.T), control.statesp._ssmatrix(P), control.statesp._ssmatrix(E.real), control.statesp._ssmatrix(MX)
