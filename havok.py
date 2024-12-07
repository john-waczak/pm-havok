import numpy as np
from scipy.linalg import expm
from scipy.signal import lsim
from tqdm import tqdm

# rows is the embedding size
def build_hankel(z, rows):
    cols = len(z) - rows + 1
    H = np.empty((rows, cols))

    for i in range(rows):
        H[i,:] = z[i:i+cols]
    return H

def build_hankel_multi(Zs, rows):
    return np.hstack([build_hankel(z, rows) for z in Zs])

def true_polys(rows, dt, r, center=False): 
    m = rows // 2

    p = np.linspace(-m*dt, m*dt, rows)
    U = []
    for j in range(r):
        if (center):
            U.append(p ** (j + 1))
        else: 
            U.append(p ** j)
    U = np.vstack(U).T

    Q = np.empty((rows, r)) 
    
    # Perform Gram-Schmidt to yield ONB

    # 1. normalize u1 to obtain e1_hat
    # 2. e2 = u2 - (e1_hat * u2)u2 then normalize to obtain e2_hat
    # 3. ... repeat via projection for rest of vectors
    for j in range(r): 
        v = U[:, j]
        for k in range(j - 1): 
            # r_jk = Q[:, k].T @ U[:, j]
            r_jk = np.dot(Q[:, k], U[:, j])
            v -= (r_jk * Q[:, k])
        r_jj = np.linalg.norm(v)
        Q[:, j] = v / r_jj
    return Q


def HAVOK(H, dt, r): 
    U, s, Vh = np.linalg.svd(H, full_matrices=False)
    V = Vh.T

    # use polynomials to correct off diagonals
    # to enforce anti-symmetry
    polys = true_polys(H.shape[0], dt, r)
    for j in range(r): 
        # if our vector and the poly are anti-aligned
        # then re-align
        if (np.dot(U[:,j], polys[:,j]) < 0.0):
            U[:, j] *= -1
            V[:, j] *= -1

    # Split up embeddings for difference computation
    V1 = V[:-1, :r]
    V2 = V[1:, :r]

    # note: we need to use pinv since slices of V may
    # not be unitary
    A = (V2.T @ np.linalg.pinv(V1.T) - np.eye(r)) / dt
    return A, U[:, :r], s[:r], V[:, :r]


def sHAVOK(H, dt, r):
    # instead of splitting V after SVD
    # we split H before SVD to maintain
    # unitarity of V slices
    H1 = H[:, :-1]
    H2 = H[:, 1:]

    U1, s1, Vh1 = np.linalg.svd(H1, full_matrices=False)
    U2, s2, Vh2 = np.linalg.svd(H2, full_matrices=False)

    # create views of V w/o adjoint
    V1 = Vh1.T 
    V2 = Vh2.T

    polys = true_polys(H.shape[0], dt, r)
    for j in range(r):
        if (np.dot(U1[:,j], polys[:,j]) < 0.0):
            U1[:,j] *= -1
            V1[:,j] *= -1        
        if (np.dot(U2[:,j], polys[:,j]) < 0.0):
            U2[:,j] *= -1
            V2[:,j] *= -1

    # now we don't need pinv sinve V1 and V2 are unitary
    # interestingly here we are only truncating *after*
    # the multiplication... will that have an impact?
    A = ((V2.T @ V1)[:r,:r] - np.eye(r)) / dt


    return A, U1[:, :r], s1[:r], V1[:, :r]


def make_expM(A, B, dt, r_model, n_control):
    r = r_model + n_control
    M = np.zeros((r,r))
    M[:r_model, :r_model] = A * dt
    M[:r_model, r_model:r] = B * dt
    expM = expm(M)

    expA = expM[:r_model, :r_model]
    expB = expM[:r_model, r_model:r]

    return expA, expB

# def take_step(v_now, f_now, expA, expB):
#     return expA @ v_now + expB @ f_now


def eval_havok(z, ts, n_embedding, r_model, n_control):
    r = r_model + n_control

    # construct Hankel matrix
    #print("...building Hankel matrix")
    H = build_hankel(z, n_embedding)

    # get true values of time series
    # after forming Hankel matrix
    dt = ts[1] - ts[0]
    z_x = H[-1, :]
    t_x = np.arange(ts[n_embedding], ts[n_embedding] + dt*len(z_x),step=dt)

    # compute sHAVOK decomp
    # M, U, s, _ = sHAVOK(H, dt, r)
    #print("...fitting HAVOK model")
    M, U, s, _ = sHAVOK(H, dt, r)

    # using our H, U, and s, reconstruct the dembeddings
    # V = H.T @ (U * (1 / s)[np.newaxis, :])
    V = H.T @ U @ np.diag(1.0 / s)

    # pick out the A and B matrices
    A = M[:r_model, :r_model]
    B = M[:r_model, r_model:r]


    # pick out forcing values
    fvals = V[:, r_model:r]

    # construct exp matrices for integration
    expA, expB = make_expM(A, B, dt, r_model, n_control)

    # set up outgoing arrays
    Vout = np.zeros((V.shape[0], r_model))

    # set initial condition
    Vout[0,:] = V[0, :r_model]

    # time-evolve the system
    #print(f"...integrating from t={t_x[0]:.3f} to t={t_x[-1]:.3f}")
    for i in range(Vout.shape[0]-1):
        Vout[i+1,:] = expA @ Vout[i,:] + expB @ fvals[i,:]

    # reconstruct Hankel matrix and get time series
    H_pred = U[:, :r_model] @ np.diag(s[:r_model]) @ Vout.T    
    z_pred = H_pred[-1,:]

    return z_x, z_pred, t_x, U, s, Vout, A, B, fvals





def eval_havok_multi(zs, ts, n_embedding, r_model, n_control):
    r = r_model + n_control

    # construct Hankel matrix
    # print("...building Hankel matrices")
    Hs = []
    dts = []
    zs_x = []
    ts_x = []
    # for i in tqdm(range(len(zs))):
    for i in range(len(zs)):
        H = build_hankel(zs[i], n_embedding)
        Hs.append(H)

        dt = ts[i][1] - ts[i][0]
        dts.append(dt)
        
        z_x = H[-1,:]
        # t_x = np.arange(ts[i][n_embedding], ts[i][n_embedding] + dt*len(z_x), step=dt)
        t_x = ts[i][n_embedding-1:]

        zs_x.append(z_x)
        ts_x.append(t_x)

    H = np.hstack(Hs)        

    assert np.allclose(dts, dts[0])
    dt = dts[0]

    # compute sHAVOK decomp
    # print("...fitting HAVOK model")
    M, U, s, _ = sHAVOK(H, dt, r)

    # using our H, U, and s, reconstruct the dembeddings
    Vs = [H.T @ U @ np.diag(1.0 / s) for H in Hs]

    # pick out forcing values
    fvals = [V[:, r_model:r] for V in Vs]

          # set up predictded embedding matrices
    Vs_out = [np.zeros((V.shape[0], r_model)) for V in Vs]

    # set initial condition
    for i in range(len(Vs)):
        Vs_out[i][0,:] = Vs[i][0,:r_model]


    # pick out the A and B matrices
    A = M[:r_model, :r_model]
    B = M[:r_model, r_model:r]


    # construct exp matrices for integration
    expA, expB = make_expM(A, B, dt, r_model, n_control)


    # time-evolve the system
    zs_pred = []

    # for k in tqdm(range(len(ts_x))):
    for k in range(len(ts_x)):
        # perform integration
        for i in range(Vs_out[k].shape[0]-1):
            Vs_out[k][i+1,:] = expA @ Vs_out[k][i,:] + expB @ fvals[k][i,:]

        H_pred = U[:, :r_model] @ np.diag(s[:r_model]) @ Vs_out[k].T
        zs_pred.append(H_pred[-1,:])


    return zs_x, zs_pred, ts_x, U, s, Vs_out, A, B, fvals




# Assuming you already have a model, integrate the existing
# model for new time series points (t, z(t))
# this is useful for checking performance on a test set
def integrate_havok(z, t, n_embedding, r_model, n_control, A, B, U, s):
    r = r_model + n_control

    H = build_hankel(z, n_embedding)
    V = H.T @ U @ np.diag(1.0 / s)

    # pick out forcing values
    fvals = V[:, r_model:r]

    # compute dt
    dt = t[1] - t[0]
    z_x = H[-1, :]
    t_x = t[n_embedding-1:]


    # construct exp matrices for integration
    expA, expB = make_expM(A, B, dt, r_model, n_control)

    # set up outgoing arrays
    Vout = np.zeros((V.shape[0], r_model))

    # set initial condition
    Vout[0,:] = V[0, :r_model]

    # time-evolve the system
    # print(f"...integrating from t={t_x[0]:.3f} to t={t_x[-1]:.3f}")
    for i in range(Vout.shape[0]-1):
        Vout[i+1,:] = expA @ Vout[i,:] + expB @ fvals[i,:]

    # reconstruct Hankel matrix and get time series
    H_pred = U[:, :r_model] @ np.diag(s[:r_model]) @ Vout.T    
    z_pred = H_pred[-1,:]

    return z_x, z_pred, t_x

