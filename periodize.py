def periodize(B):
    Bp = np.zeros((B.shape[0]*2-1, B.shape[1]*2-1))
    Bp[:B.shape[0],:B.shape[1]] = B
    Bp[:B.shape[0],B.shape[1]-1:] = B[:,::-1]
    Bp[B.shape[0]-1:,:B.shape[1]] = B[::-1,:]
    Bp[B.shape[0]-1:,B.shape[1]-1:] = B[:,::-1]
    return Bp
