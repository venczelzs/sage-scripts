# test.sage

def multi_indices(q):
    """
    Returns all the multiindices of multivariate polynomial of degree
    deg on n variables
    """
    n = len(q.variables())
    deg = q.degree()
    psfs = partitions_of_fixed_size(deg, n)
    #
    #
    #
    return psfs
            
def partitions_of_fixed_size(deg, n):
    ps = filter(lambda p: len(p) <= n, Partitions(deg))
    lps = []
    S3 = Permutations(3)
    for p in ps:
        lp = []
        k = len(p)
        for i in range(k):
            lp += [p[i]]
        if (k < n):
            lp += ([0] * (n - k))
        lps.append(lp)
    print lps[0]
    return lps

def permute_multi_index(i):
    return i
