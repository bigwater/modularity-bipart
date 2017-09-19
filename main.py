

import numpy as np
import networkx as nx



def calcModularityMatrix(adjMatrix, N, M, degs):
    modularityMat = np.zeros(shape=(N, N))
    
    for i in range(N):
    	for j in range(N):
    		modularityMat[i][j] = adjMatrix[i][j] - (0.0 + degs[i] * degs[j]) / (2 * M)

    return modularityMat


def isclose(a, b, rel_tol=1e-08, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def calcQ(s, degs, A, N, M):
	Q = 0
	verifyS = 0
	for i in range(N):
		for j in range(N):
			verifyS += A[i][j]

	print 'sum of G equals 2M ? = ', verifyS == 2*M
	print 'sum of deg equals 2M ? = ', sum(degs) == 2*M

	for i in range(N):
		assert A[i][i] == 0

	s = map(int, s)
	#print s
	assert len(s) == N

	for i in range(N):
		for j in range(N):
			Q = Q + ( A[i][j] - (0.0 + degs[i] * degs[j]) / (2 * M) ) * (s[i] * s[j])

	Q1 = 0
	for i in range(N):
		for j in range(N):
			Q1 = Q1 + ( A[i][j] - (0.0 + degs[i] * degs[j]) / (2 * M) ) * (s[i] * s[j] + 1)

	#print Q, ' ', Q1
	assert isclose(Q, Q1)
	#Q = float(Q) / (4 * M)
	return Q / (4.0 * M)


def calcQ1(s, modularityMat, N, M):
	st = np.matrix(s)
	#print st.shape, ' ', modularityMat.shape, ' ', st.T.shape
	Q1 = 1.0 / (4 * M) * st * modularityMat * st.T
	#print Q1.shape

	return Q1.item((0,0))
	#return 1.0 / (4 * M) * st * modularityMat * st.T


def calcQ2(modularityMat, s, betas, vecs, M):
    Q = 0
    st = np.matrix(s)
    for i in range(0, len(s)):
        bi = betas[i]
        ui = np.matrix(vecs[:,i])
        #print ui.shape, ' ', ui.T.shape
        #print 'aa', (ui * st.T).item(0,0) ** 2 * bi
        Q += (ui * st.T) ** 2 * bi
    
    Q = Q / (4 * M)
    return Q.item(0,0)


def calcWithNetX(edgepairs):
    nG = nx.Graph(edgepairs)
    B = nx.linalg.modularity_matrix(nG)
    #print 'B'
    #print B
    #print nx.modularity_spectrum(nG)

    return B


if __name__ == '__main__':

    N = 0

    pairs = []
    with open('graph.txt', 'r') as handle:

        for line in handle:
            if not line.strip():
                continue  # This skips blank lines

            values = map(int, line.split())
            if (max(values) > N):
                N = max(values)
            pairs.append(values)


    print 'N = ', N
    G = np.zeros(shape=(N, N))

    for p in pairs:
    	a = p[0] - 1
    	b = p[1] - 1
    	G[a][b] = 1
    	G[b][a] = 1

    degs = np.zeros(N)
    M = 0
    for i in range(N):
    	d1 = 0
        for j in range(N):
        	if (i != j and G[i][j] == 1):
        	    degs[i] += 1
        	    M += 1

    M = M / 2
    print 'M = ', M
    #print degs
    print 'sum of degrees = ', sum(degs)

    B1 = calcWithNetX(pairs)

    modularityMat = calcModularityMatrix(G, N, M, degs)
    print 'compare my modularityMat with networkx = ', (modularityMat == B1).all()

    t1 = np.ones(shape=(1,N))

    print 'ModularityMat is symmetric ? = ', (modularityMat.T == modularityMat).all()

    for i in range(N):
    	#print sum(modularityMat[i])
        assert isclose(sum(modularityMat[i]), 0.0, 1e-8, 1e-9)

    #print "largest eigenvalue beta1 = ", max(np.linalg.eigvals(modularityMat))
    eiValues, eiVecs = np.linalg.eigh(modularityMat)

    #print len(np.linalg.eig(modularityMat))
    print "largest eigenvalue beta1 = ", max(eiValues)

    u1ind = np.argmax(eiValues)

    b1 = eiValues[u1ind]
    u1 = eiVecs[:, u1ind]

    print 'b1 = ', b1
    print 'u1 = ', u1

    groupone = [1,2,3,4,5,6,7,8,11,12,13,14,17,18,20,22] #fact
   
    s = np.ones(N)
    for t in enumerate(u1):
    	if (t[1] > 0):
    		s[t[0]] = -1
    	elif isclose(t[1],0):
    		s[t[0]] = -1

    for i in range(N):
    	sp = 'circle'
    	if (i+1) in groupone:
    		sp = 'square'

    	if (s[i] == 1):
    		print "%d [shape=%s, style=filled, fillcolor=red]" % (i+1, sp)
    	else:
    		print "%d [shape=%s, style=filled, fillcolor=green]" % (i+1, sp)

    Q = calcQ (s, degs, G, N, M)
    Q1 = calcQ1 (s, modularityMat, N, M)
    Q2 = calcQ2(modularityMat, s, eiValues, eiVecs, M)

    print 'Q = ', Q
    print 'Q equals to matrix representation of Q ? = ', isclose(Q, Q1)
    print 'Q2 = ', Q2

















