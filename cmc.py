import numpy as np
import pandas as pd

#ns = [0.01, 0.02, 0.03]
#Q = np.array([[0, 0.001, 0], [0.001, 0, 0.001], [0, 0.001, 0]])

def Q2QQ(Q, ns):
  N = Q.shape[0]
  pops = [str(i) for i in range(N)]
  states = pops + [x+y for x in pops for y in pops]
  
  QQ = pd.DataFrame(np.zeros([N*N+N, N*N+N]))
  QQ.index = QQ.columns = states
  
  for pop in pops:
    QQ.loc[pop+pop, pop] = ns[int(pop)]
    QQ.loc[pop+pop, pop+pop] = - ns[int(pop)] - 2 * Q[int(pop)].sum()
    for pop2 in set(pops).difference(pop):
      QQ.loc[pop+pop, pop+pop2] = Q[int(pop), int(pop2)]
      QQ.loc[pop+pop, pop2+pop] = Q[int(pop), int(pop2)]
      QQ.loc[pop2+pop, pop+pop] = Q[int(pop2), int(pop)]
      QQ.loc[pop+pop2, pop+pop] = Q[int(pop2), int(pop)]
      QQ.loc[pop+pop2, pop+pop2] = - Q[int(pop)].sum() - Q[int(pop2)].sum()
      for pop3 in set(pops).difference([pop, pop2]):
        QQ.loc[pop+pop2, pop+pop3] = Q[int(pop2), int(pop3)]
        QQ.loc[pop+pop2, pop3+pop2] = Q[int(pop), int(pop3)]
    
    return QQ


#P = np.array([[0,0.2,0.8],[0,1,0],[0,0,1]])

def P2PP(P):
  N = P.shape[0]
  pops = [str(i) for i in range(N)]
  states = pops + [x+y for x in pops for y in pops]
  
  PP = pd.DataFrame(np.zeros([N*N+N, N*N+N]))
  PP.index = PP.columns = states
  
  for pop in pops:
    for pop2 in pops:
      for pop3 in pops:
        for pop4 in pops:
          PP.loc[pop+pop2, pop3+pop4] = P[int(pop), int(pop3)] * P[int(pop2), int(pop4)]
    
    return PP





