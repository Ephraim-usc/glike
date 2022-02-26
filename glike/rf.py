
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn import metrics

import random
import tskit
import numpy as np

def BFSencode(tree):
    ns = []
    ts = []
    queue = [tree.root]
    
    while queue:
        node = queue.pop(0)
        ns.append(len(list(tree.samples(node))))
        ts.append(tree.time(node))
        children = tree.children(node)
        
        if len(children) == 2:
            n1 = len(list(tree.samples(children[0])))
            n2 = len(list(tree.samples(children[1])))
            if n1 >= n2:
                queue.append(children[0])
                queue.append(children[1])
            else:
                queue.append(children[1])
                queue.append(children[0])
    
    zs = [int(n == 1) for n in ns]
    return np.array(ns + ts)


def BFRencode(tree):
    ns = []
    ts = []
    queue = [tree.root]
    
    while len(ns) <= 510: # 1+2+...+256
        node = queue.pop(0)
        
        if node == -1:
            ns.append(0)
            ts.append(0)
            queue.append(-1)
            queue.append(-1)
            continue
        
        ns.append(len(list(tree.samples(node))))
        ts.append(tree.time(node))
        children = tree.children(node)
        
        if len(children) == 0:
            queue.append(-1)
            queue.append(-1)
            continue
        
        if len(children) == 2:
            n1 = len(list(tree.samples(children[0])))
            n2 = len(list(tree.samples(children[1])))
            if n1 >= n2:
                queue.append(children[0])
                queue.append(children[1])
            else:
                queue.append(children[1])
                queue.append(children[0])
    
    zs = [int(n == 1) for n in ns]
    return np.array(ns + ts)


def data_generate(size = 100):
    X = []
    Y = []
    for i in range(size):
        y = float(random.randint(0, 1))
        if y == 1.0:
            demo = twoway_admixture_demography(t = 10, r = 0.5, N_ab = 1000, N_a = 10000, N_b = 10000, m = 1e-6)
        else:
            demo = twoway_admixture_demography(t = 10, r = 0.7, N_ab = 1000, N_a = 10000, N_b = 10000, m = 1e-6)
        tree = msprime.sim_ancestry(demography = demo, samples = {"AB":1000}, ploidy = 1).first()
        X.append(BFSencode(tree))
        Y.append(y)
    
    X = np.array(X)
    Y = np.array(Y)
    return X, Y


model = RandomForestClassifier(n_estimators=100)

X, Y = data_generate(1000)
model.fit(X, Y)

for epoch in range(100):
  X, Y = data_generate(1000)
  Y_pred = model.predict(X)
  print("Accuracy:",metrics.accuracy_score(Y_pred, Y))
  model.fit(X, Y)


model = RandomForestClassifier(warm_start = True, n_estimators = 1)

for i in range(100):
    X, Y = data_generate(1000)
    if i == 0:
        model.fit(X, Y)
    Y_pred = model.predict(X)
    print("Accuracy:",metrics.accuracy_score(Y_pred, Y), flush = True)
    model.n_estimators += 1
    model.fit(X, Y)
    

    
    
    
def data_generate(size = 100):
    X = []
    Y = []
    for i in range(size):
        y = random.random() * 20000
        demo = twoway_admixture_demography(t = 10, r = 0.7, N_ab = 1000, N_a = y, N_b = 10000, m = 1e-6)
        tree = msprime.sim_ancestry(demography = demo, samples = {"AB":1000}, ploidy = 1).first()
        X.append(BFSencode(tree))
        Y.append(y)
    
    X = np.array(X)
    Y = np.array(Y)
    return X, Y

def data_generate(size = 100):
    X = []
    Y = []
    for i in range(size):
        y = random.random()
        demo = twoway_admixture_demography(t = 10, r = y, N_ab = 1000, N_a = 20000, N_b = 5000, m = 1e-6)
        tree = msprime.sim_ancestry(demography = demo, samples = {"AB":1000}, ploidy = 1).first()
        X.append(BFSencode(tree))
        Y.append(y)
    
    X = np.array(X)
    Y = np.array(Y)
    return X, Y

def data_generate(size = 100):
    X = []
    Y = []
    for i in range(size):
        y = random.random() * 0.5
        demo = twoway_admixture_demography(t = 10, r = 0.5 + y, N_ab = 1000, N_a = 10000, N_b = 10000, m = 1e-4)
        tree = msprime.sim_ancestry(demography = demo, samples = {"AB":1000}, ploidy = 1).first()
        X.append(BFSencode(tree))
        Y.append(y)
    
    X = np.array(X)
    Y = np.array(Y)
    return X, Y


model = RandomForestRegressor(warm_start = True, n_estimators = 10)

for i in range(100):
    X, Y = data_generate(1000)
    if i == 0:
        model.fit(X, Y)
    else:
        Y_pred = model.predict(X)
        print("MAE:", metrics.mean_absolute_error(Y_pred, Y), flush = True)
        model.n_estimators += 10
        model.fit(X, Y)

    
    
    
    
    
    
