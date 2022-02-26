import contextlib
import io
import sys


import random
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.utils import clip_grad_norm_


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
    return np.array(ns + ts + zs)




class gNN(nn.Module):
    def __init__(self, N):
        super(gNN, self).__init__()
        size = 6 * N - 3
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(size, 64),
            nn.ReLU(),
            nn.Linear(64, 4),
            nn.ReLU(),
            nn.Linear(4, 1),
        )
    
    def forward(self, trees):
        outputs = []
        for tree in trees:
            input = torch.tensor([BFSencode(tree)]).float()
            output = self.linear_relu_stack(input)
            outputs.append(output)
        predictions = torch.sigmoid(torch.cat(outputs, 0))
        return predictions
    
    def getLoss(self, trees, targets):
        predictions = self.forward(trees)
        criterion = nn.BCELoss()
        loss = criterion(predictions, targets)
        return predictions, loss



model = gNN(N = 100)
optimizer = torch.optim.SGD(model.parameters(), lr=0.01, momentum=0.9, dampening=0.0)

for epoch in range(100):
    trees = []
    targets = []
    for i in range(100):
        target = float(random.randint(0, 1))
        if target == 1.0:
            demo = twoway_admixture_demography(t = 10, r = 0.5, N_ab = 2000, N_a = 10000, N_b = 10000, m = 1e-6)
        else:
            demo = twoway_admixture_demography(t = 10, r = 0.8, N_ab = 2000, N_a = 10000, N_b = 10000, m = 1e-6)
        tree = msprime.sim_ancestry(demography = demo, samples = {"AB":100}, ploidy = 1).first()
        trees.append(tree.copy())
        targets.append([target])
    
    targets = torch.tensor(targets)
    predictions, loss = model.getLoss(trees, targets)
    tags = torch.round(predictions)
    accuracy = (tags == targets).float().mean()
    
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()
    
    print(accuracy, flush = True)





model = gNN(N = 100)
optimizer = torch.optim.SGD(model.parameters(), lr = 0.01)

for epoch in range(100):
    trees = []
    targets = []
    for i in range(100):
        target = float(random.randint(0, 1))
        if target == 1.0:
            demo = single_demography(10000)
        else:
            demo = single_demography(20000)
        tree = msprime.sim_ancestry(demography = demo, samples = {"A":100}, ploidy = 1).first()
        trees.append(tree.copy())
        targets.append([target])
    
    targets = torch.tensor(targets)
    predictions, loss = model.getLoss(trees, targets)
    tags = torch.round(predictions)
    accuracy = (tags == targets).float().mean()
    
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()
    
    print(accuracy, flush = True)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



def MAPE_loss(prediction, target):
    loss = torch.mean(torch.abs((prediction - target) / target))
    return loss

class gRNN(nn.Module):
    def __init__(self, size = 100, num_params = 5):
        super(gRNN, self).__init__()
        self.Wf = nn.Linear(1, size, bias = True)
        self.Wi = nn.Linear(1, size, bias = True)
        self.Wu = nn.Linear(1, size, bias = True)
        self.Wo = nn.Linear(1, size, bias = True)
        self.Uf = nn.Linear(size, size, bias = True)
        self.Ui = nn.Linear(size, size, bias = True)
        self.Uu = nn.Linear(size, size, bias = True)
        self.Uo = nn.Linear(size, size, bias = True)
        self.projection = nn.Linear(size, num_params, bias = True)
    
    # f-forget, i-input, u-update, o-output, c-state, h-hidden
    def traverse(self, tree, node):
        input = torch.tensor([[tree.time(node)]])
        children = tree.children(node)
        if len(children) == 0:
            i = torch.sigmoid(self.Wi(input))
            u = F.relu(self.Wu(input))
            c = i * u
            o = torch.sigmoid(self.Wo(input))
            h = o * F.relu(c)
        if len(children) == 1:
            print("Error: unary node.", flush = True)
        if len(children) == 2:
            c1, h1 = self.traverse(tree, children[0])
            c2, h2 = self.traverse(tree, children[1])
            f1 = torch.sigmoid(self.Uf(h1))
            f2 = torch.sigmoid(self.Uf(h2))
            i = torch.sigmoid(self.Wi(input) + self.Ui(h1) + self.Ui(h2))
            u = torch.relu(self.Wu(input) + self.Uu(h1) + self.Uu(h2))
            c = i * u + f1 * c1 + f2 * c2
            o = torch.sigmoid(self.Wo(input) + self.Uo(h1) + self.Uo(h2))
            h = o * F.relu(c)
        return c, h
    
    def forward(self, trees):
        predictions = []
        for tree in trees:
          c, h = self.traverse(tree, tree.root)
          prediction = self.projection(h)
          predictions.append(prediction)
        results = torch.cat(predictions, 0)
        return results
    
    def getLoss(self, trees, target):
        results = self.forward(trees)
        loss = MAPE_loss(results, target)
        return results, loss
    









model = gRNN(num_params = 1)
optimizer = torch.optim.SGD(model.parameters(), lr=0.01, momentum=0.9, dampening=0.0)

for epoch in range(100):
    trees = []
    params = []
    for i in range(100):
        N = random.random() * 20000
        demo = single_demography(N)
        tree = msprime.sim_ancestry(demography = demo, samples = {"A":20}, ploidy = 1).first()
        trees.append(tree.copy())
        params.append([N])
    
    params = torch.tensor(params)
    results, loss = model.getLoss(trees, params)
    optimizer.zero_grad()
    loss.backward()
    
    save_stdout = sys.stdout
    sys.stdout = io.open('trash', 'w')
    clip_grad_norm_(model.parameters(), 5, norm_type=2.)
    sys.stdout = save_stdout
    
    optimizer.step()
    print(loss, flush = True)
 
