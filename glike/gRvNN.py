import contextlib
import io
import sys


import random
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.utils import clip_grad_norm_

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
        self.activation = F.relu
        self.projection = nn.Linear(size, num_params, bias = True)
    
    def traverse(self, tree, node):
        input = torch.tensor([[tree.time(node)]])
        children = tree.children(node)
        if len(children) == 0:
            i = F.sigmoid(self.Wi(input))
            u = F.relu(self.Wu(input))
            c = i * u
            o = F.sigmoid(self.Wo(input))
            h = o * F.relu(c)
        if len(children) == 1:
            print("Error: unary node.", flush = True)
        if len(children) == 2:
            left = children[0]
            right = children[1]
            output = self.activation(self.W(self.traverse(tree, left) + self.traverse(tree, right)) + self.Q(torch.tensor([[time]])))
        return c, h
    
    def forward(self, trees):
        predictions = []
        for tree in trees:
          output = self.traverse(tree, tree.root)
          prediction = self.projection(output)
          predictions.append(prediction)
        results = torch.cat(predictions, 0)
        return results
    
    def getLoss(self, trees, target):
        results = self.forward(trees)
        loss = MAPE_loss(results, target)
        return results, loss
    



class gRNN(nn.Module):
    def __init__(self, size = 100, num_params = 5):
        super(gRNN, self).__init__()
        self.Q = nn.Linear(1, size, bias = True)
        self.W = nn.Linear(size, size, bias = True)
        self.activation = F.relu
        self.projection = nn.Linear(size, num_params, bias = True)
    
    def traverse(self, tree, node):
        time = tree.time(node)
        children = tree.children(node)
        if len(children) == 0:
            output = self.activation(self.Q(torch.tensor([[time]])))
        if len(children) == 1:
            print("Error: unary node.", flush = True)
        if len(children) == 2:
            left = children[0]
            right = children[1]
            output = self.activation(self.W(self.traverse(tree, left) + self.traverse(tree, right)) + self.Q(torch.tensor([[time]])))
        return output
    
    def forward(self, trees):
        predictions = []
        for tree in trees:
          output = self.traverse(tree, tree.root)
          prediction = self.projection(output)
          predictions.append(prediction)
        results = torch.cat(predictions, 0)
        return results
    
    def getLoss(self, trees, target):
        results = self.forward(trees)
        loss = MAPE_loss(results, target)
        return results, loss


model = gRNN(num_params = 5)
optimizer = torch.optim.SGD(model.parameters(), lr=0.01, momentum=0.9, dampening=0.0)

for epoch in range(100):
    trees = []
    params = []
    for i in range(1000):
        t = random.random()
        r = random.random()
        N_ab = random.random()
        N_a = random.random()
        N_b = random.random()
        demo = twoway_admixture_demography(t * 1000, r, N_ab * 20000, N_a * 20000, N_b * 20000)
        tree = msprime.sim_ancestry(demography = demo, samples = {"AB":20}, ploidy = 1).first()
        trees.append(tree.copy())
        params.append([t, r, N_ab, N_a, N_b])
    
    params = torch.tensor(params) * 10000
    results, loss = model.getLoss(trees, params)
    optimizer.zero_grad()
    loss.backward()
    
    save_stdout = sys.stdout
    sys.stdout = io.open('trash', 'w')
    clip_grad_norm_(model.parameters(), 5, norm_type=2.)
    sys.stdout = save_stdout
    
    optimizer.step()
    print(loss, flush = True)



model = gRNN(num_params = 1)
optimizer = torch.optim.SGD(model.parameters(), lr=0.01, momentum=0.9, dampening=0.0)

for epoch in range(100):
    trees = []
    params = []
    for i in range(1000):
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
 
