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
            print("Error: unary node " + str(node), flush = True)
        if len(children) == 2:
            left = children[0]
            right = children[1]
            output = self.activation(self.W(self.traverse(tree, left) + self.traverse(tree, right)) + self.Q(torch.tensor([[time]])))
        return output
    
    def forward(self, tree):
        output = self.traverse(tree, tree.root)
        results = self.projection(output)
        return results
    
    def getLoss(self, tree, target):
        results = self.forward(tree)
        loss = torch.mean(torch.abs((target - results) / target))
        return results, loss


demo = single_demography(10000)
tree = msprime.sim_ancestry(demography = demo, samples = 10, ploidy = 1).first()
print(tree.draw_text())

model = gRNN(num_params = 1)
optimizer = torch.optim.SGD(model.parameters(), lr=0.01, momentum=0.9, dampening=0.0)
results, loss = model.getLoss(tree, torch.tensor([[10000]]))
optimizer.zero_grad()
loss.backward()
clip_grad_norm_(model.parameters(), 5, norm_type=2.)
optimizer.step()




for epoch in range(max_epochs):
  print("Epoch %d" % epoch)
  pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(trn)).start()
  for step, tree in enumerate(trn):
     predictions, loss = model.getLoss(tree)
     optimizer.zero_grad()
     loss.backward()
     clip_grad_norm(model.parameters(), 5, norm_type=2.)
     optimizer.step()
     pbar.update(step)
 
