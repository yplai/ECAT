import sys,random,copy
from collections import Counter

# a node is a collection of branches (or a leaf, with a label)
# each node has a branchlength

class Newick:
  def __init__(self,s,i=0):
    if s=="": return
    self.s,self.start = s,i
    #print s[i],i,l
    if s[i]=="(":
      self.isleaf,self.subtrees = False,[]
      i += 1
      #self.intlabel = l
      #l = l+1
      while True:
        node = Newick(s,i) 
        self.subtrees.append(node)
        j = node.next
        if s[j]==",": 
          i = j+1
        elif s[j]==")":
          self.branchlength,self.next = self.scan_branchlen(s,j+1)
          return
        else: print "error: unexpected format at",s[j]; sys.exit(0)
    else: # leaf node
      j = s.find(':',i)
      branchlen,next = self.scan_branchlen(s,j)
      self.isleaf,self.label,self.branchlength,self.next = True,s[i:j],branchlen,next

  # input: j = index of ":" char
  # returns: (float len,int next)
  def scan_branchlen(self,s,j):
    if j==len(self.s) or s[j]==';': return 0,j+1
    elif s[j]==":":
      j += 1; k = j
      while s[k] in "0123456789.": k += 1
      branchlen = float(s[j:k])
      return branchlen,k
    else: print "error: unexpected format2 at",s[j]; sys.exit(0)

  def indented_print(self,i=0):
    print ' '*i,
    if self.isleaf==True: print self.branchlength,self.label,self.nid
    else:
      print self.branchlength,self.nid
      #print self.branchlength
      for sub in self.subtrees: sub.indented_print(i+1)

  def get_strains(self):
    if self.isleaf==True: return [self.label]
    else:
      strains = []
      for sub in self.subtrees: strains += sub.get_strains()
      return strains

  def get_all_nodes(self):
    nodes = [self]
    if self.isleaf==False:
      for sub in self.subtrees: nodes += sub.get_all_nodes()
      #nodes.append(self)
    return nodes


  def label_all_nodes(self,i=1):
    #print len(self.get_all_nodes())
    #print len(self.get_strains())
    for node in self.get_all_nodes():
      node.nid = i
      i += 1
      #if node.isleaf==False:
      #  print node.isleaf,node.nid
      #if node.isleaf==True:
      #  print node.isleaf,node.nid,node.label


  def remove(self,node):
    if self.isleaf==True: return # if node is a leaf, should be removed from parent's list of subtrees
    for i in range(len(self.subtrees)): self.subtrees[i].remove(node)
    for i in range(len(self.subtrees)): # could leave it empty! (Which is bad) or with 1 substree (which is OK)
      if self.subtrees[i]==node: del self.subtrees[i]; return

  # oldnode is a subtree or leaf (detached)
  # newnode is place to attach it to
  # split the branch to newnode (i.e. as subtree of parent), which avoids 3-way splits and works for leafs too
  # I need the parent
  def insert(self,oldnode,newnode):
    if self.isleaf==True: return 
    for i in range(len(self.subtrees)): 
      if self.subtrees[i]==newnode: 
        new_br_pt = Newick("")
        new_br_pt.isleaf = False
        new_br_pt.next = oldnode.next
        new_br_pt.branchlength = 0.5*newnode.branchlength 
        new_br_pt.subtrees = [newnode,oldnode]
        newnode.branchlength -= new_br_pt.branchlength # so they sum up to original branchlen
        self.subtrees[i] = new_br_pt
      else: self.subtrees[i].insert(oldnode,newnode)

  def toString(self):
    if self.isleaf: return "%s:%s" % (self.label,self.branchlength)
    else:
      s,first = "",True
      for sub in self.subtrees:
        if first==True: first = False
        else: s += ","
        s += sub.toString()
      #return "(%s):%s" % (s,self.branchlength)
      return "(%s)%s:%s" % (s,self.nid,self.branchlength)

  def sankoff(self,h): # h is a hash table from strain names to nucs
    if self.isleaf==True:
      changes = {"A":999,"G":999,"C":999,"T":999} # infinity
      changes[h[self.label]] = 0
      return changes
    else:
      changes = {"A":0,"G":0,"C":0,"T":0}
      for t in self.subtrees:
        ch = t.sankoff(h)
        for a in "AGCT":
          vals = []
          for b in "AGCT": vals.append((1 if a!=b else 0)+ch[b])
          changes[a] += min(vals)
      return changes

  def val_sank(self,h):
     for node in self.get_all_nodes():
       val = node.sankoff(h)
       node.skval = val
       #print node.nid,node.skval,node.isleaf

  def root_leaf(self,path=[1],allpath=[]):
    if self.isleaf==True:
      path.append(self.nid)
      allpath.append(path[:-1])
      path.pop()
      return
    for t in self.subtrees:
      path.append(t.nid)
      t.root_leaf(path,allpath)
      path.pop()
    #print path
    #print allpath
    return allpath

  def h_parent(self,allpath):
     parent = {}
     for path in allpath:
       for i,p in enumerate(path):
         if i==0:
            parent[p] = 'root'
         else:
            parent[p] = path[i-1]
     return parent


###########################################################################

def read_newick(fname):
  s = ""
  for line in open(fname): s += line.rstrip()
  #return Newick(s)
  tree = Newick(s)
  tree.label_all_nodes()
  return tree

if __name__=="__main__":
  tree = read_newick(sys.argv[1])
  #tree2 = tree.label_all_nodes()
  strains = tree.get_strains()
  for s in strains: print "#",s
  nodes = tree.get_all_nodes()
  for n in nodes:
    if n.isleaf==True:
      print n.isleaf,n.nid,n.label
    else:
      print n.isleaf,n.nid

  tree.indented_print()
  print "input tree:"
  print tree.toString()

