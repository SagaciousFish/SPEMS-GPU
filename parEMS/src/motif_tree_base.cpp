#include <iostream>
#include <string>
#include <bitset>
#include <vector>
#include "utils.h"
#include "motif_tree_base.hpp"

using namespace std;

MotifTreeBase::MotifTreeBase(uint32_t _max_depth, Motifs& _motifs, const char* _name) : 
  domain_size(domain.size()), root(TreeNodeFast::allocateNode()), max_depth(_max_depth), x(max_depth, ' '), 
  mask((1<<domain_size) - 1), name(_name), motifs(_motifs) 
{
  ;
}
  
MotifTreeBase::~MotifTreeBase() 
{
  ;
}

void MotifTreeBase::setDomain(const string& _domain) 
{ 
  domain = _domain; 
  domain_size= domain.size(); 
  mask = (1 << domain_size) - 1; 
}
  
void MotifTreeBase::traverse() 
{
  motifs.clear();
  (this)->traverseRecursive(root, 0);
}

void MotifTreeBase::traverseOut() 
{
  cout << "######## Traverse " << name << " #############" << endl;
  (this)->traverse();
  cout << "count=" << motifs.size() << endl;
  for (uint32_t i=0; i<motifs.size(); i++) 
  {
    cout << motifs[i] << endl;
  }
}

void MotifTreeBase::insert(const std::string& motif) 
{
  (this)->insertRecursive(root, motif, 0);
}

void MotifTreeBase::intersect(const MotifTreeBase* other) 
{
  (this)->intersectRecursive(root, other->root, 0);
}  