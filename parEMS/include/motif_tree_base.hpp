#pragma once __MOTIF_TREE_BASE__

#include <iostream>
#include <string>
#include "utils.h"

#define POOL_SIZE 1024

struct TreeNodeFast
{
  union
  {
    struct
    {
      TreeNodeFast* children[4]; // TODO make dynamic
      size_t sharing_info;
    };

    TreeNodeFast* next;
  };

  static TreeNodeFast* free_head;

  static TreeNodeFast* allocateNode()
  {
    if (!free_head)
    {
      free_head = new TreeNodeFast[POOL_SIZE];
      TreeNodeFast* head = free_head;
      for (uint32_t i = 1; i < POOL_SIZE; i++)
      {
        head->next = head + 1;
        head++;
      }
      head->next = 0;
    }
    TreeNodeFast* head = free_head;
    free_head = head->next;
    memset((char*)head, 0, sizeof(TreeNodeFast));
    return head;
  }

  static void deallocateNode(TreeNodeFast* node)
  {
    node->next = free_head;
    free_head = node;
  }

  void print()
  {
    std::cout << "Node " << this << ", sharing info=" << sharing_info << std::endl;
    // for (size_t j=0; j<domain_size; j++) 
    // {
    //   std::cout << " child " << j << " ptr=" << children[j] << std::endl;
    // }
  }
};

class MotifTreeBase
{
protected:
  string domain;
  size_t domain_size;
  TreeNodeFast* root;

  uint32_t max_depth;
  std::string x;
  uint32_t mask;
  std::string name;
  Motifs& motifs;

public:
  virtual void traverseRecursive(TreeNodeFast* node, uint32_t depth)
  {
  }

  virtual void insertRecursive(TreeNodeFast* node, const std::string& motif, uint32_t depth)
  {
  }

  virtual void intersectRecursive(TreeNodeFast* to, const TreeNodeFast* from, uint32_t depth)
  {
  }

  MotifTreeBase(uint32_t _max_depth, Motifs& _motifs, const char* _name);
  virtual ~MotifTreeBase();

  virtual void setDomain(const string& _domain);
  virtual void traverse();
  virtual void traverseOut();
  virtual void insert(const std::string& motif);
  virtual void intersect(const MotifTreeBase* other);
};

