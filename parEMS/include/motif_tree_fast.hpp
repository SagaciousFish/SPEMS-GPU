#pragma once __MOTIF_TREE_FAST__

#include "motif_tree_base.hpp"

class MotifTreeFast final : public MotifTreeBase
{
  void print_recursive(const TreeNodeFast* t);
  void deleteNode(TreeNodeFast* node);
  bool emptyNode(const TreeNodeFast* node);
  void copy(TreeNodeFast* to, const TreeNodeFast* from);
  void copyUnion(TreeNodeFast* to, const TreeNodeFast* from, const std::string& motif, size_t depth);
  void copyIntersection(TreeNodeFast* to, const TreeNodeFast* from, const TreeNodeFast* other, size_t depth);
  void insertCommonRecursive(TreeNodeFast* to, const TreeNodeFast* from, const std::string& motif, size_t depth);
  bool hasIntersect(TreeNodeFast* node, const std::string& motif, size_t depth);

public:
  MotifTreeFast(uint32_t _max_depth, Motifs& _motifs, const char* _name) : MotifTreeBase(_max_depth, _motifs, _name)
  {
  }

  ~MotifTreeFast() override
  {
    deleteNode(root);
  }

  void print();
  void traverseRecursive(TreeNodeFast* node, size_t depth);
  void insertRecursive(TreeNodeFast* node, const std::string& motif, size_t depth);
  void insertRecursiveNew(TreeNodeFast* node, const std::string& motif, size_t depth);
  void intersectRecursive(TreeNodeFast* to, const TreeNodeFast* from, size_t depth);

  void setDomain(const string& _domain) override;
  void traverse() override;
  void traverseOut() override;
  void insert(const std::string& motif) override;
  void intersect(const MotifTreeFast* other);
};

