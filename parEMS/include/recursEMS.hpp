#pragma once __RECURS_EMS_H_

#include "ems.hpp"
#include "ems2.hpp"
#include "motif_tree_base.hpp"
#include "motif_tree_fast.hpp"
#include "utils.h"

class recursEms final : public MotifFinder
{
protected:
  int m_num_threads;

protected:
  //void taskAssignment(Task* assigned_task);
  static void* callFunc(void* args);

public:
  recursEms(const std::string& input, int l, int d, Params& params);
  ~recursEms() override;
  void search() override;

  [[nodiscard]] int getL() const { return l; }
  [[nodiscard]] int getD() const { return d; }
  string getDomain() { return domain; }
  [[nodiscard]] Params& getParams() const { return params; }
  string getSeq(const int id) { return reads[id]; }
  string getInputFile() { return input; }
};

class worker_recursEms
{
protected:
  Motifs m_mix_motifs;
  recursEms* m_caller;
  int m_cur_seq_id;
  MotifTreeFast* main_tree;

protected:
  void gen_nbrhood3(int start, int alpha);
  void gen_nbrhood2(int start, int sigma, int alpha);
  void gen_nbrhood(int start, int delta, int sigma, int alpha);
  void gen_all(std::string& seq, int l, int d);


  int leftmost;
  int rightmost;
  std::string kmer;
  int stars;
  int m_domain_size;

public:
  void run();
  void setTask(int _seq_id, recursEms* _caller);
  void setTree(int l);
  int getMotifsSize();
  Motifs& getMotifs();

  worker_recursEms();
  explicit worker_recursEms(int l);
  ~worker_recursEms();
};
