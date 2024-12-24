#pragma once  //__PAR_EMS_H_

#include <vector>

#include <tbb/tbb.h>

#include "ems.hpp"

typedef tbb::concurrent_hash_map<string, string> StringTable;

class ParEMS final : public MotifFinder
{
protected:
  int m_num_threads;

  // protected:
  //   //void taskAssignment(Task* assigned_task);
  //   static void* callFunc(void* args);

public:
  ParEMS(const std::string& input, int l, int d, Params& params);
  ~ParEMS() override;
  void search() override;

  [[nodiscard]] int getL() const { return l; }
  [[nodiscard]] int getD() const { return d; }
  string getDomain() { return domain; }
  [[nodiscard]] Params& getParams() const { return params; }
  string getSeq(const int id) { return reads[id]; }
  string getInputFile() { return input; }
};

class worker_parEms
{
protected:
  ParEMS* m_caller{};
  int m_cur_seq_id;

protected:
  void gen_nbrhood3(int start, int alpha);
  void gen_nbrhood2(int start, int sigma, int alpha);
  void gen_nbrhood(int start, int delta, int sigma, int alpha);

  int leftmost{};
  int rightmost{};
  std::string kmer;
  int stars{};
  int m_domain_size{};

public:
  void gen_all(const std::string& seq, int motif_len, int motif_distance);
  void setDomainSize(int _domain_size);

  worker_parEms();
  explicit worker_parEms(int _id);
  ~worker_parEms();
};

