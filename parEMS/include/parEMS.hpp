#ifndef __PAR_EMS_H_
#define __PAR_EMS_H_
#include <iostream>
#include <vector>
#include <cmath>
#include <sched.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <bitset>
#include <sys/resource.h>
#include <sys/time.h>

#include <tbb/tbb.h>

#include "omp.h"
#include "ems.hpp"
#include "ems2.hpp"
#include "motif_tree_base.hpp"
#include "motif_tree_fast.hpp"
#include "utils.h"

typedef tbb::concurrent_hash_map<string,string> StringTable;

class ParEMS : public MotifFinder
{
protected:
  int m_num_threads;

// protected:
//   //void taskAssignment(Task* assigned_task);
//   static void* callFunc(void* args);

public:

  ParEMS(const std::string &input, int l, int d, Params &params);
  ~ParEMS();
  void search();

  int getL() {return l;} 
  int getD() {return d;} 
  string getDomain() {return domain; }
  Params& getParams() {return params;}
  string getSeq(int id) {return reads[id];}
  string getInputFile() {return input;}

};

class worker_parEms
{
protected:
  ParEMS* m_caller;
  int m_cur_seq_id;

protected:
  void gen_nbrhood3(int start, int alpha);
  void gen_nbrhood2(int start, int sigma, int alpha);
  void gen_nbrhood(int start, int delta, int sigma, int alpha);

  int leftmost;
  int rightmost;
  std::string kmer;
  int stars;
  int m_domain_size;

public:
  void gen_all(std::string& seq, int l, int d);
  void setDomainSize(int _domain_size);

  worker_parEms();
  worker_parEms(int _id);
  ~worker_parEms();
};


#endif // __PAR_EMS_H_
