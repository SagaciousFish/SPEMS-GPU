#ifndef __EMS2_H_
#define __EMS2_H_

#include "ems.hpp"
#include "motif_tree_base.hpp"
#include "motif_tree_fast.hpp"


class Ems2 : public MotifFinder
{
protected:

  MotifTreeFast *main_tree;
  MotifTreeFast *curr_tree;

  size_t count = 0;
  int leftmost;
  int rightmost;
  std::string kmer;
  std::vector<std::string> nbrs;
  int stars;
  double gen_time=0.0;

  void gen_nbrhood3(int start, int alpha);
  void gen_nbrhood2(int start, int sigma, int alpha);
  void gen_nbrhood(int start, int delta, int sigma, int alpha);
  void gen_all(std::string& seq, int l, int d);

public:

  Ems2(const std::string &input, int l, int d, Params &params);
  Ems2(const std::string& _name, Reads _cur_seqs, 
    std::string _alphabet,  int _l, int _d, Params &_params);
  ~Ems2();
  void search();
  void search1(std::string& seq, int l, int d);
  
  void write_out_motifs(Motifs &_cur_motif);
};

#endif // __EMS2_H_

