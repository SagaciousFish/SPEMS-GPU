#pragma once __EMS_HPP__

#include <utility>

#include "utils.h"

struct Params
{
  int l;
  int d;
  int num_threads;
};

class MotifFinder
{
protected:
  std::string name; // the current motif method
  std::string domain;
  size_t domain_size{};

  std::string input; // the input file name
  Reads reads; // "vector<string> Reads" which is defined in utils.h, real input
  // seq
  int l, d;
  Params& params;

  Motifs
  motifs; // "vector<string> Motifs" which is defined in utils.h, out seqs

public:
  MotifFinder(std::string _name, Reads _reads, const int _l, const int _d,
              Params& _params)
    : name(std::move(_name)), reads(std::move(_reads)), l(_l), d(_d),
      params(_params)
  {
  }

  MotifFinder(std::string _name, std::string _input, const int _l, const int _d,
              Params& _params)
    : name(std::move(_name)), input(std::move(_input)), l(_l), d(_d),
      params(_params)
  {
    read_file(input.c_str(), reads);
    domain = getAlphabet(reads);
    domain_size = domain.size();
    encodeStrings(reads, domain);
  }

  MotifFinder(std::string _name, const Reads& _cur_seqs,
              const std::string& _alphabet, const int _l, const int _d,
              Params& _params)
    : name(std::move(_name)), l(_l), d(_d), params(_params)
  {
  }

  virtual ~MotifFinder()
  {
    reads.clear();
    motifs.clear();
  }

  virtual void search()
  {
  }

  Motifs& searchGetMotifs();
  void searchWriteMotifs(const Params& params);
};
