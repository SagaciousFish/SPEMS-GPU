#ifndef __EMS_HPP__
#define __EMS_HPP__

#include <sys/time.h>
#include <sys/resource.h>
#include <iostream>
#include <fstream>

#include "omp.h"
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
  size_t domain_size;

  std::string input; // the input file name
  Reads reads; // "vector<string> Reads" which is defined in utils.h, real input seq
  int l, d;
  Params &params;

  Motifs motifs; // "vector<string> Motifs" which is defined in utils.h, out seqs

public:

  MotifFinder(const std::string& _name, const Reads &_reads, 
    int _l, int _d, Params &_params) : name(_name), reads(_reads), 
      l(_l), d(_d), params(_params) { }

  MotifFinder(const std::string& _name, const std::string &_input, 
    int _l, int _d, Params &_params) : name(_name), input(_input), 
      l(_l), d(_d), params(_params) 
  {
    read_file(input.c_str(), reads);
    domain = getAlphabet(reads);
    domain_size = domain.size();
    encodeStrings(reads, domain);
  }

  MotifFinder(const std::string& _name, Reads _cur_seqs, std::string _alphabet,  
    int _l, int _d, Params &_params) : name(_name), l(_l), d(_d), params(_params)
  {}

  ~MotifFinder() 
  {
    reads.clear();
    motifs.clear();
  }

  virtual void search(){}
  Motifs& searchGetMotifs();
  void searchWriteMotifs(Params &params);
};

#endif // __EMS_HPP__



