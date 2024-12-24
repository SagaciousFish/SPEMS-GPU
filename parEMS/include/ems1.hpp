#pragma once __EMS1_H_

#include "ems.hpp"

class Ems1 final : public MotifFinder
{
protected:
  uint8_t *A1, *A2;
  uint8_t current_seq;
  uint64_t num_lmer;

  void printString(string& s);
  string to_str(string& s);
  std::string to_str(uint64_t id, int l);
  uint64_t to_int(std::string s);

  void gen_nbrhood(std::string& x, int errors);
  void gen_all(std::string& x);

public:
  Ems1(const std::string& input, int l, int d, Params& params);
  ~Ems1() override;
  void search() override;
};

