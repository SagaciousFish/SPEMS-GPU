#include <iostream>
#include "ems2.hpp"
#include "motif_tree_base.hpp"
#include "motif_tree_fast.hpp"

#include "utils.h"
using namespace std;

Ems2::Ems2(const std::string &input, int l, int d, Params &params) : 
  MotifFinder("ems2", input, l, d, params)
{
  main_tree = new MotifTreeFast(l, motifs, "main");
  curr_tree = main_tree; 
  main_tree->setDomain(domain);
}

Ems2::Ems2(const std::string& _name, Reads _cur_seqs, std::string _alphabet,  int _l, int _d, Params &_params) :
  MotifFinder("ems2", _cur_seqs, _alphabet, _l, _d, _params)
{
  for (int i = 0; i < _cur_seqs.size(); i++)
    reads.push_back(_cur_seqs[i]);
  
  domain = _alphabet;
  domain_size = domain.size();
  //encodeStrings(reads, domain);  

  main_tree = new MotifTreeFast(l, motifs, "main");
  curr_tree = main_tree; 
  main_tree->setDomain(domain);  
}


Ems2::~Ems2()
{
  // delete main_tree;
  // delete curr_tree;
}

void Ems2::write_out_motifs(Motifs &_cur_motif)
{
  for (size_t i=0; i<motifs.size(); ++i) 
    _cur_motif.push_back(motifs[i]); 
}

void Ems2::search() 
{
  cout << "Processing sequence " << 0 << "..." << endl;

  // for (int j = 0; j < reads.size(); j++)
  //   std::cout << reads[j].size() << std::endl;

  std::cout.flush();
  gen_all(reads[0], l, d);

  for (size_t i=1; i<reads.size(); i++) 
  {
    cout << "Processing sequence " << i << "..." << endl;
    std::cout.flush();
    Motifs tmp;
    MotifTreeFast tmp_tree(l, tmp, "tmp");
    tmp_tree.setDomain(domain);
    curr_tree = &tmp_tree;
    gen_all(reads[i], l, d);
    main_tree->intersect(curr_tree);
  }
  main_tree->traverse();
}

void Ems2::search1(std::string& _seq, int _l, int _d) 
{
  gen_all(_seq, _l, _d);
  //main_tree->intersect(curr_tree);
  main_tree->traverse();
}

void Ems2::gen_nbrhood3(int start, int alpha) 
{
  int len = kmer.size();
  if (alpha > 0) 
  {
    int j;
    for (j=start; j<len; j++) 
    {
      if (kmer[j] == '*') 
        continue;
      if (kmer[j] == '-') 
        continue;
      if ((j>0) && (kmer[j-1] == '-')) 
        continue;
      kmer.insert(j,1,'*');
      gen_nbrhood3(j+1, alpha-1); 
      kmer.erase(j,1);
    }
    if (rightmost && (kmer[j-1] != '-')) 
    {
      kmer.insert(j,1,'*');
      gen_nbrhood3(j+1, alpha-1); 
      kmer.erase(j,1);
    }
  } 
  else 
  {
    std::string t(kmer);
    int i=t.size()-1;
    while (i>=0) 
    {
      if (t[i] == '-') 
        t.erase(i,1);
      if (t[i] == '*') 
        t[i] = domain_size;
      i--;
    }
    curr_tree->insert(t);
  }
}

void Ems2::gen_nbrhood2(int start, int sigma, int alpha) {
  int len = kmer.size();
  if (sigma > 0) 
  {
    for (int j=start; j<len; j++) 
    {
      if (kmer[j] == '-') 
      { 
        while (kmer[++j] == '-'); 
          continue; 
      }
      if (!rightmost && (stars+1 == j) && (kmer[j+1] == '-')) 
        continue;
      char t = kmer[j];
      kmer[j] = '*';
      int old_stars = stars;
      if (stars+1 == j) 
        stars++;
      gen_nbrhood2(j+1, sigma-1, alpha); 
      kmer[j] = t;
      stars = old_stars;
    }
  } 
  else 
  {
    int new_start = leftmost ? 0 : stars+2;
    gen_nbrhood3(new_start, alpha);
  }
}

void Ems2::gen_nbrhood(int start, int delta, int sigma, int alpha) 
{
  int len = kmer.size();  
  if (delta > 0) 
  {
    for (int j=start; j<len; j++) 
    {
      char t = kmer[j];
      //kmer.erase(j,1);
      kmer[j] = '-';
      gen_nbrhood(j+1, delta-1,sigma,alpha);
      //kmer.insert(j,1,t);
      kmer[j] = t;
    }
  } 
  else 
  {
    gen_nbrhood2(0, sigma, alpha);
  }
}

void Ems2::gen_all(std::string& seq, int l, int d) 
{
  int m = seq.size();
  for (int q=-d; q<=+d; q++) 
  {
    int k = l+q;
    for (int delta = std::max(0,q); delta <= (d+q)/2; delta++) 
    {
      int alpha = delta - q;
      int sigma = d - alpha - delta;
      for (int i=0; i<m-k+1; i++) 
      {
        kmer = seq.substr(i, k);
        int start = 1; 
        stars = -1;
        leftmost = rightmost = 0;
        if (i == 0) 
          leftmost = 1; 
        if (i+k == m) 
        { 
          start = 0; 
          rightmost = 1; 
        }
        gen_nbrhood(start, delta, sigma, alpha);
      }
    }

  }
}
