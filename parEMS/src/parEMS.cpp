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
#include <pthread.h>

#include <tbb/tbb.h>
#include "tbb/task_group.h"
#include "tbb/task_arena.h"


#include "parEMS.hpp"
#include "motif_tree_base.hpp"
#include "motif_tree_fast.hpp"

#include "utils.h"
using namespace std;
using namespace tbb;

// This can be auto generated when d > 3
string mapDE_PAR [64] = {"AAA", "AAT","AAG","AAC", "ATA", "ATT","ATG","ATC",
"AGA", "AGT","AGG","AGC", "ACA", "ACT","ACG","ACC", 
"TAA", "TAT","TAG","TAC", "TTA", "TTT","TTG","TTC",
"TGA", "TGT","TGG","TGC", "TCA", "TCT","TCG","TCC",
"GAA", "GAT","GAG","GAC", "GTA", "GTT","GTG","GTC",
"GGA", "GGT","GGG","GGC", "GCA", "GCT","GCG","GCC",
"CAA", "CAT","CAG","CAC", "CTA", "CTT","CTG","CTC",
"CGA", "CGT","CGG","CGC", "CCA", "CCT","CCG","CCC"};

string ALLZERO;
string ALLONE;

thread_local int my_cpu = -1;
class pinning_observer : public tbb::task_scheduler_observer 
{
public:
    pinning_observer() { observe(true); } // activate the observer}

    void on_scheduler_entry( bool worker ) override 
    {
      cpu_set_t *mask;
      auto number_of_slots = tbb::this_task_arena::max_concurrency();
      cout << number_of_slots << endl;
      mask = CPU_ALLOC(number_of_slots);
      auto mask_size = CPU_ALLOC_SIZE(number_of_slots);

      auto slot_number = tbb::this_task_arena::current_thread_index();
      CPU_ZERO_S(mask_size, mask);
      CPU_SET_S(slot_number, mask_size, mask);
      if (sched_setaffinity(0, mask_size, mask))
        printf("Error in sched_setaffinity. \n");
      my_cpu = sched_getcpu();
    }

};

int GLOBAL_DOMAIN;
StringTable global_Hash;
vector<string> originalData;

worker_parEms::worker_parEms()
{
  m_cur_seq_id = -1;
}

worker_parEms::worker_parEms(int _id)
{
  m_cur_seq_id = _id;
}

worker_parEms::~worker_parEms()
{
  ;
}

void worker_parEms::gen_nbrhood3(int start, int alpha) 
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
    int num_possible = 1;
    int num_star = 0;
    while (i>=0) 
    {
      if (t[i] == '-') 
        t.erase(i,1);
      if (t[i] == '*') 
      {
        num_possible = num_possible * m_domain_size;
        num_star++;
      }
      //   t[i] = 'm_domain_size';
      i--;
    }
    //curr_tree->insert(t);
    if (num_possible != 1)
    {
      for (int counter = 0; counter < num_possible; counter++)
      {
        string temp = t;
        int startIndex = num_star-1;
        for (int i = 0; i < temp.size(); i++)
        {
          if (temp.at(i) == '*')
          {
            temp.at(i) = mapDE_PAR[counter].at(startIndex);
            startIndex--;
          }          
        }
        StringTable::accessor b;
        if (global_Hash.find(b, temp) )
          (b->second).at(m_cur_seq_id)='1';
          //b->second = (b->second) | (1<<m_cur_seq_id);
        else
        {
          global_Hash.insert(b, temp);
          // b->second = (1<<m_cur_seq_id);
          b->second = (ALLZERO);
          (b->second).at(m_cur_seq_id)='1';
        }
      }          
    }
    else
    {
      StringTable::accessor b;
      if (global_Hash.find(b, t))
        (b->second).at(m_cur_seq_id)='1';
        // b->second = (b->second) | (1<<m_cur_seq_id);
      else
      {
        global_Hash.insert(b, t);
        // b->second = (1<<m_cur_seq_id);
          b->second = (ALLZERO);
          (b->second).at(m_cur_seq_id)='1';
      }      
    }
  }
}

void worker_parEms::gen_nbrhood2(int start, int sigma, int alpha) 
{
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

void worker_parEms::setDomainSize(int _domain_size)
{
  m_domain_size = _domain_size;
} 

void worker_parEms::gen_nbrhood(int start, int delta, int sigma, int alpha) 
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

void worker_parEms::gen_all(std::string& seq, int l, int d) 
{
  int m = seq.size();
  //cout << "Start running on " << m_cur_seq_id << " " << seq << " " << l << " " << d << endl;
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

void runString(string* p, int l, int d)
{ 
  int cur_id;
  auto itrV = find(originalData.begin(), originalData.end(), *p);
  if (itrV != originalData.end())
    cur_id = itrV - originalData.begin();
  else
  {
    cout << "Something wrong ... " << endl;
    cur_id = 0;
  }
  printf("Start Working on %d (on core %d)\n", cur_id, sched_getcpu());
  // cout << "Start Working on " << cur_id << endl;
  worker_parEms* obj;
  obj = new worker_parEms(cur_id);
  obj->setDomainSize(GLOBAL_DOMAIN);
  obj->gen_all(*p, l, d);
};

struct Task 
{
    int l;
    int d;
    Task(int l_, int d_) : l(l_), d(d_) {}

    void operator()( const blocked_range<string*> range ) const 
    {
        for( string* p=range.begin(); p!=range.end(); ++p ) 
          runString(p, l, d);
    }
};

//Class ParEMS -----------------------------
ParEMS::ParEMS(const std::string &input, int l, int d, Params &params) : 
  MotifFinder("ParEMS", input, l, d, params)
{
  m_num_threads = params.num_threads;
  decodeStrings(reads, domain);
}

ParEMS::~ParEMS()
{
  // delete main_tree;
  // delete curr_tree;
} 

void ParEMS::search() 
{
  GLOBAL_DOMAIN = getDomain().size();
  int num_seqs = reads.size();
  ALLZERO = string(num_seqs, '0');
  ALLONE = string(num_seqs, '1');
  string* Data;
  Data = new string [num_seqs]();  
  int loop_num = num_seqs;
  // int loop_num = m_num_threads; // num_seqs;
  // if (m_num_threads > num_seqs)
  //   loop_num = num_seqs;
  for (int i = 0; i < loop_num; i++)
  {
    Data[i] = reads.at(i);  
    originalData.push_back(reads.at(i));
  } 

  tbb::global_control control(tbb::global_control::max_allowed_parallelism, m_num_threads); 
  static affinity_partitioner ap;
  tbb::parallel_for(tbb::blocked_range<string*>(Data, Data+loop_num),
                Task(getL(), getD()), ap);

  cout << "global tbb hash table size is " << global_Hash.size() * getL()  << endl;
  // cout << "Should times " << num_seqs *1.0 / loop_num << endl;

  //int number = pow(2, num_seqs) - 1;
  for( StringTable::iterator i=global_Hash.begin(); i!=global_Hash.end(); ++i)
  {
    //cout << i->second << endl;
        if (i->second == ALLONE)
          motifs.push_back(i->first);    
  }
}