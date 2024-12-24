#include <iostream>
#include <vector>
#include <cmath>
#include <sched.h>
#include <fstream>
#include <cstdio>
#include <unistd.h>
#include <sys/time.h>
#include <pthread.h>

#include <tbb/tbb.h>
#include "tbb/task_arena.h"


#include "parEMS.hpp"

#include "utils.h"
using namespace std;
using namespace tbb;

// This can be auto generated when d > 3
string mapDE_PAR[64] = {
  "AAA", "AAT", "AAG", "AAC", "ATA", "ATT", "ATG", "ATC",
  "AGA", "AGT", "AGG", "AGC", "ACA", "ACT", "ACG", "ACC",
  "TAA", "TAT", "TAG", "TAC", "TTA", "TTT", "TTG", "TTC",
  "TGA", "TGT", "TGG", "TGC", "TCA", "TCT", "TCG", "TCC",
  "GAA", "GAT", "GAG", "GAC", "GTA", "GTT", "GTG", "GTC",
  "GGA", "GGT", "GGG", "GGC", "GCA", "GCT", "GCG", "GCC",
  "CAA", "CAT", "CAG", "CAC", "CTA", "CTT", "CTG", "CTC",
  "CGA", "CGT", "CGG", "CGC", "CCA", "CCT", "CCG", "CCC"
};

string ALLZERO;
string ALLONE;

thread_local int my_cpu = -1;

class pinning_observer final : public tbb::task_scheduler_observer
{
public:
  pinning_observer() { observe(true); } // activate the observer

  void on_scheduler_entry(bool worker) override
  {
    const auto number_of_slots = tbb::this_task_arena::max_concurrency();
    cout << number_of_slots << endl;
    cpu_set_t* mask = CPU_ALLOC(number_of_slots);
    const auto mask_size = CPU_ALLOC_SIZE(number_of_slots);

    const auto slot_number = tbb::this_task_arena::current_thread_index();
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

worker_parEms::worker_parEms(const int _id)
{
  m_cur_seq_id = _id;
}

worker_parEms::~worker_parEms()
{
  ;
}

void worker_parEms::gen_nbrhood3(const int start, const int alpha)
{
  printf("Entered gen_nbrhood3, start = %d, alpha = %d\n", start, alpha);
  const int len = static_cast<int>(kmer.size());;
  if (alpha > 0)
  {
    int j;
    for (j = start; j < len; j++)
    {
      if (kmer[j] == '*')
        continue;
      if (kmer[j] == '-')
        continue;
      if ((j > 0) && (kmer[j - 1] == '-'))
        continue;
      kmer.insert(j, 1, '*');
      gen_nbrhood3(j + 1, alpha - 1);
      kmer.erase(j, 1);
    }
    if (rightmost && (kmer[j - 1] != '-'))
    {
      kmer.insert(j, 1, '*');
      gen_nbrhood3(j + 1, alpha - 1);
      kmer.erase(j, 1);
    }
  }
  else
  {
    std::string t(kmer);
    int i = static_cast<int>(t.size() - 1);
    int num_possible = 1;
    int num_star = 0;
    while (i >= 0)
    {
      if (t[i] == '-')
        t.erase(i, 1);
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
        int startIndex = num_star - 1;
        for (char& letter : temp)
        {
          if (letter == '*')
          {
            printf("Special operation for letter '*' has been run\n");
            letter = mapDE_PAR[counter].at(startIndex);
            startIndex--;
          }
        }
        if (StringTable::accessor b; global_Hash.find(b, temp))
        {
          printf("Found temp in the global_Hash, b->second = %s\t", b->second.c_str());
          (b->second).at(m_cur_seq_id) = '1';
          printf("After setting, b->second = %s\n", b->second.c_str());
          //b->second = (b->second) | (1<<m_cur_seq_id);
        }
        else
        {
          printf("Not found temp in the global_Hash \t");
          global_Hash.insert(b, temp);
          // b->second = (1<<m_cur_seq_id);
          b->second = (ALLZERO);
          (b->second).at(m_cur_seq_id) = '1';
          printf("After setting, b->second = %s\n", b->second.c_str());
        }
      }
    }
    else
    {
      StringTable::accessor b;
      if (global_Hash.find(b, t))
      {
        printf("(alpha <= 0) Found temp in the global_Hash, b->second = %s\t", b->second.c_str());
        (b->second).at(m_cur_seq_id) = '1';
        printf("After setting, b->second = %s\n", b->second.c_str());
        // b->second = (b->second) | (1<<m_cur_seq_id);
      }
      else
      {
        global_Hash.insert(b, t);
        // b->second = (1<<m_cur_seq_id);

        printf("(alpha <= 0) Not found temp in the global_Hash\t");
        b->second = (ALLZERO);
        (b->second).at(m_cur_seq_id) = '1';
        printf("After setting, b->second = %s\n", b->second.c_str());
      }
    }
  }
  printf("Exiting gen_nbrhood3\n");
}

void worker_parEms::gen_nbrhood2(const int start, const int sigma, const int alpha)
{
  printf("Entered gen_nbrhood2, start = %d, sigma = %d, alpha = %d\n", start, sigma, alpha);
  const int len = static_cast<int>(kmer.size());
  if (sigma > 0)
  {
    for (int j = start; j < len; j++)
    {
      if (kmer[j] == '-')
      {
        while (kmer[++j] == '-' && j < len)
        {
        }
        continue;
      }
      if (!rightmost && (stars + 1 == j) && (kmer[j + 1] == '-'))
        continue;
      const char t = kmer[j];
      kmer[j] = '*';
      const int old_stars = stars;
      if (stars + 1 == j)
        stars++;
      gen_nbrhood2(j + 1, sigma - 1, alpha);
      kmer[j] = t;
      stars = old_stars;
    }
  }
  else
  {
    const int new_start = leftmost ? 0 : stars + 2;
    gen_nbrhood3(new_start, alpha);
  }
  printf("Exiting gen_nbrhood2\n");
}

void worker_parEms::setDomainSize(const int _domain_size)
{
  m_domain_size = _domain_size;
}

void worker_parEms::gen_nbrhood(const int start, const int delta, const int sigma, const int alpha)
{
  printf("Entered gen_nbrhood, start = %d, delta = %d, sigma = %d, alpha = %d\n", start, delta, sigma, alpha);
  const int len = static_cast<int>(kmer.size());
  if (delta > 0)
  {
    for (int j = start; j < len; j++)
    {
      const char t = kmer[j];
      //kmer.erase(j,1);
      kmer[j] = '-';
      gen_nbrhood(j + 1, delta - 1, sigma, alpha);
      //kmer.insert(j,1,t);
      kmer[j] = t;
    }
  }
  else
  {
    gen_nbrhood2(0, sigma, alpha);
  }
  printf("Exiting gen_nbrhood\n");
}

void worker_parEms::gen_all(const std::string& seq, const int motif_len, const int motif_distance)
{
  const int seq_len = static_cast<int>(seq.size());
  cout << "Start running on " << m_cur_seq_id << " " << seq << " " << motif_len << " " << motif_distance << endl;
  for (int motif_search_ptr = -motif_distance; motif_search_ptr <= +motif_distance; motif_search_ptr++)
  // TODO: Better name and interpretation required
  {
    const int k = motif_len + motif_search_ptr;
    for (int delta = std::max(0, motif_search_ptr); delta <= (motif_distance + motif_search_ptr) / 2; delta++)
    {
      const int alpha = delta - motif_search_ptr;
      const int sigma = motif_distance - alpha - delta;
      for (int i = 0; i < seq_len - k + 1; i++) // for each k-mer in seq
      {
        kmer = seq.substr(i, k); // seq[i:i+k]
        int start = 1;
        stars = -1;
        leftmost = rightmost = 0;
        if (i == 0)
          leftmost = 1;
        if (i + k == seq_len)
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
  if (const auto itrV = ranges::find(originalData, *p); itrV != originalData.end())
    cur_id = static_cast<int>(itrV - originalData.begin());
  else
  {
    cout << "Something wrong ... " << endl;
    cur_id = 0;
  }
  printf("Start Working on %d (on core %d)\n", cur_id, sched_getcpu());
  // cout << "Start Working on " << cur_id << endl;
  const auto obj = std::make_unique<worker_parEms>(cur_id);
  obj->setDomainSize(GLOBAL_DOMAIN);
  obj->gen_all(*p, l, d);
};

struct Task
{
  int l;
  int d;

  Task(const int l_, const int d_) : l(l_), d(d_)
  {
  }

  void operator()(const blocked_range<string*>& range) const
  {
    for (string* p = range.begin(); p != range.end(); ++p)
      runString(p, l, d);
  }
};

//Class ParEMS -----------------------------
ParEMS::ParEMS(const std::string& input, int l, int d, Params& params) :
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
  GLOBAL_DOMAIN = static_cast<int>(getDomain().size());
  const int num_seqs = static_cast<int>(reads.size());
  ALLZERO = string(num_seqs, '0');
  ALLONE = string(num_seqs, '1');
  auto* Data = new string [num_seqs]();
  const int loop_num = num_seqs;
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
  tbb::parallel_for(tbb::blocked_range<string*>(Data, Data + loop_num),
                    Task(getL(), getD()), ap);

  cout << "global tbb hash table size is " << global_Hash.size() * getL() << endl;
  // cout << "Should times " << num_seqs *1.0 / loop_num << endl;

  //int number = pow(2, num_seqs) - 1;
  for (auto& [fst, snd] : global_Hash)
  {
    //cout << i->second << endl;
    if (snd == ALLONE)
      motifs.push_back(fst);
  }
}
