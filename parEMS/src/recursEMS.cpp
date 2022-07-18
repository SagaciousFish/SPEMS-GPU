#include <iostream>
#include <pthread.h>
#include <string>
#include <sched.h>
#include<vector>
#include "omp.h"

#include "recursEMS.hpp"
#include "motif_tree_base.hpp"
#include "motif_tree_fast.hpp"

#include "utils.h"
using namespace std;

string mapDE[64] = {"AAA", "AAT","AAG","AAC", "ATA", "ATT","ATG","ATC",
"AGA", "AGT","AGG","AGC", "ACA", "ACT","ACG","ACC", 
"TAA", "TAT","TAG","TAC", "TTA", "TTT","TTG","TTC",
"TGA", "TGT","TGG","TGC", "TCA", "TCT","TCG","TCC",
"GAA", "GAT","GAG","GAC", "GTA", "GTT","GTG","GTC",
"GGA", "GGT","GGG","GGC", "GCA", "GCT","GCG","GCC",
"CAA", "CAT","CAG","CAC", "CTA", "CTT","CTG","CTC",
"CGA", "CGT","CGG","CGC", "CCA", "CCT","CCG","CCC"};

void worker_recursEms::run()
{
  double begin = omp_get_wtime(); 

  m_domain_size = m_caller->getDomain().size();
  Motifs tmp;
  MotifTreeFast tmp_tree(m_caller->getL(), tmp, "tmp");
  tmp_tree.setDomain(m_caller->getDomain());
  // cout << m_caller->getDomain() << endl;
  main_tree = &tmp_tree;
  // main_tree->setDomain(m_caller->getDomain());
  string cur = m_caller->getSeq(m_cur_seq_id);
  printf("Thread %d start working ...\n", sched_getcpu());
  gen_all(cur, m_caller->getL(), m_caller->getD());
  if (sched_getcpu() == 0)
    main_tree->print();
  printf("Thread %d start traversing ...\n", sched_getcpu());
  main_tree->traverse();
  copy(tmp.begin(), tmp.end(), back_inserter(m_mix_motifs));
  // main_tree->print();
  // std::sort(m_mix_motifs.begin(), m_mix_motifs.end());
  // m_mix_motifs.erase(unique( m_mix_motifs.begin(), m_mix_motifs.end() ), m_mix_motifs.end() );
  printf("The current sides of candidates is %d\n", m_mix_motifs.size());
  double end = omp_get_wtime();
  double elapsed = end-begin;
  printf("Time on Thread %d runs seq %d is %f secs. \n", sched_getcpu(), m_cur_seq_id, elapsed);
  // printf("Thread %d works on %d\n", sched_getcpu(), m_cur_seq_id);
  //cout << "Thread " << sched_getcpu() << " works on " << m_cur_seq_id << endl;

}

worker_recursEms::worker_recursEms()
{
  m_cur_seq_id = -1;
  // main_tree = NULL;
}

worker_recursEms::worker_recursEms(int l)
{
  m_cur_seq_id = -1;
  // main_tree = new MotifTreeFast(l, m_mix_motifs, "main");
}

worker_recursEms::~worker_recursEms()
{
  if (m_cur_seq_id != -1)
    m_mix_motifs.clear();
  // if (main_tree != NULL)
  //   delete main_tree;
}

void worker_recursEms::setTask(int _seq_id, recursEms* _caller)
{
  m_cur_seq_id = _seq_id;
  m_caller = _caller;
}


void worker_recursEms::setTree(int l)
{
  ;
  //main_tree = new MotifTreeFast(l, m_mix_motifs, "main");
}

// void worker_recursEms::gen_nbrhood3(int start, int alpha) 
// {
//   int len = kmer.size();
//   if (alpha > 0) 
//   {
//     int j;
//     for (j=start; j<len; j++) 
//     {
//       if (kmer[j] == '*') 
//         continue;
//       if (kmer[j] == '-') 
//         continue;
//       if ((j>0) && (kmer[j-1] == '-')) 
//         continue;
//       kmer.insert(j,1,'*');
//       gen_nbrhood3(j+1, alpha-1); 
//       kmer.erase(j,1);
//     }
//     if (rightmost && (kmer[j-1] != '-')) 
//     {
//       kmer.insert(j,1,'*');
//       gen_nbrhood3(j+1, alpha-1); 
//       kmer.erase(j,1);
//     }
//   } 
//   else 
//   {
//     std::string t(kmer);
//     int i=t.size()-1;
//     int num_possible = 1;
//     int num_star = 0;
//     while (i>=0) 
//     {
//       if (t[i] == '-') 
//         t.erase(i,1);
//       if (t[i] == '*') 
//       {
//         num_possible = num_possible * m_domain_size;
//         num_star++;
//       }
//       //   t[i] = 'm_domain_size';
//       i--;
//     }
//     //curr_tree->insert(t);
//     if (num_possible != 1)
//     {
//       for (int counter = 0; counter < num_possible; counter++)
//       {
//         string temp = t;
//         int startIndex = num_star-1;
//         for (int i = 0; i < temp.size(); i++)
//         {
//           if (temp.at(i) == '*')
//           {
//             temp.at(i) = mapDE[counter].at(startIndex);
//             startIndex--;
//           }          
//         }
//         m_mix_motifs.push_back(temp);
//       }          
//     }
//     else
//       m_mix_motifs.push_back(t);
//     // m_mix_motifs.push_back(t);
//   }
// }

// void worker_recursEms::gen_nbrhood2(int start, int sigma, int alpha) 
// {
//   int len = kmer.size();
//   if (sigma > 0) 
//   {
//     for (int j=start; j<len; j++) 
//     {
//       if (kmer[j] == '-') 
//       { 
//         while (kmer[++j] == '-'); 
//           continue; 
//       }
//       if (!rightmost && (stars+1 == j) && (kmer[j+1] == '-')) 
//         continue;
//       char t = kmer[j];
//       kmer[j] = '*';
//       int old_stars = stars;
//       if (stars+1 == j) 
//         stars++;
//       gen_nbrhood2(j+1, sigma-1, alpha); 
//       kmer[j] = t;
//       stars = old_stars;
//     }
//   } 
//   else 
//   {
//     int new_start = leftmost ? 0 : stars+2;
//     gen_nbrhood3(new_start, alpha);
//   }
// }

// void worker_recursEms::gen_nbrhood(int start, int delta, int sigma, int alpha) 
// {
//   int len = kmer.size();  
//   if (delta > 0) 
//   {
//     for (int j=start; j<len; j++) 
//     {
//       char t = kmer[j];
//       //kmer.erase(j,1);
//       kmer[j] = '-';
//       gen_nbrhood(j+1, delta-1,sigma,alpha);
//       //kmer.insert(j,1,t);
//       kmer[j] = t;
//     }
//   } 
//   else 
//   {
//     gen_nbrhood2(0, sigma, alpha);
//   }
// }

// void worker_recursEms::gen_all(std::string& seq, int l, int d) 
// {
//   int m = seq.size();
//   for (int q=-d; q<=+d; q++) 
//   {
//     int k = l+q;
//     for (int delta = std::max(0,q); delta <= (d+q)/2; delta++) 
//     {
//       int alpha = delta - q;
//       int sigma = d - alpha - delta;
//       for (int i=0; i<m-k+1; i++) 
//       {
//         kmer = seq.substr(i, k);
//         int start = 1; 
//         stars = -1;
//         leftmost = rightmost = 0;
//         if (i == 0) 
//           leftmost = 1; 
//         if (i+k == m) 
//         { 
//           start = 0; 
//           rightmost = 1; 
//         }
//         gen_nbrhood(start, delta, sigma, alpha);
//       }
//     }

//   }
// }

void worker_recursEms::gen_nbrhood3(int start, int alpha) 
{
  int len = kmer.size();
  if (alpha > 0) 
  {
  int j;
  for (j=start; j<len; j++) {
    if (kmer[j] == '*') continue;
    if (kmer[j] == '-') continue;
    if ((j>0) && (kmer[j-1] == '-')) continue;
    kmer.insert(j,1,'*');
      gen_nbrhood3(j+1, alpha-1); 
    kmer.erase(j,1);
  }
  if (rightmost && (kmer[j-1] != '-')) {
    kmer.insert(j,1,'*');
    gen_nbrhood3(j+1, alpha-1); 
    kmer.erase(j,1);
    }
  } else {
    std::string t(kmer);
    int i=t.size()-1;
    while (i>=0) {
      if (t[i] == '-') t.erase(i,1);
      if (t[i] == '*') t[i] = m_domain_size;
      i--;
    }
    main_tree->insert(t);
    //std::cout << "Insert one " << std::endl;
  }
}

void worker_recursEms::gen_nbrhood2(int start, int sigma, int alpha) 
{
  int len = kmer.size();
  if (sigma > 0) {
  for (int j=start; j<len; j++) {
    if (kmer[j] == '-') { while (kmer[++j] == '-'); continue; }
    if (!rightmost && (stars+1 == j) && (kmer[j+1] == '-')) continue;
    char t = kmer[j];
    kmer[j] = '*';
    int old_stars = stars;
    if (stars+1 == j) stars++;
    gen_nbrhood2(j+1, sigma-1, alpha); 
    kmer[j] = t;
    stars = old_stars;
    }
  } else {
  int new_start = leftmost ? 0 : stars+2;
  gen_nbrhood3(new_start, alpha);
  }
}

void worker_recursEms::gen_nbrhood(int start, int delta, int sigma, int alpha) 
{
  int len = kmer.size();
    if (delta > 0) {
    for (int j=start; j<len; j++) {
      char t = kmer[j];
      //kmer.erase(j,1);
      kmer[j] = '-';
      gen_nbrhood(j+1, delta-1,sigma,alpha);
      //kmer.insert(j,1,t);
      kmer[j] = t;
      }
    } else {
    gen_nbrhood2(0, sigma, alpha);
    }
}

void worker_recursEms::gen_all(std::string& seq, int l, int d) 
{
  int m = seq.size();
  for (int q=-d; q<=+d; q++) 
  {
    int k = l+q;
    for (int delta = std::max(0,q); delta <= (d+q)/2; delta++) {
      int alpha = delta - q;
      int sigma = d - alpha - delta;
      for (int i=0; i<m-k+1; i++) {
      kmer = seq.substr(i, k);
      int start = 1; stars = -1;
      leftmost = rightmost = 0;
      if (i == 0) { leftmost = 1; }
      if (i+k == m) { start = 0; rightmost = 1; }
      gen_nbrhood(start, delta, sigma, alpha);
      }
    }
  }
}


int worker_recursEms::getMotifsSize()
{
  return m_mix_motifs.size();
}

Motifs& worker_recursEms::getMotifs()
{
  return m_mix_motifs;
}

//Class recursEms -----------------------------
recursEms::recursEms(const std::string &input, int l, int d, Params &params) : 
  MotifFinder("recursEms", input, l, d, params)
{
  m_num_threads = params.num_threads;
  decodeStrings(reads, domain);
}

recursEms::~recursEms()
{
  // delete curr_tree;
}


void* recursEms::callFunc(void* args)
{
  worker_recursEms* Arg = (worker_recursEms*)args;
  Arg->run();
  pthread_exit((void*)args);
}  

void recursEms::search() 
{
  int num_seqs = reads.size();
  int parallel_num_seqs = ((int)(num_seqs/m_num_threads)) * m_num_threads;
  int rest = num_seqs - parallel_num_seqs;
  worker_recursEms* workers;
  workers = new worker_recursEms[num_seqs];  

  for (int c1 = 0; c1 < parallel_num_seqs; c1+= m_num_threads)
  {
    //---------------PTHREAD 
    pthread_t* threads;     
    threads=(pthread_t*)malloc(sizeof(pthread_t)* m_num_threads);
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setstacksize(&attr, 5504628000);
    cpu_set_t cpus; 

    for(int i = 0; i < m_num_threads; i++)
    {
      workers[c1+i].setTask(c1+i, this);
      // workers[c1+i].setTree(getL());
      CPU_ZERO(&cpus);
      CPU_SET(i, &cpus);  
      pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpus);
      size_t stksize;
      if (pthread_attr_getstacksize(&attr, &stksize) == -1) 
      {                         
        perror("error in pthread_attr_getstackstate()");                           
        exit(2);                                                                   
      }                                                                            
      pthread_create(&threads[i], &attr, callFunc, (void*)(&workers[c1+i]));         
    }
    for(int i=0; i< m_num_threads; i++)
        pthread_join(threads[i], NULL); 
    //---------------PTHREAD  
  }

  if (rest != 0)
  {
    //---------------PTHREAD 
    pthread_t* threads;     
    threads=(pthread_t*)malloc(sizeof(pthread_t)* rest);
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setstacksize(&attr, 5504628000);
    cpu_set_t cpus; 

    for(int i = 0; i < rest; i++)
    {
      workers[i+parallel_num_seqs].setTask(i+parallel_num_seqs, this);
      CPU_ZERO(&cpus);
      CPU_SET(i, &cpus);  
      pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpus);
      size_t stksize;
      if (pthread_attr_getstacksize(&attr, &stksize) == -1) 
      {                         
        perror("error in pthread_attr_getstackstate()");                           
        exit(2);                                                                   
      }                                                                            
      pthread_create(&threads[i], &attr, callFunc, (void*)(&workers[i+parallel_num_seqs]));         
    }
    for(int i = 0; i < rest; i++)
        pthread_join(threads[i], NULL); 
    //---------------PTHREAD      
  }

  Motifs curMotif0 = workers[0].getMotifs();
  Motifs curMotif1 = workers[1].getMotifs();
  std::set_intersection(curMotif0.begin(), curMotif0.end(), curMotif1.begin(), curMotif1.end(), 
        back_inserter(motifs));
  Motifs temp_motif1;
  std::set_intersection((workers[2].getMotifs()).begin(), (workers[2].getMotifs()).end(), 
    (workers[3].getMotifs()).begin(), (workers[3].getMotifs()).end(), 
        back_inserter(temp_motif1));
  Motifs temp_motif2;
  std::set_intersection(motifs.begin(), motifs.end(), 
   temp_motif1.begin(), temp_motif1.end(), 
        back_inserter(temp_motif2));
   std::swap(motifs, temp_motif2);  

  // omp_set_num_threads(m_num_threads);
  // #pragma omp parallel for
  for (int x = 4; x < num_seqs; ++x)
  {
    Motifs temp_motif;
    std::set_intersection(motifs.begin(), motifs.end(), (workers[x].getMotifs()).begin(), (workers[x].getMotifs()).end(), 
        back_inserter(temp_motif));  
    std::swap(motifs, temp_motif);  
    temp_motif.clear();
    if (motifs.size() == 0)
      break;     
  }

  // std::sort(motifs.begin(), motifs.end());
  // motifs.erase( unique( motifs.begin(), motifs.end() ), motifs.end() );
   

   //cout << mapDE[0].at(2) << mapDE[1].at(2) << mapDE[2].at(2) << endl;
  // for (int j = 0; j < motifs.size(); j++)
  //   cout << motifs[j].size() << " " << motifs[j] << endl;

}
