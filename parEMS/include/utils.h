#ifndef __UTILS_
#define __UTILS_

#include<string>
#include<vector>
#include<list>
#include<algorithm>
#include<sstream>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<csignal>
#include<cstring>
#include<cstdint>
#include<cmath>
#include<ctime>
#include<cassert>

using namespace std;

typedef std::vector<std::string> Reads;
typedef std::vector<std::string> Motifs;

typedef unsigned char uchar;
typedef unsigned int uint32;
typedef unsigned long long uint64;

//const std::string domain="ACGT";
//const size_t domain_size = domain.size();


// Read strings from file and put them into vector @reads
void read_file(const char *fname, Reads &reads);

string getAlphabet(vector<string>& strings);

void encodeString(string& s, string &sigma);

void encodeStrings(vector<string>& strings, string &sigma);

void decodeStrings(vector<string>& strings, string &sigma);
void decodeString(string& s, string &sigma);

bool diff_motifs(Motifs & m1, Motifs & m2);

void write_to_file(int l, int d, Motifs &m, std::string fname);

void write_to_file(Motifs &m, std::string fname);

int edist(std::string& s1, std::string& s2);

bool has_overlap(std::string x, std::string y, int common);

bool found_in_seq(std::string& candi, const std::string& seq, int l, int d);

inline double diffclock(clock_t clock1,clock_t clock2);

inline void show_progress(size_t done, size_t todo, clock_t start_clk);

std::string removeExtension(const std::string & filename);

std::string get_out_file(const std::string& input, int l, int d, const std::string & prefix="");

std::string get_out_file(const std::string& input, int l, int d, int ldash, const std::string & prefix="");

template<typename T>
void printList(const std::vector<std::vector<T>* >& l, const std::string& msg="", int m=-1) 
{
    std::cout << "================ " << msg << " ==================" << std::endl;
    //if (m==-1) m = l.size();
    for (auto &v : l) 
    {
      for (auto const &i : *v) 
      {
        std::cout << std::setw(4) << (int)i;
      }
      std::cout << std::endl;
      //if (--m == 0) return;
    }
}

#endif // __UTILS_
