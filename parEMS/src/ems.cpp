#include <sys/time.h>
#include <sys/resource.h>
#include <iostream>
#include <fstream>

#include "omp.h"
#include "ems.hpp"
#include "utils.h"

using namespace std;

Motifs& MotifFinder::searchGetMotifs()
{
    search();
    return motifs;
}

void MotifFinder::searchWriteMotifs(const Params& params)
{
    const std::string output = get_out_file(input, l, d, name);
    std::cout << "l = " << l << ", d = " << d << ", t = " << params.num_threads << std::endl;
    cout << "Input stream number is " << reads.size() << ". " << endl;
    cout << "Domain includes:" << domain << ", and size is " << domain_size << endl;
    std::cout << "input File = " << input << std::endl;
    std::cout << "output File = " << output << std::endl;
    // clock_t begin=clock();
    // omp_set_num_threads(4);

    const double begin = omp_get_wtime();

    std::cout << "Start Measurement ... " << std::endl;

    this->search();
    // clock_t end=clock();

    const double end = omp_get_wtime();
    //double elapsed = diffclock(end,begin);

    const double elapsed = end - begin;
    struct rusage usage{};
    getrusage(RUSAGE_SELF, &usage);

    // ofstream myFile;
    // myFile.open("emsTimeMemory",ios::app);
    // myFile << name << ": (" << l << "," << d << ") Edited Motifs found using ";
    // myFile << params.num_threads << " threads:(in "<< elapsed << " sec, using ";
    // myFile << (size_t)usage.ru_maxrss << " KB): " << motifs.size() << "       " << std::endl;
    // myFile.close();    

    cout << endl;
    cout << name << ": (" << l << "," << d << ") Edited Motifs found is " << motifs.size() << endl;
    cout << "TOTAL Time spent is " << elapsed << " sec" << endl;
    //cout << "TOTAL Memory spent is " << (size_t)usage.ru_maxrss << " KB" << endl;
    const double MEM_GB = static_cast<double>(usage.ru_maxrss) * 1.0 / 1000000.0;
    cout << "TOTAL Memory spent is " << MEM_GB << " GB" << endl;

    // for (size_t i=0; i<motifs.size(); ++i) 
    // {
    //   cout << motifs[i] << endl;
    // }

    // std::ofstream out(output);
    // for (size_t i=0; i<motifs.size(); ++i) 
    // {
    //   out << motifs[i] << std::endl;
    // }
    // out.close();
}
