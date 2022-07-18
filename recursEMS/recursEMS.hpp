#ifndef __RECURS_EMS_H_
#define __RECURS_EMS_H_

#include <iostream>
#include <pthread.h>
#include <string>
#include <sched.h>
#include <vector>
#include <unistd.h>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/mman.h>
#include <chrono>
#include <thread>

#include "omp.h"
#include "ems2.hpp"
#include "ems.hpp"
#include "omp.h"
#include "motif_tree.hpp"
#include "motif_tree_fast.hpp"


using namespace std;

#define MAX_SHM_ELEMENT 800000

static char* sharedMemory;
static int* ptr;
static int* passNumMotifs;
#define errExit(msg) do { perror(msg); exit(EXIT_FAILURE); \
                               } while (0)


class worker_recursEms
{
public:
	worker_recursEms(){}
	~worker_recursEms()
	{
		m_mix_motifs.clear();
	}

	void setID(int i)
	{
		seq_id = i;
	}

	void setTask(const std::string &input, int l, int d, Params &params)
	{
		mL = l;
		mD = d;
		file_name = input;
		mParam = params;		
	}	

	void run()
	{
		// generating the motifs via EMS2
		// printf("Current Thread is %d, running target is %d \n", sched_getcpu(), seq_id);
		Ems2<MotifTreeFast> ems(file_name, mL, mD, mParam);
		ems.setTarget(seq_id);
		ems.searchWriteMotifs(mParam);
		cout << "Current Processor " << seq_id << " Finds motifs " << (ems.getMotifs()).size() << endl;

		// Collecting the memory usage for each process
		while (passNumMotifs[0] != (seq_id/2))
			std::this_thread::sleep_for(std::chrono::nanoseconds(1));
    	
    	//cout << "Current Processor " << seq_id << " Finds motifs " << (ems.getMotifs()).size() << endl;
    	//cout << "Start Checking" << sched_getcpu() << endl;
    	Motifs existing;
    	if (passNumMotifs[0] != 0) // the current id
    	{
    		if (passNumMotifs[1] != 0)
    		{
    			if (passNumMotifs[1] <= MAX_SHM_ELEMENT)
    			{
		    		for (int counter = 0; counter < passNumMotifs[1]; counter++)
		    			existing.push_back(std::string(&sharedMemory[counter*mL], &sharedMemory[(counter+1)*mL]));
    			}
    			else
    			{
    				ifstream inTemp;
    				inTemp.open("uselessTemp.txt");
		    		for (int counter = 0; counter < passNumMotifs[1]; counter++)
		    		{
		    			string temp_str;
		    			inTemp >> temp_str; 
		    			existing.push_back(temp_str);
		    		}
    				inTemp.close();
    			}

	    		Motifs intersections;
	    		std::set_intersection((ems.getMotifs()).begin(),(ems.getMotifs()).end(),
	                          existing.begin(),existing.end(), back_inserter(intersections));  
	    		// cout << intersections.size() << endl;
	    		if (intersections.size() > MAX_SHM_ELEMENT)
	    		{
	    			ofstream outFile_temp2;
	    			outFile_temp2.open("uselessTemp.txt");
	    			for (int counter3 = 0; counter3 < intersections.size(); counter3++)
	    				outFile_temp2 << intersections[counter3] << endl;
	    			outFile_temp2.close();
	    		}
	    		else
	    		{
		    		if (intersections.size() != 0)
		    		{
		    			for (int counter2 = 0; counter2 < intersections.size(); counter2++)
		    				strcpy(&sharedMemory[counter2 * mL], intersections[counter2].c_str());
		    		}	    			
	    		}

	    		cout << "Intersection is " << intersections.size() << endl;
				passNumMotifs[1] = intersections.size();
    		}
    	}
    	else
    	{
    		if ((ems.getMotifs()).size() > MAX_SHM_ELEMENT)
    		{
    			ofstream outFile_temp;
    			outFile_temp.open("uselessTemp.txt");
    			for (int counter1 = 0; counter1 < (ems.getMotifs()).size(); counter1++)
    				outFile_temp << (ems.getMotifs())[counter1] << endl;
    			// cout << "PayAttention, something may be wrong ... " << endl;   
    			outFile_temp.close();
    		}
    		else
    		{
	     		if ((ems.getMotifs()).size() != 0)
	    		{
	    			for (int counter1 = 0; counter1 < (ems.getMotifs()).size(); counter1++)
	    				strcpy(&sharedMemory[counter1 * mL], (ems.getMotifs())[counter1].c_str());
	    		}			
    		}
    		passNumMotifs[1] = (ems.getMotifs()).size();   

    	}
		passNumMotifs[0]++;		
		//cout << "Done work " << sched_getcpu() << endl;
	}

	int mL;
	int mD;
	string file_name;
	Params mParam;
	int seq_id;
	Motifs m_mix_motifs;
	long memoryUsed;
};

class recursEms
{
public:
	recursEms(const std::string &input, int l, int d, Params &params);
	~recursEms();
	void run();
	int getNumSeqs(string file_name);

protected:
	int mL;
	int mD;
	string file_name;
	Params mParam;
	int m_num_threads; 
	void callFunc(void* args);	
};

recursEms::recursEms(const std::string &input, int l, int d, Params &params)
{
	mL = l;
	mD = d;
	file_name = input;
	mParam = params;
}

recursEms::~recursEms()
{
	;
}

int recursEms::getNumSeqs(string file_name) 
{
	std::ifstream f1(file_name);
	if (!f1.is_open()) 
	{
		std::cerr << "ERROR: could not open file " << file_name << "\n";
		exit(-1);
	}
	std::string p1;
	int counter = 0;
	while(std::getline(f1,p1)) 
	{
		if ((p1.size() > 0) && (p1[0] != '>')) 
			counter++;
	}
	f1.close();
	return counter;
}

void recursEms::run()
{
	struct rusage usage;
	double begin = omp_get_wtime();

    m_num_threads = mParam.num_threads;
    
    int num_seqs = getNumSeqs(file_name);
    int num_seqs_pair = num_seqs / 2;
    if (num_seqs % 2 != 0)
    {
    	num_seqs_pair++;
    }

  	worker_recursEms* workers;
  	workers = new worker_recursEms[num_seqs];      
  
  	for (int cc = 0; cc < num_seqs_pair; cc++)
  	{
  		workers[cc].setTask(file_name, mL, mD, mParam);
  		workers[cc].setID(cc*2);
  		workers[cc].memoryUsed = 0;
  	}

    ptr = static_cast<int*>(mmap(NULL,(num_seqs_pair+1)*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0,0));    
    if(ptr == MAP_FAILED)
    {
    	printf("Mapping Failed\n");
    	exit(-1);
    }
    for(int i=0; i < num_seqs_pair+1; i++)
    	ptr[i] = 0;

    passNumMotifs = static_cast<int*>(mmap(NULL,2*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0,0));    
    if(passNumMotifs == MAP_FAILED)
    {
    	printf("Pass Counter Mapping Failed\n");
    	exit(-1);
    }
    for(int i=0; i < 2; i++)
    	passNumMotifs[i] = 0;    

    sharedMemory = static_cast<char*>(mmap(NULL,(MAX_SHM_ELEMENT * mL) * sizeof(char), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0,0));    
    if(sharedMemory == MAP_FAILED)
    {
    	printf("Shared Memory Mapping Failed\n");
    	exit(-1);
    }

    // Start Parallel Computations for Motif Search
	for (int c1 = 0; c1 < num_seqs_pair; c1+= m_num_threads)
	{
		pid_t* pids;
		pids = new pid_t [m_num_threads];		
		int i;
		int n = m_num_threads;
		
		for (i = 0; i < m_num_threads; i++) 
		{

			if (c1+i >= num_seqs_pair)
			{
				n = i+1;
				break;
			}

			int gpuID = i;
			if ((pids[i] = fork()) < 0) 
			{
				perror("fork");
				abort();
			} 
			else if (pids[i] == 0) 
			{
				cpu_set_t set;
				CPU_ZERO(&set);
				CPU_SET(gpuID, &set);
				if (sched_setaffinity(getpid(), sizeof(set), &set) == -1)
					errExit("sched_setaffinity");
				
				callFunc((void*)&workers[c1+i]);
				exit(0);
			}
		}

		int status;
		pid_t pid;
		while (n > 0) 
		{
		  pid = wait(&status);
		  --n;  
		}  
	}

    double end = omp_get_wtime();
    double elapsed = end-begin;
   	long MEM_GB = getrusage( RUSAGE_SELF, &usage);
   	for (int j = 1; j < num_seqs_pair+1; j++)
   	{
   		//cout << "Current spent is "<< j << " " << ptr[j]/1000000.0 << " GB" << endl;
    	MEM_GB += ptr[j];
   	}
   	MEM_GB /= ((int)(num_seqs/ m_num_threads )+1);
   	cout << mL << " " << mD << " " << file_name << endl;
    cout << "TOTAL Time spent is " << elapsed << " sec"<< endl;
    cout << "TOTAL Memory spent is " << MEM_GB/1000000.0 << " GB" << endl;
    cout << "Final Motif Number is  " << passNumMotifs[1] << endl << endl;

    int err = munmap(ptr, (num_seqs_pair+1) * sizeof(int));
    if(err != 0)
    {
    	printf("UnMapping Failed\n");
    	exit(-1);
    }    
    int err2 = munmap(sharedMemory, (MAX_SHM_ELEMENT * mL) * sizeof(char));
    if(err2 != 0)
    {
    	printf("UnMapping SHM Failed\n");
    	exit(-1);
    } 
    int err3 = munmap(passNumMotifs, 2 * sizeof(int));
    if(err3 != 0)
    {
    	printf("UnMapping Pass Number Failed\n");
    	exit(-1);
    }        
}

void recursEms::callFunc(void* args)
{
	struct rusage usage;
	worker_recursEms* Arg = (worker_recursEms*)args;
	Arg->run();
	getrusage( RUSAGE_SELF, &usage );

	// Collecting the memory usage for each process
	while (ptr[0] != ((Arg->seq_id)/2))
		std::this_thread::sleep_for(std::chrono::nanoseconds(1));
	ptr[((Arg->seq_id)/2)+1] = usage.ru_maxrss;
	ptr[0]++;
} 

#endif