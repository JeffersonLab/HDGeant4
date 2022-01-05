//
// samplesep - computes the Zech-Aslan two-sample "energy" distribution
//             that can be used to test if multi-dimensional samples A
//             and B might have come from the same parent distribution.
//
// author: richard.t.jones at uconn.edu
// version: october 27, 2021
//
// reference: G. Zech and B. Aslan, "A Multivariate Two-Sample Test 
//            Based on the Concept of Minimum Energy", PHYSTAT2003,
//            SLAC, Stanford, California, September 8-11, 2003.

#define USE_OPENMP_MULTITHREADING 1

#include <string.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <math.h>

#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TH1D.h>

#include <omp.h>

#define sqr(x) ((x)*(x))

#define GCC_VERSION (__GNUC__ * 10000 \
                     + __GNUC_MINOR__ * 100 \
                     + __GNUC_PATCHLEVEL__)

/* Test for GCC > 6.3.0 */
#if GCC_VERSION < 60300
#error Minimum required GNU compiler version is 6.3
#else

std::map<std::string, std::vector<std::vector<double> > > sample;

unsigned int Npartitions(1024);
unsigned int random_seed(0);
unsigned int work_sets(1);
unsigned int work_set(0);

void usage()
{
   std::cout << "Usage: samplesep [options] <input_file1>:<sample1> "
             << " [<input_file2>]:<sample2>" << std::endl
             << "   where <input_file1> is a ROOT file containing" << std::endl
             << "   a tree named <sample1> containing the first" << std::endl
             << "   sample, and <input_file2> is a second root file" << std::endl
             << "   [may be same as <input_file1>] containing the" << std::endl
             << "   second tree <sample2> to compared. Whatever" << std::endl
             << "   else they contain, these two trees must have" << std::endl
             << "   a single vector d[] of type double, length n+1" << std::endl
             << "   where n is the sample dimension." << std::endl
             << "       d[0:n] : sample n-dimensional vector data" << std::endl
             << "       d[n] : sample weight for this vector" << std::endl
             << "Options can be any of the following:" << std::endl
             << "   -f : use single precision floats to compute the" << std::endl
             << "        sample interaction energy, default double" << std::endl
             << "   -g : use a CUDA kernel to compute the sample" << std::endl
             << "        interaction energy on a gpu, implies -f" << std::endl
             << "   -s <seed> : set random number generator <seed>" << std::endl
             << "        when setting up the data set partitions," << std::endl
             << "        randomized by default, unsigned int value" << std::endl
             << "   -i <ilist> : include only the data components" << std::endl
             << "        listed in <ilist> when comparing the two" << std::endl
             << "        samples, where <ilist> contains a comma-" << std::endl
             << "        separated list of component indices i or" << std::endl
             << "        index ranges i-j, eg. -i 0,1,5-7,9,10" << std::endl
             << "        default is to include all components" << std::endl
             << "   -w <n>/><N> : split the work into <N> equal" << std::endl
             << "        parts and only do the work of part <n>," << std::endl
             << "        with <n> in [0,<N>), default is -w 0/1" << std::endl
             << "Partial results are saved in output ROOT file" << std::endl
             << "samplesep.root in a histogram named \"energy\"" << std::endl;
   exit(1);
}

void read_sample(const char *infile, const char *sample_name,
                 const std::vector<int> &include_list)
{
   TFile fin(infile);
   TTree *tsample = (TTree*)fin.Get(sample_name);
   if (tsample == NULL) {
      std::cerr << "samplesep.read_sample error - sample named " 
                << sample_name << " not found in " << infile
                << std::endl;
      exit(9);
   }
   double dsample[99999];
   TLeaf *lsample = tsample->GetLeaf("dsample", "d");
   int dlen = lsample->GetLenStatic();
   tsample->SetBranchAddress("dsample", dsample);
   if (sample.find(sample_name) != sample.end())
      sample[sample_name].clear();
   for (int i=0; tsample->GetEntry(i); ++i) {
      std::vector<double> d;
      for (int k=0; k < dlen - 1; ++k) {
         if (include_list.size() == 0)
            d.push_back(dsample[k]);
         else if (k < (int)include_list.size() && include_list[k] > 0)
            d.push_back(dsample[k]);
      }
      d.push_back(dsample[dlen - 1]);
      sample[sample_name].push_back(d);
   }
}

void compute_energy(const char *sample1, const char *sample2,
                    std::vector<double> &energy)
{
   if (sample.find(sample1) == sample.end())
   {
      std::cerr << "samplesep.compute_energy error - sample named " 
                << sample1 << " not found in memory"
                << std::endl;
      exit(9);
   }
   else if (sample.find(sample2) == sample.end()) {
      std::cerr << "samplesep.compute_energy error - sample named " 
                << sample2 << " not found in memory"
                << std::endl;
      exit(9);
   }
   if (strlen(sample1) == 0 || strlen(sample2) == 0)
      return;
   if (sample[sample1][0].size() != sample[sample2][0].size()) {
      std::cerr << "samplesep.compute_energy error - samples " 
                << sample1 << " and " << sample2
                << " cannot be compared, different dimensions "
                << sample[sample1][0].size() << " != "
                 << sample[sample2][0].size()
                << std::endl;
      exit(9);
   }
   int sample_dim = sample[sample1][0].size() - 1;
   if (sample_dim == 0) {
      std::cerr << "samplesep.compute_energy error - "
                << "zero sample size, nothing to do."
                << std::endl;
      return;
   }
 
   std::vector<std::vector<double> > dMoments;
   std::vector<double> z3 = {0,0,0};
   for (int k=0; k < sample_dim; ++k) {
      dMoments.push_back(z3);
   }
   for (int i=0; i < (int)sample[sample1].size(); ++i) {
      for (int k=0; k < sample_dim; ++k) {
         double a = sample[sample1][i][k];
         double w = sample[sample1][i].back();
         dMoments[k][0] += w;
         dMoments[k][1] += w * a;
         dMoments[k][2] += w * a * a;
      }
   }
   for (int i=0; i < (int)sample[sample2].size(); ++i) {
      for (int k=0; k < sample_dim; ++k) {
         double a = sample[sample2][i][k];
         double w = sample[sample2][i].back();
         dMoments[k][0] += w;
         dMoments[k][1] += w * a;
         dMoments[k][2] += w * a * a;
      }
   }
   std::vector<double> stdev;
   for (int k=0; k < sample_dim; ++k) {
      if (dMoments[k][0] == 0) {
         std::cout << "bad moments in component " 
                   << k << ":" << dMoments[k][0]
                   << std::endl;
         exit(9);
      }
      double mean = dMoments[k][1] / dMoments[k][0];
      double vari = dMoments[k][2] / dMoments[k][0] - mean * mean;
      if (vari == 0) {
         std::cout << "bad variance component " << k 
                   << ":" << "mean,vari=" << mean << "," << vari
                   << std::endl;
         exit(9);
      }
      stdev.push_back(sqrt(fabs(vari)));
   }

   std::vector<std::vector<double> > dsample;
   for (int i=0; i < (int)sample[sample1].size(); ++i) {
      std::vector<double> d;
      std::vector<unsigned int> p;
      for (int k=0; k < sample_dim; ++k) {
         d.push_back(sample[sample1][i][k] / stdev[k]);
      }
      d.push_back(sample[sample1][i].back());
      dsample.push_back(d);
   }
   for (int i=0; i < (int)sample[sample2].size(); ++i) {
      std::vector<double> d;
      std::vector<unsigned int> p;
      for (int k=0; k < sample_dim; ++k) {
         d.push_back(sample[sample2][i][k] / stdev[k]);
      }
      d.push_back(sample[sample2][i].back());
      dsample.push_back(d);
   }

   std::vector<std::vector<char> > partition(dsample.size());

   std::vector<double> aCharge(Npartitions);
   std::vector<double> bCharge(Npartitions);

   if (random_seed == 0)
      std::srand(std::time(0));
   else
      std::srand(random_seed);

   std::vector<int> p0;
   for (int i=0; i < (int)dsample.size(); ++i)
      p0.push_back(i);
   for (int n=0; n < (int)Npartitions; ++n) {
      std::vector<int> p = p0;
      if (n > 0)
         std::random_shuffle(p.begin(), p.end());
      for (int i=0; i < (int)dsample.size(); ++i) {
         if (p[i] < (int)sample[sample1].size()) {
            partition[i].push_back(0);
            aCharge[n] += dsample[i].back();
         }
         else {
            partition[i].push_back(1);
            bCharge[n] += dsample[i].back();
         }
      }
   }
   std::vector<double> aaNorm(Npartitions);
   std::vector<double> bbNorm(Npartitions);
   std::vector<double> abNorm(Npartitions);
   for (int n=0; n < (int)Npartitions; ++n) {
      aaNorm[n] = 1 / (aCharge[n] * aCharge[n]);
      bbNorm[n] = 1 / (bCharge[n] * bCharge[n]);
      abNorm[n] = 1 / (aCharge[n] * bCharge[n]);
   }

   std::vector<double> Ewsum(Npartitions);
   double *Ewsum_ = Ewsum.data();

   //            T H E   I N N E R   L O O P
   //
   // This is the inner loop of the energy sum calculation,
   // which sums over all pairs i<j in the joint sample.
   // The work is divided into a N x N grid of blocks of
   // dimension M x M covering the upper-diagonal triangle
   // of the matrix. Blocks are numbered sequentially
   // starting from zero in the upper-left corner of the
   // matrix, going left-to-right across the matrix until
   // reaching the last column, then jumping to the diagonal
   // block on the next row and continuing until all 
   // n(n+1)/2 blocks have been covered.

   const int M(1024);
   const int N = (dsample.size() + M - 1) / M;
   int block_count = N * (N + 1) / 2;
   int block_start = block_count * work_set / work_sets;
   int block_end = block_count * (work_set + 1) / work_sets;
   for (int block = block_start; block < block_end; ++block) {
      int block_idy = (N + 0.5) - sqrt(sqr(N + 0.5) - 2 * block);
      int block_idx = block - (N - 1) * block_idy + 
                      block_idy * (block_idy - 1) / 2;
      int i0 = block_idy * M;
      int i1 = (i0 + M < (int)dsample.size())? i0 + M : dsample.size();

#if USE_OPENMP_MULTITHREADING
  #pragma omp parallel for reduction(+:Ewsum_[:1024])
#endif
      for (int i = i0; i < i1; ++i) {
         int j0 = block_idx * M;
         int j1 = (j0 + M < (int)dsample.size())? j0 + M : dsample.size();
         j0 = (j0 > i)? j0 : i + 1;
         for (int j = j0; j < j1; ++j) {
            double r2(1e-4);
            for (int k=0; k < sample_dim; ++k)
               r2 += sqr(dsample[i][k] - dsample[j][k]);
            double E = 1 / r2;
            double Ew = E * dsample[i].back() * dsample[j].back();
            for (int n=0; n < (int)Npartitions; ++n) {
               if (partition[i][n] == 0 && partition[j][n] == 0)
                  Ewsum_[n] += Ew * aaNorm[n];
               else if (partition[i][n] == 1 && partition[j][n] == 1)
                  Ewsum_[n] += Ew * bbNorm[n];
               else
                  Ewsum_[n] -= Ew * abNorm[n];
            }
         }
      }

      std::cout << "completed block " << block 
                << "(" << block_idx << "," << block_idy << ")"
                << "\r" << std::flush;
   }
   for (int n=0; n < (int)Npartitions; ++n) {
      energy.push_back(Ewsum[n]);
   }
   std::cout << std::endl;
}

void compute_energy(const char *sample1, const char *sample2,
                    std::vector<float> &energy)
{
   if (sample.find(sample1) == sample.end())
   {
      std::cerr << "samplesep.compute_energy error - sample named " 
                << sample1 << " not found in memory"
                << std::endl;
      exit(9);
   }
   else if (sample.find(sample2) == sample.end()) {
      std::cerr << "samplesep.compute_energy error - sample named " 
                << sample2 << " not found in memory"
                << std::endl;
      exit(9);
   }
   if (strlen(sample1) == 0 || strlen(sample2) == 0)
      return;
   if (sample[sample1][0].size() != sample[sample2][0].size()) {
      std::cerr << "samplesep.compute_energy error - samples " 
                << sample1 << " and " << sample2
                << " cannot be compared, different dimensions "
                << sample[sample1][0].size() << " != "
                 << sample[sample2][0].size()
                << std::endl;
      exit(9);
   }
   int sample_dim = sample[sample1][0].size() - 1;
   if (sample_dim == 0) {
      std::cerr << "samplesep.compute_energy error - "
                << "zero sample size, nothing to do."
                << std::endl;
      return;
   }
 
   std::vector<std::vector<double> > dMoments;
   std::vector<double> z3 = {0,0,0};
   for (int k=0; k < sample_dim; ++k) {
      dMoments.push_back(z3);
   }
   for (int i=0; i < (int)sample[sample1].size(); ++i) {
      for (int k=0; k < sample_dim; ++k) {
         double a = sample[sample1][i][k];
         double w = sample[sample1][i].back();
         dMoments[k][0] += w;
         dMoments[k][1] += w * a;
         dMoments[k][2] += w * a * a;
      }
   }
   for (int i=0; i < (int)sample[sample2].size(); ++i) {
      for (int k=0; k < sample_dim; ++k) {
         double a = sample[sample2][i][k];
         double w = sample[sample2][i].back();
         dMoments[k][0] += w;
         dMoments[k][1] += w * a;
         dMoments[k][2] += w * a * a;
      }
   }
   std::vector<double> stdev;
   for (int k=0; k < sample_dim; ++k) {
      if (dMoments[k][0] == 0) {
         std::cout << "bad moments in component " 
                   << k << ":" << dMoments[k][0]
                   << std::endl;
         exit(9);
      }
      double mean = dMoments[k][1] / dMoments[k][0];
      double vari = dMoments[k][2] / dMoments[k][0] - mean * mean;
      if (vari == 0) {
         std::cout << "bad variance component " << k 
                   << ":" << "mean,vari=" << mean << "," << vari
                   << std::endl;
         exit(9);
      }
      stdev.push_back(sqrt(fabs(vari)));
   }

   std::vector<std::vector<float> > fsample;
   for (int i=0; i < (int)sample[sample1].size(); ++i) {
      std::vector<float> d;
      std::vector<unsigned int> p;
      for (int k=0; k < sample_dim; ++k) {
         d.push_back(sample[sample1][i][k] / stdev[k]);
      }
      d.push_back(sample[sample1][i].back());
      fsample.push_back(d);
   }
   for (int i=0; i < (int)sample[sample2].size(); ++i) {
      std::vector<float> d;
      std::vector<unsigned int> p;
      for (int k=0; k < sample_dim; ++k) {
         d.push_back(sample[sample2][i][k] / stdev[k]);
      }
      d.push_back(sample[sample2][i].back());
      fsample.push_back(d);
   }

   std::vector<std::vector<char> > partition(fsample.size());

   std::vector<double> aCharge(Npartitions);
   std::vector<double> bCharge(Npartitions);

   if (random_seed == 0)
      std::srand(std::time(0));
   else
      std::srand(random_seed);

   std::vector<int> p0;
   for (int i=0; i < (int)fsample.size(); ++i)
      p0.push_back(i);
   for (int n=0; n < (int)Npartitions; ++n) {
      std::vector<int> p = p0;
      if (n > 0)
         std::random_shuffle(p.begin(), p.end());
      for (int i=0; i < (int)fsample.size(); ++i) {
         if (p[i] < (int)sample[sample1].size()) {
            partition[i].push_back(0);
            aCharge[n] += fsample[i].back();
         }
         else {
            partition[i].push_back(1);
            bCharge[n] += fsample[i].back();
         }
      }
   }
   std::vector<float> aaNorm(Npartitions);
   std::vector<float> bbNorm(Npartitions);
   std::vector<float> abNorm(Npartitions);
   for (int n=0; n < (int)Npartitions; ++n) {
      aaNorm[n] = 1 / (aCharge[n] * aCharge[n]);
      bbNorm[n] = 1 / (bCharge[n] * bCharge[n]);
      abNorm[n] = 1 / (aCharge[n] * bCharge[n]);
   }

   std::vector<float> Ewsum(Npartitions);
   float *Ewsum_ = Ewsum.data();

#if USE_OPENMP_MULTITHREADING
  #pragma omp parallel for reduction(+:Ewsum_[:1024])
#endif
   for (int i=0; i < (int)fsample.size(); ++i) {
      for (int j=i+1; j < (int)fsample.size(); ++j) {
         float r2(1e-4);
         for (int k=0; k < sample_dim; ++k)
            r2 += sqr(fsample[i][k] - fsample[j][k]);
         float E = 1 / r2;
         float Ew = E * fsample[i].back() * fsample[j].back();
         for (int n=0; n < (int)Npartitions; ++n) {
            if (partition[i][n] == 0 && partition[j][n] == 0)
               Ewsum_[n] += Ew * aaNorm[n];
            else if (partition[i][n] == 1 && partition[j][n] == 1)
               Ewsum_[n] += Ew * bbNorm[n];
            else
               Ewsum_[n] -= Ew * abNorm[n];
         }
      }
      std::cout << "outer loop index i=" << i << "\r" << std::flush;
   }
   for (int n=0; n < (int)Npartitions; ++n) {
      energy.push_back(Ewsum[n]);
   }
}

#ifdef __CUDACC__
void compute_energy_gpu(const char *sample1, const char *sample2,
                        std::vector<float> &energy)
{
}

#else
void compute_energy_gpu(const char *sample1, const char *sample2,
                        std::vector<float> &energy)
{
   std::cout << "gpu not implemented" << std::endl;
}
#endif

void parse_ilist(const char *ilist, std::vector<int> &include_list)
{
   std::stringstream s0(ilist);
   std::string irange;
   while (getline(s0, irange, ',')) {
      std::stringstream s1(irange);
      std::string index;
      std::vector<int> limits;
      for (int i=0; getline(s1, index, '-'); ++i)
         limits.push_back(atoi(index.c_str()));
      if (limits.size() == 0)
         continue;
      for (int i=include_list.size(); i <= limits.front(); ++i)
         include_list.push_back(0);
      for (int i=limits.front(); i <= limits.back(); ++i) {
         if (i >= (int)include_list.size())
            include_list.push_back(1);
         else
            include_list[i] = 1;
      }
   }
}

int main(int argc, char *argv[]) {
   int singleprec(0);
   int usegpu(0);
   int nfile(0);
   char *infile[2] = {0,0};
   char *sname[2] = {0,0};
   std::vector<int> include_list;
   for (int i=1; i < argc; ++i) {
      if (strstr(argv[i], "-f") == argv[i])
         singleprec = 1;
      else if (strstr(argv[i], "-c") == argv[i])
         usegpu = 1;
      else if (strstr(argv[i], "-i") == argv[i]) {
         if (strlen(argv[i]) > 2)
            parse_ilist(argv[i] + 2, include_list);
         else
            parse_ilist(argv[++i], include_list);
      }
      else if (strstr(argv[i], "-s") == argv[i]) {
         if (strlen(argv[i]) > 2)
            random_seed = strtoul(argv[i] + 2, 0, 0);
         else
            random_seed = strtoul(argv[++i], 0, 0);
      }
      else if (strstr(argv[i], "-w") == argv[i]) {
         char *wset;
         if (strlen(argv[i]) > 2)
            wset = argv[i] + 2;
         else
            wset = argv[++i];
         char *wsets = strstr(wset, "/");
         if (wsets == 0)
            usage();
         *(wsets++) = 0;
         work_set = strtol(wset, 0, 0);
         work_sets = strtol(wsets, 0, 0);
      }
      else if (argv[i][0] == '-')
         usage();
      else
         infile[nfile++] = argv[i];
   }
   if (nfile != 2)
      usage();
   sname[0] = strstr(infile[0], ":");
   sname[1] = strstr(infile[1], ":");
   if (sname[0] == 0 || sname[1] == 0)
      usage();
   else if (sname[0] == infile[0])
      usage();
   else if (sname[1] == infile[1])
      infile[1] = infile[0];
   *(sname[0]++) = 0;
   *(sname[1]++) = 0;

   read_sample(infile[0], sname[0], include_list);
   read_sample(infile[1], sname[1], include_list);

   std::vector<float> fenergy;
   std::vector<double> denergy;

   if (usegpu) {
      compute_energy_gpu(sname[0], sname[1], fenergy);
      if (fenergy.size() > 0) {
         double moment[3] = {0,0,0};
         for (int i=1; i < (int)fenergy.size(); ++i) {
            moment[0] += 1;
            moment[1] += fenergy[i];
            moment[2] += fenergy[i] * fenergy[i];
         }
         double mean = moment[1] / moment[0];
         double vari = moment[2] / moment[0] - mean * mean;
         std::cout << "got back energy " << fenergy[0] << ", "
                   << "mean=" << mean << ", sigma=" << sqrt(vari) 
                   << std::endl;
         std::cout << "significance is " << (fenergy[0] - mean) / sqrt(vari)
                   << " sigma" << std::endl;
      }
   }
   else if (singleprec) {
      compute_energy(sname[0], sname[1], fenergy);
      if (fenergy.size() > 0) {
         double moment[3] = {0,0,0};
         for (int i=1; i < (int)fenergy.size(); ++i) {
            moment[0] += 1;
            moment[1] += fenergy[i];
            moment[2] += fenergy[i] * fenergy[i];
         }
         double mean = moment[1] / moment[0];
         double vari = moment[2] / moment[0] - mean * mean;
         std::cout << "got back energy " << fenergy[0] << ", "
                   << "mean=" << mean << ", sigma=" << sqrt(vari) 
                   << std::endl;
         std::cout << "significance is " << (fenergy[0] - mean) / sqrt(vari)
                   << " sigma" << std::endl;
      }
   }
   else {
      compute_energy(sname[0], sname[1], denergy);
      if (denergy.size() > 0) {
         double moment[3] = {0,0,0};
         for (int i=1; i < (int)denergy.size(); ++i) {
            moment[0] += 1;
            moment[1] += denergy[i];
            moment[2] += denergy[i] * denergy[i];
         }
         double mean = moment[1] / moment[0];
         double vari = moment[2] / moment[0] - mean * mean;
         std::cout << "got back energy " << denergy[0] << ", "
                   << "mean=" << mean << ", sigma=" << sqrt(vari) 
                   << std::endl;
         std::cout << "significance is " << (denergy[0] - mean) / sqrt(vari)
                   << " sigma" << std::endl;
      }
   }

   int nenergy = denergy.size() + fenergy.size();
   if (nenergy > 0) {
      TFile fout("samplesep.root", "update");
      std::stringstream title;
      title << "energy of interaction between samples " 
            << sname[0] << " and " << sname[1];
      if (work_set != 0 or work_sets != 1)
         title << " part " << work_set << "/" << work_sets;
      TH1D hen("energy", title.str().c_str(), nenergy, 0, nenergy);
      for (int i=0; i < (int)fenergy.size(); ++i)
         hen.SetBinContent(i+1, fenergy[i]);
      for (int i=0; i < (int)denergy.size(); ++i)
         hen.SetBinContent(i+1, denergy[i]);
      hen.Write();
   }
   return 0;
}
#endif
