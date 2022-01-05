//
// beamtree.cc : main program for GlueXBremsstrahlung beam generation
//
// author: richard.t.jones at uconn.edu
// version: february 4, 2017
//

#include <GlueXBremsstrahlungGenerator.hh>
#include <unistd.h>
#include <iostream>
#include <TFile.h>

void usage()
{
   std::cout << std::endl
          << "Usage: beamtree [options]" << std::endl
          << " where options are:" << std::endl
          << "   -n # : generate # events, default is 1000" << std::endl
          << "   -r # : set the random seed to #, default is 1" << std::endl
          << "   -o <outfile.root> : write to outfile.root" << std::endl
          << "                       (default is beamtree.root)"
          << std::endl;
   exit(9);
}

int main(int argc,char** argv)
{
   int nevents = 1000;
   long int seed = 1;
   TFile *outfile = 0;

   char c;
   while ((c = getopt(argc, argv, "r:n:o:")) != -1) {
      if (c == 'r') {
         seed = atoi(optarg);
      }
      else if (c == 'n') {
         nevents = atoi(optarg);
      }
      else if (c == 'o') {
         outfile = new TFile(optarg, "recreate");
      }
      else {
         usage();
      }
   }

   GlueXBremsstrahlungGenerator gen(outfile);
   gen.fRandom.SetSeed(seed);
   gen.GenerateBeamPhotons(nevents);
   return 0;
}
