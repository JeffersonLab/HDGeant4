//
// AdaptiveSampler - provide an adaptive algorithm for improved 
//                   importance sampling in Monte Carlo integration
//
// file: AdaptiveSampler.hh (see AdaptiveSampler.cc for implementation)
// author: richard.t.jones at uconn.edu
// original release: february 23, 2016
//
// new feature release: may 30,2 2020
//    - see new note 4 below
//    - richard.t.jones at uconn.edu
//
// usage:
// 1) Express your Monte Carlo integral in the form
//       result = Int[0,1]^N { integrand(u[N]) d^N u }
//    where you have included the Jacobian that relates the measure
//    of your problem's integration coordinates to these rescaled
//    coordinates u[i], i=1..N that vary from 0 to 1. In this form,
//    a simple MC estimator for result is given by
//      result_MC = (1 / nMC) sum[i=1,nMC] { integrand(u_i) }
//    where the vectors u_i in the above sum are uniformly sampled
//    within the unit hypercube [0,1]^N. The idea of importance
//    sampling stems from the fact that the statistical noise on this
//    sum can be extremely large, even for very large nMC, if the
//    integrand is sharply peaked around specific points in the unit
//    hypercube. Convergence can be greatly improved by replacing the
//    result_MC estimator with the importance-sampled estimator:
//      result_IS = (1 / nMC) sum[i=1,nMC] { integrand(u_i) w_i }
//    where the u_i are non-uniformly sampled over the unit hypercube
//    and the weights w_i compensate for the non-uniformity by
//    multiplying the integrand by 1 / the_oversampling_factor.
//    By preferentially visiting the regions in the hypercube where
//    the integrand is large, evaluation of the integral to a given
//    statistical precision can be accomplished with a smaller nMC.
//    In the large nMC limit, it can be proved that the statistical
//    variance on this estimator is
//      var_IS = (1 / nMC^2) sum[i=1,nMC] { integrand(u_i) w_i }^2
//
// 2) The efficiency of the MC estimation procedure is measured by the
//    number of events that would be needed to obtain an equivalent
//    statistical precision over any small subset of the full domain
//    of integration for a perfect sampler, divided by the actual
//    number required by this sampler. A perfect sampler for a given
//    integrand is one for which the product { integrand(u_i) w_i } is
//    the same for all i. For a perfect sampler, var_IS is exactly 0
//    because every sample returns the same value: the true result.
//    Actually achieving a perfect sampler is not a practical goal
//    in a problem solving because it requires one to already know the
//    answer in order to compute it. However, getting reasonably close
//    to unity is all that is needed to ensure rapid convergence.
//    In the large N limit, the efficiency can be measured as
//      efficiency = (S_1 S_1) / (S_0 S_2)
//    where
//      S_P = sum[i=1,nMC] { integrand(u_i) w_i }^P
//    Thus the above formula for var_IS should thus be augmented by the
//    factor sqrt(1 - efficiency), although this makes very little
//    difference in most practical situations where efficiencies close
//    to unity are very difficult to achieve.
//
// 3) The AdaptiveSampler class generates a sequence of (u,w) pairs
//    which a user application uses to compute the integrand(u_i) whose
//    sum gives result_IS. If the user feeds back the integrand(u_i) to
//    AdaptiveSampler, the AdaptiveSampler can use this information to
//    internally adjust its weighting scheme for the following (u,w)
//    pairs that it generates, thus improving the rate at which var_IS
//    decreases with increasing nMC. In the following example, D is the
//    dimension of the Monte Carlo integral, ie. the number of independent
//    u variables ~Unif(0,1) are needed to completely specify a point in
//    the space of integration.
//
//    #include <TRandom.h>
//    #include <AdaptiveSampler.h>
//
//    // This is just an example,
//    // implement my_randoms as you wish.
//    TRandom randoms(0);
//    void my_randoms(int n, double *u) { randoms.RndmArray(n,u); }
//
//    AdaptiveSampler sampler(D, my_randoms);
//    double sum[3] = {0,0,0};
//    for (int i=0; i < nMC; i++) {
//       double u[D];
//       double w = sampler.sample(u);
//       double I = ... compute your integrand here ...
//       double wI = w*I;
//       sum[0] += 1;
//       sum[1] += wI;
//       sum[2] += wI*wI;
//       sampler.feedback(u,wI); // if you want AdaptiveSampler to adapt
//       if (i % 1000 == 0)
//          sampler.adapt();     // adapt based on feedback received
//    }
//    double mu = sum[1] / sum[0];
//    double effic = sum[1] * sum[1] / (sum[0] * sum[2]);
//    double sigma = sqrt((1 - effic) * sum[2]) / sum[0];
//    std::cout << "result_IS is " << mu << " +/- " << sigma << std::endl;
//
// 4) NEW FEATURE: Support for Parametric Models
//    Prior to this introduction of this new feature, it was assumed that the
//    user's integrand function wI(u) is a deterministic function of the random
//    variables u[D]. However this is sometimes not the case in Monte Carlo
//    integration, where some parameters for a process might be read in from
//    an external source, and the function wI(u) might change from one event to
//    the next, depending on the values of those parameters. This can ruin the
//    convergence the adaptation algorithm, as the variation of the parameters
//    will be wrongly interpreted as variation of wI(u) with u, leading to cell
//    adaptation that produces no improvement in the efficiency. To address this
//    situation, I now allow the user to specify a dimension D of the hypercube
//    that is larger than the u-space being sampled by sample(), with the extra
//    dimensions taken up by the parameters. The parameters must be mapped onto
//    the interval [0,1] and supplied as input to sample() each time it is called.
//    If parameters are present, they are stored in the first P elements of u[D]
//    when it is passed to sample(), with their count given by nfixed, the new
//    last argument to sample(). The default value of nfixed is 0, so code that
//    is based on earlier versions of AdaptiveSampler will continue to compile
//    and run as before, without modification.
//
//    In summary, if sample() is called with argument nfixed > 0, then the first
//    nfixed elements of u[D] will not be modified, and only the last N - nfixed
//    elements will be overwritten with the new hypercube vector for this event.
//    The user-supplied parameters are treated on the same basis as the variables
//    that sample() generates: as Monte Carlo integration variables that together
//    define the integral being estimated. This effectively means that the result
//    for the integral is an average over the parameters, with the weight in the
//    parametric average determined externally by the user's parameter generator.
//    

#include <vector>
#include <string>
#include <map>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <assert.h>

using Uniform01 = void (*)(int n, double *randoms);

class AdaptiveSampler {
 public:
   AdaptiveSampler(int dim, Uniform01 user_generator, int nfixed=0);
   AdaptiveSampler(const AdaptiveSampler &src);
   ~AdaptiveSampler();
   AdaptiveSampler operator=(const AdaptiveSampler &src);

 public:
   // action methods
   double sample(double *u);
   void feedback(const double *u, double wI);
   void optimize_tree();
   int adapt();
   void reset_stats();

   // getter methods
   int getNdim() const;
   int getNfixed() const;
   int getNcells() const;
   long int getNsample() const;
   double getWItotal() const;
   double getWI2total(bool optimized=false);
   double getWI4total(bool optimized=false);
   double getEfficiency(bool optimized=false);
   double getResult(double *error=0, double *error_uncertainty=0);
   double getReweighted(double *error=0, double *error_uncertainty=0);
   double getAdaptation_sampling_threshold() const;
   double getAdaptation_efficiency_target() const;
   int getAdaptation_maximum_depth() const;
   int getAdaptation_maximum_cells() const;
   static int getVerbosity();
   
   // setter methods
   void setAdaptation_sampling_threshold(double threshold);
   void setAdaptation_efficiency_target(double target);
   void setAdaptation_maximum_depth(int depth);
   void setAdaptation_maximum_cells(int ncells);
   static void setVerbosity(int verbose);

   // persistency methods
   int saveState(const std::string filename, bool optimized=false) const;
   int mergeState(const std::string filename);
   int restoreState(const std::string filename);

   // diagnostic methods
   void display_tree(bool optimized=false);
   int check_subsets(bool optimized=false);

 protected:
   int fNdim;
   int fNfixed;
   Uniform01 fRandom;
   unsigned int fMaximum_depth;
   unsigned int fMaximum_cells;
   double fSampling_threshold;
   double fEfficiency_target;
   double fMinimum_sum_wI2_delta;
   static int verbosity;

   class Cell;
   Cell *fTopCell;

   Cell *findCell(double ucell, int &depth, 
                  double *u0, double *u1,
                  const double *u);
   Cell *findCell(const double *u, int &depth, 
                  double *u0, double *u1) const;
   int recursively_update(std::vector<int> index);
   double display_tree(Cell *cell, double subset, std::string id,
                       double *u0, double *u1, bool optimized=false);

   // internal weighting tables

   class Cell {
    public:
      int ndim;
      int nfixed;
      int divAxis;
      long int nhit;
      double sum_wI;
      double sum_wI2;
      double sum_wI4;
      double sum_wI2s;
      double *sum_wI2d;
      double *sum_wI4d;
      double subset;
      // optimized transforms of the above statistics
      long int opt_nhit;
      // sum_wI is an invariant;
      double opt_wI2;
      double opt_wI4;
      // sum_wI2s is an invariant;
      // sum_wI2d is input for optimization;
      // sum_wI4d is input for optimization;
      double opt_subset;

      Cell *subcell[3];

      Cell(int dim, int fixed=0)
       : ndim(dim),
         nfixed(fixed),
         divAxis(-1),
         nhit(0),
         sum_wI(0),
         sum_wI2(0),
         sum_wI4(0),
         sum_wI2s(0),
         subset(0),
         opt_nhit(0),
         opt_wI2(0), 
         opt_wI4(0),
         opt_subset(0)
      {
         sum_wI2d = new double[3*ndim];
         sum_wI4d = new double[3*ndim];
         std::fill(sum_wI2d, sum_wI2d + 3*ndim, 0);
         std::fill(sum_wI4d, sum_wI4d + 3*ndim, 0);
         subcell[0] = 0;
         subcell[1] = 0;
         subcell[2] = 0;
      }
      Cell(const Cell &src) {
         ndim = src.ndim;
         nfixed = src.nfixed;
         divAxis = src.divAxis;
         nhit = src.nhit;
         sum_wI = src.sum_wI;
         sum_wI2 = src.sum_wI2;
         sum_wI4 = src.sum_wI4;
         sum_wI2s = src.sum_wI2s;
         sum_wI2d = new double[3*ndim];
         sum_wI4d = new double[3*ndim];
         for (int i=0; i < 3*ndim; ++i) {
            sum_wI2d[i] = src.sum_wI2d[i];
            sum_wI4d[i] = src.sum_wI4d[i];
         }
         subset = src.subset;
         opt_nhit = src.opt_nhit;
         opt_wI2 = src.opt_wI2;
         opt_wI4 = src.opt_wI4;
         opt_subset = src.opt_subset;
         if (src.subcell[0] != 0) {
            subcell[0] = new Cell(*src.subcell[0]);
            subcell[1] = new Cell(*src.subcell[1]);
            subcell[2] = new Cell(*src.subcell[2]);
         }
         else {
            subcell[0] = 0;
            subcell[1] = 0;
            subcell[2] = 0;
         }
      }
      ~Cell() {
         if (subcell[0] != 0) {
            delete subcell[0];
            delete subcell[1];
            delete subcell[2];
         }
         delete [] sum_wI2d;
         delete [] sum_wI4d;
      }
      Cell operator=(const Cell &src) {
         return Cell(src);
      }
      void reset_stats() {
         nhit = 0;
         sum_wI = 0;
         sum_wI2 = 0;
         sum_wI4 = 0;
         sum_wI2s = 0;
         std::fill(sum_wI2d, sum_wI2d + 3*ndim, 0);
         std::fill(sum_wI4d, sum_wI4d + 3*ndim, 0);
         if (divAxis > -1) {
            subcell[0]->reset_stats();
            subcell[1]->reset_stats();
            subcell[2]->reset_stats();
         }
      }
      int sum_stats()
      {
         int count = 1;
         if (divAxis > -1) {
            nhit = 0;
            sum_wI = 0;
            sum_wI2 = 0;
            sum_wI4 = 0;
            sum_wI2s = 0;
            for (int n=0; n < 3; ++n) {
               count += subcell[n]->sum_stats();
               nhit += subcell[n]->nhit;
               sum_wI += subcell[n]->sum_wI;
               sum_wI2 += subcell[n]->sum_wI2;
               sum_wI4 += subcell[n]->sum_wI4;
               sum_wI2s += subcell[n]->sum_wI2s;
            }
         }
         if (divAxis < nfixed) {
            sum_wI2s = sqrt(nhit * sum_wI2);
         }
         return count;
      }
      void optimize() {
         // Assume opt_nhit and opt_subset already set upon entry,
         // task is to assign opt_wI2, opt_wI4 for this cell, and
         // all opt_* parameters for child nodes in the tree.
         if (divAxis > -1) {
            if (divAxis < nfixed) {
               opt_wI2 = 0;
               opt_wI4 = 0;
               for (int n=0; n < 3; ++n) {
                  subcell[n]->opt_nhit = subcell[n]->nhit;
                  subcell[n]->opt_subset = opt_subset;
                  subcell[n]->optimize();
                  opt_wI2 += subcell[n]->opt_wI2;
                  opt_wI4 += subcell[n]->opt_wI4;
               }
            }
            else {
               opt_wI2 = 0;
               opt_wI4 = 0;
               for (int n=0; n < 3; ++n) {
                  double f0 = 5e-3;  // minimum subset fraction
                  double r = (subcell[n]->sum_wI2s / sum_wI2s + f0)
                                                   / (1 + 3*f0);
                  r = (nhit > 100)? r : 1;
                  subcell[n]->opt_nhit = nhit * r;
                  subcell[n]->opt_subset = opt_subset * r;
                  subcell[n]->optimize();
                  opt_wI2 += subcell[n]->opt_wI2;
                  opt_wI4 += subcell[n]->opt_wI4;
               }
            }
         }
         else if (nhit > 0) {
            double r = (opt_nhit + 1e-99) / nhit;
            opt_wI2 = sum_wI2 / r;
            opt_wI4 = sum_wI4 / pow(r,3);
         }
      }
      int serialize(std::ofstream &ofs, bool optimized=false) {
         ofs << "divAxis=" << divAxis << std::endl;
         if (nhit != 0)
            ofs << "nhit=" << nhit << std::endl;
         if (sum_wI != 0)
            ofs << "sum_wI=" << sum_wI << std::endl;
         if (sum_wI2 != 0)
            ofs << "sum_wI2=" << sum_wI2 << std::endl;
         if (sum_wI4 != 0)
            ofs << "sum_wI4=" << sum_wI4 << std::endl;
         for (int i=0; i < 3*ndim; ++i) {
            if (sum_wI2d[i] != 0)
               ofs << "sum_wI2d[" << i << "]=" << sum_wI2d[i] << std::endl;
         }
         for (int i=0; i < 3*ndim; ++i) {
            if (sum_wI4d[i] != 0)
               ofs << "sum_wI4d[" << i << "]=" << sum_wI4d[i] << std::endl;
         }
         if (optimized)
            ofs << "subset=" << std::setprecision(20) << opt_subset << std::endl;
         else
            ofs << "subset=" << std::setprecision(20) << subset << std::endl;
         ofs << "=" << std::endl;
         int count = 1;
         if (divAxis > -1) {
            count += subcell[0]->serialize(ofs, optimized);
            count += subcell[1]->serialize(ofs, optimized);
            count += subcell[2]->serialize(ofs, optimized);
         }
         return count;
      }
      int deserialize(std::ifstream &ifs, double subset_multiplier=1)
      {
         std::map<std::string,double> keyval;
         while (true) {
            std::string key;
            std::getline(ifs, key, '=');
            if (key.size() == 0) {
               std::getline(ifs, key);
               break;
            }
            ifs >> keyval[key];
            std::getline(ifs, key);
         }
         divAxis = keyval.at("divAxis");
         nhit += keyval["nhit"];
         sum_wI += keyval["sum_wI"];
         sum_wI2 += keyval["sum_wI2"];
         sum_wI4 += keyval["sum_wI4"];
         subset = keyval.at("subset") * subset_multiplier;
         std::map<std::string,double>::iterator it;
         for (it = keyval.begin(); it != keyval.end(); ++it) {
            if (it->first.substr(0,8) == "sum_wI2u") {
               int n3dim=1;
               for (int n=0; n < ndim; ++n)
                  n3dim *= 3;
               int j;
               if (sscanf(it->first.c_str(), "sum_wI2u[%d]=", &j) == 1) {
                  assert (j < n3dim);
                  for (int n=0; n < ndim; ++n) {
                     double wI2 = it->second;
                     sum_wI2d[n*3+(j%3)] += wI2;
                     sum_wI4d[n*3+(j%3)] += wI2*wI2 * n3dim/nhit;
                     j /= 3;
                  }
               }
            }
            else if (it->first.substr(0,8) == "sum_wI2d") {
               int j;
               if (sscanf(it->first.c_str(), "sum_wI2d[%d]=", &j) == 1) {
                  assert (j < 3*ndim);
                  sum_wI2d[j] += it->second;
               }
            }
            else if (it->first.substr(0,8) == "sum_wI4d") {
               int j;
               if (sscanf(it->first.c_str(), "sum_wI4d[%d]=", &j) == 1) {
                  assert (j < 3*ndim);
                  sum_wI4d[j] += it->second;
               }
            }
         }
         int count = 1;
         if (divAxis > -1) {
            //subset_multiplier *= (divAxis < 1)? 3 : 1;
            for (int i=0; i < 3; ++i) {
               if (subcell[i] == 0) {
                  subcell[i] = new Cell(ndim, nfixed);
               }
               count += subcell[i]->deserialize(ifs, subset_multiplier);
            }
         }
         return count;
      }
      int check_subsets(bool optimized=false, std::string id="0") {
         if (divAxis < 0) {
            return 0;
         }

         // test trinomial statistics at each node
         int warnings = 0;
         long int Ntot = subcell[0]->nhit + subcell[1]->nhit + subcell[2]->nhit;
         if (Ntot != nhit) {
            std::cerr << "Error in Cell::check_subsets - "
                      << "nhit consistency check #1 failed for cell " << id
                      << std::endl
                      << "  split cell nhit=" << nhit
                      << ", sum of subcell nhit=" << Ntot
                      << std::endl;
            return 999;
         }
         if (divAxis < nfixed) {
            if (subcell[0]->subset != subset ||
                subcell[1]->subset != subset ||
                subcell[2]->subset != subset)
            {
               std::cerr << "Error in Cell::check_subsets - "
                         << "subset consistency check #1 failed for cell " << id
                         << std::endl
                         << "  split fixed cell subset=" << subset
                         << ", subcell subsets="
                         << subcell[0]->subset << ","
                         << subcell[1]->subset << ","
                         << subcell[2]->subset
                         << std::endl;
               return 999;
            }
         }
         else {
            double subsetsum = subcell[0]->subset + 
                               subcell[1]->subset +
                               subcell[2]->subset;
            if (fabs(subsetsum - subset) > subset * 1e-12) {
               std::cerr << "Error in Cell::check_subsets - "
                         << "subset consistency check #2 failed for cell " << id
                         << std::endl
                         << "  split cell subset=" << subset
                         << ", sum of subcell subsets=" << subsetsum
                         << std::endl;
               return 999;
            }
            for (int i=0; i < 3; ++i) {
               double p = subcell[i]->subset / (subset + 1e-99);
               double mu = nhit * p;
               double sigma = sqrt(nhit * p * (1-p));
               double p1 = subcell[(i+1)%3]->subset / (subset + 1e-99);
               double mu1 = nhit * p1;
               double sigma1 = sqrt(nhit * p1 * (1-p1));
               double p2 = subcell[(i+2)%3]->subset / (subset + 1e-99);
               double mu2 = nhit * p2;
               double sigma2 = sqrt(nhit * p2 * (1-p2));
               if (subcell[i]->nhit < 30) {
                  double prob = 1;
                  for (int k=1; k <= nhit; ++k) {
                     if (k <=subcell[i]->nhit)
                        prob *= p * double(nhit - k + 1) / k;
                     else
                        prob *= 1 - p;
                  }
                  if (prob < 1e-6) {
                     std::cerr << "Warning in Cell::check_subsets - "
                               << "nhit - subset probability < 1e-6, cell "
                               << id << ", subcell " << i 
                               << " has nhit=" << subcell[i]->nhit
                               << ", expected " << mu << " +/- " << sigma
                               << ", P-value " << prob
                               << std::endl;
                     std::cerr << "  This cell is splits " << nhit
                               << " events along axis " << divAxis
                               << std::endl;
                     std::cerr << "  Other branch of this cell:"
                               << "  subcell " << (i+1)%3
                               << " with nhit=" << subcell[(i+1)%3]->nhit
                               << ", expected " << mu1 << " +/- " << sigma1
                               << ", " << (subcell[(i+1)%3]->nhit - mu1) / sigma1
                               << " sigma." << std::endl;
                     std::cerr << "  Other branch of this cell:"
                               << "  subcell " << (i+2)%3
                               << " with nhit=" << subcell[(i+2)%3]->nhit
                               << ", expected " << mu2 << " +/- " << sigma2
                               << ", " << (subcell[(i+2)%3]->nhit - mu2) / sigma2
                               << " sigma." << std::endl;
                     warnings++;
                  }
               }
               else if (fabs(subcell[i]->nhit - mu) > 5 * sigma) {
                  std::cerr << "Warning in Cell::check_subsets - "
                            << "nhit - subset mismatch > 5 sigma, cell "
                            << id << ", subcell " << i 
                            << " has nhit=" << subcell[i]->nhit
                            << ", expected " << mu << " +/- " << sigma
                            << ", " << (subcell[i]->nhit - mu) / sigma
                            << " sigma!" << std::endl;
                  std::cerr << "  This cell is splits " << nhit
                            << " events along axis " << divAxis
                            << std::endl;
                  std::cerr << "  Other branch of this cell:"
                            << "  subcell " << (i+1)%3
                            << " with nhit=" << subcell[(i+1)%3]->nhit
                            << ", expected " << mu1 << " +/- " << sigma1
                            << ", " << (subcell[(i+1)%3]->nhit - mu1) / sigma1
                            << " sigma." << std::endl;
                  std::cerr << "  Other branch of this cell:"
                            << "  subcell " << (i+2)%3
                            << " with nhit=" << subcell[(i+2)%3]->nhit
                            << ", expected " << mu2 << " +/- " << sigma2
                            << ", " << (subcell[(i+2)%3]->nhit - mu2) / sigma2
                            << " sigma." << std::endl;
                  warnings++;
               }
            }
         }
         for (int i=0; i < 3; ++i) {
            warnings += subcell[i]->check_subsets(optimized,
                                                  id + std::to_string(i));
         }
         return warnings;
      }
   };
};
