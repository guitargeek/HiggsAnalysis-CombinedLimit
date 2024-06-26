#ifndef HiggsAnalysis_CombinedLimit_FitterAlgoBase_h
#define HiggsAnalysis_CombinedLimit_FitterAlgoBase_h
/** \class FitterAlgoBase
 *
 * Do a ML fit of the data with background and signal+background hypothesis and print out diagnostics plots 
 *
 * \author Giovanni Petrucciani (UCSD)
 *
 *
 */
#include "HiggsAnalysis/CombinedLimit/interface/LimitAlgo.h"
#include "HiggsAnalysis/CombinedLimit/interface/Significance.h"
class RooFitResult;
class RooMinimizer;
class RooCmdArg;
class RooAbsReal;
class RooArgList;
class CascadeMinimizer;
#include <RooArgSet.h>

class FitterAlgoBase : public LimitAlgo {
public:
  FitterAlgoBase(const char *title="<FillMe> specific options") ;

  void applyOptionsBase(const boost::program_options::variables_map &vm) ;

  // configures the minimizer and then calls runSpecific
  virtual bool run(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint);

protected:
  //static std::string minimizerAlgo_, 
  static std::string minimizerAlgoForMinos_;
  //static float       minimizerTolerance_, 
  static float 	     minimizerToleranceForMinos_;
  static float 	     crossingTolerance_;
  //static int         minimizerStrategy_, 
  static int 	     minimizerStrategyForMinos_;

  static float preFitValue_;

  static bool robustFit_, do95_, forceRecreateNLL_;
  static float stepSize_;
  static int   maxFailedSteps_;

  enum ProfilingMode { ProfileAll, ProfileUnconstrained, ProfilePOI, NoProfiling };
  static ProfilingMode profileMode_;
  RooArgSet parametersToFreeze_;

  static bool  saveNLL_, keepFailures_, protectUnbinnedChannels_;
  static std::string autoBoundsPOIs_, autoMaxPOIs_;
  RooArgSet autoBoundsPOISet_, autoMaxPOISet_;
  static double nllValue_, nll0Value_;
  std::unique_ptr<RooAbsReal> nll;

  RooArgSet allParameters_;

  // method that is implemented in the subclass
  virtual bool runSpecific(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) = 0;

  // utilities
  /// Fit data with pdf, with parameters of interest in r, and specified constraint
  /// If ndim = 1, errors on each parameter are from a 1-dim chisquare, as for a single parameter fit
  /// If ndim > 1, errors on each parameter are from a n-dim chisquare, as for a joint estimation of N parameters 
  RooFitResult *doFit(RooAbsPdf &pdf, RooAbsData &data, RooRealVar &r,  const RooCmdArg &constrain, bool doHesse=true, int ndim=1,bool reuseNLL=false, bool saveFitResult=true) ;
  RooFitResult *doFit(RooAbsPdf &pdf, RooAbsData &data, const RooArgList &rs, const RooCmdArg &constrain, bool doHesse=true, int ndim=1,bool reuseNLL=false, bool saveFitResult=true) ;
  double findCrossing(CascadeMinimizer &minim, RooAbsReal &nll, RooRealVar &r, double level, double rStart, double rBound) ;
  double findCrossingNew(CascadeMinimizer &minim, RooAbsReal &nll, RooRealVar &r, double level, double rStart, double rBound) ;

  void optimizeBounds(const RooWorkspace *w, const RooStats::ModelConfig *mc) ;
  void restoreBounds(const RooWorkspace *w, const RooStats::ModelConfig *mc) ;
};


#endif
