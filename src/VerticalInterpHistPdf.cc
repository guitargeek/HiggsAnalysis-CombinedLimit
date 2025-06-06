#include "../interface/VerticalInterpHistPdf.h"

#include <cassert>
#include <memory>

#include "RooFit.h"
#include "Riostream.h"

#include "RooRealVar.h"
#include "RooMsgService.h"
#include "RooAbsData.h"

//#define TRACE_CALLS
#ifdef TRACE_CALLS
#include "../interface/ProfilingTools.h"
#define TRACEME()   PerfCounter::add( __PRETTY_FUNCTION__ );
#else
#define TRACEME() 
#endif


#define PATCH_FOR_HZZ_TEMPLATES
#ifdef PATCH_FOR_HZZ_TEMPLATES
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "../interface/utils.h"
namespace {
    std::unique_ptr<TH1> safeCreateHist2D(RooAbsPdf *pdf, const RooRealVar &x, const RooRealVar &y, bool conditional) {
        if (!pdf->getAttribute("safeCreateHist2D:ok") && typeid(*pdf) == typeid(RooHistPdf)) {
            RooHistPdf *hpdf = static_cast<RooHistPdf *>(pdf);
            RooDataHist &dataHist = hpdf->dataHist();
            bool ok = true;
            for (int i = 0; i < 2; ++i) {
                const RooRealVar *v = (i ? &y : &x);
                RooRealVar* lvarg = dynamic_cast<RooRealVar*>(dataHist.get()->find(v->GetName()));
                const RooAbsBinning* binning = lvarg->getBinningPtr(0);
                if (binning->numBins() != lvarg->numBins() || 
                    binning->binLow(0) != lvarg->getMin()  ||
                    binning->binHigh(binning->numBins()-1) != lvarg->getMax()) {
                    std::cout << "ERROR: inconsistent binning of RooDataHist " << dataHist.GetName() << ", var  " << lvarg->GetName() << std::endl;
                    std::cout << "  bins: " << binning->numBins() << " (binning) vs " << lvarg->numBins() << " (var)" << std::endl;
                    std::cout << "  min:  " << binning->binLow(0) << " (binning) vs " << lvarg->getMin() << " (var)" << std::endl;
                    std::cout << "  max:  " << binning->binHigh(binning->numBins()-1) << " (binning) vs " << lvarg->getMax() << " (var)" << std::endl;
                    ok = false;
                }
            }
            if (!ok) {
                std::cout << "BINNED DATASET: " << std::endl;
                utils::printRDH(&dataHist);
                const RooAbsBinning* xbinning = x.getBinningPtr(0);
                const RooAbsBinning* ybinning = y.getBinningPtr(0);
                if (x.numBins() != xbinning->numBins()) assert(0);
                if (y.numBins() != ybinning->numBins()) assert(0);
                double xdelta = x.getMin() - xbinning->binLow(0);
                double ydelta = y.getMin() - ybinning->binLow(0);
                double * xarray = xbinning->array(),  * yarray = ybinning->array();
                std::vector<double> xbins(xarray, xarray+(x.numBins()+1));
                std::vector<double> ybins(yarray, yarray+(y.numBins()+1));
                for (double &xe : xbins) xe += xdelta;
                for (double &ye : ybins) ye += ydelta;
                std::unique_ptr<TH1> hist(new TH2F("","",x.numBins(),&xbins[0],y.numBins(),&ybins[0]));
                TH2F *h2d = static_cast<TH2F*>(hist.get()); h2d->SetDirectory(0);
                TAxis *xaxis = h2d->GetXaxis(), * yaxis = h2d->GetYaxis();
                for (unsigned int id = 0, nd = dataHist.numEntries(); id < nd; ++id) {
                    const RooArgSet *point = dataHist.get(id);
                    double weight = dataHist.weight();
                    double xval = point->getRealValue(x.GetName());
                    double yval = point->getRealValue(y.GetName());
                    int bx = xaxis->FindBin(xval);
                    int by = yaxis->FindBin(yval);
                    if (bx == 0 || bx > x.numBins() || fabs(xaxis->GetBinCenter(bx)-xval) > 1e-4*std::max(1.0,std::abs(xval))) {
                        std::cout << "ERROR: dataset entry inconsistent with bin center along X" << std::endl;
                        point->Print("V");
                        assert(0); 
                    }
                    if (by == 0 || by > y.numBins() || fabs(yaxis->GetBinCenter(by)-yval) > 1e-4*std::max(1.0,std::abs(yval))) {
                        std::cout << "ERROR: dataset entry inconsistent with bin center along Y" << std::endl;
                        point->Print("V");
                        assert(0); 
                    }
                    h2d->Fill(xval,yval,weight); 
                }
                std::cout << "RECOVERED TEMPLATE FROM SLOW FILL" << std::endl;
                return hist;
            } else {
                pdf->setAttribute("safeCreateHist2D:ok");
            }
        }
        const RooCmdArg &cond = conditional ? RooFit::ConditionalObservables(RooArgSet(x)) : RooCmdArg::none();
        return std::unique_ptr<TH1>(pdf->createHistogram("", x, RooFit::YVar(y), cond));
    }
}
#endif

ClassImp(VerticalInterpHistPdf)


//_____________________________________________________________________________
VerticalInterpHistPdf::VerticalInterpHistPdf(const char *name, const char *title, const RooRealVar &x, const RooArgList& inFuncList, const RooArgList& inCoefList, Double_t smoothRegion, Int_t smoothAlgo) :
  RooAbsPdf(name,title),
  _x("x","Independent variable",this,const_cast<RooRealVar&>(x)),
  _funcList("funcList","List of pdfs",this),
  _coefList("coefList","List of coefficients",this), // we should get shapeDirty when coefficients change
  _smoothRegion(smoothRegion),
  _smoothAlgo(smoothAlgo)
{ 

  if (inFuncList.getSize()!=2*inCoefList.getSize()+1) {
    coutE(InputArguments) << "VerticalInterpHistPdf::VerticalInterpHistPdf(" << GetName() 
			  << ") number of pdfs and coefficients inconsistent, must have Nfunc=1+2*Ncoef" << std::endl ;
    assert(0);
  }

  for (RooAbsArg *func : inFuncList) {
    RooAbsPdf *pdf = dynamic_cast<RooAbsPdf*>(func);
    if (!pdf) {
      coutE(InputArguments) << "ERROR: VerticalInterpHistPdf::VerticalInterpHistPdf(" << GetName() << ") function  " << func->GetName() << " is not of type RooAbsPdf" << std::endl;
      assert(0);
    }
    std::unique_ptr<RooArgSet> params{pdf->getParameters(RooArgSet(x))};
    if (params->getSize() > 0) {
      coutE(InputArguments) << "ERROR: VerticalInterpHistPdf::VerticalInterpHistPdf(" << GetName() << ") pdf  " << func->GetName() << " has some parameters." << std::endl;
      assert(0);
    }
    _funcList.add(*func) ;
  }

  for (RooAbsArg *coef : inCoefList) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      coutE(InputArguments) << "ERROR: VerticalInterpHistPdf::VerticalInterpHistPdf(" << GetName() << ") coefficient " << coef->GetName() << " is not of type RooAbsReal" << std::endl;
      assert(0);
    }
    _coefList.add(*coef) ;    
  }
}




//_____________________________________________________________________________
VerticalInterpHistPdf::VerticalInterpHistPdf(const VerticalInterpHistPdf& other, const char* name) :
  RooAbsPdf(other,name),
  _x("x",this,other._x),
  _funcList("funcList",this,other._funcList),
  _coefList("coefList",this,other._coefList),
  _smoothRegion(other._smoothRegion),
  _smoothAlgo(other._smoothAlgo)
{
  // Copy constructor
}


//_____________________________________________________________________________
Double_t VerticalInterpHistPdf::evaluate() const 
{
  if (!_cacheTotal) setupCaches();
#if 0
  printf("Evaluate called for x %f, cache status %d: ", _x.arg().getVal(), _sentry.good());
  int ndim = _coefList.getSize();
  for (int i = 0; i < ndim; ++i) { 
    RooAbsReal *rar = dynamic_cast<RooAbsReal *>(_coefList.at(i));
    printf("%s = %+6.4f  ", rar->GetName(), rar->getVal());
  }
  std::cout << std::endl;
#endif

  if (!_sentry.good()) syncTotal();
  int nbin = _cacheTotal->GetNbinsX(), ibin = _cacheTotal->FindBin(_x);
  if (ibin < 1) ibin = 1;
  else if (ibin > nbin) ibin = nbin;
  return _cacheTotal->GetBinContent(ibin);
}


void VerticalInterpHistPdf::syncComponent(int i) const {
    RooAbsPdf *pdfi = dynamic_cast<RooAbsPdf *>(_funcList.at(i));
    _cacheSingle[i] = std::unique_ptr<TH1>{pdfi->createHistogram("",dynamic_cast<const RooRealVar &>(_x.arg()))};
    _cacheSingle[i]->SetDirectory(0);
    if (_cacheSingle[i]->Integral("width")) { _cacheSingle[i]->Scale(1.0/_cacheSingle[i]->Integral("width")); }
    if (i > 0) {
        for (int b = 1, nb = _cacheSingle[i]->GetNbinsX(); b <= nb; ++b) {
            double y  = _cacheSingle[i]->GetBinContent(b);
            double y0 = _cacheSingle[0]->GetBinContent(b);
            if (_smoothAlgo < 0) {
                if (y > 0 && y0 > 0) {
                    double logk = log(y/y0);
                    // odd numbers correspond to up variations, even numbers to down variations,
                    // and down variations need -log(kappa) instead of log(kappa)
                    _cacheSingle[i]->SetBinContent(b, logk);
                } else {
                    _cacheSingle[i]->SetBinContent(b, 0);
                }
            } else {
                _cacheSingle[i]->SetBinContent(b, y - y0);
            }
        }
    }
    _cacheSingleGood[i] = true;
}

void VerticalInterpHistPdf::syncTotal() const {
    int ndim = _coefList.getSize();
    for (int i = 0; i < 2*ndim+1; ++i) {
        if (!_cacheSingleGood[i]) syncComponent(i);
    }
    for (int b = 1, nb = _cacheTotal->GetNbinsX(); b <= nb; ++b) {
        double val = _cacheSingle[0]->GetBinContent(b);
        int i = 0;
        for (RooAbsArg *coef : _coefList) {
            double dhi = _cacheSingle[2*i+1]->GetBinContent(b);
            double dlo = _cacheSingle[2*i+2]->GetBinContent(b);
            double x = dynamic_cast<RooAbsReal *>(coef)->getVal();
            double alpha = x * 0.5 * ((dhi-dlo) + (dhi+dlo)*smoothStepFunc(x));
            // alpha(0) = 0
            // alpha(+1) = dhi 
            // alpha(-1) = dlo
            // alpha(x >= +1) = |x|*dhi
            // alpha(x <= -1) = |x|*dlo
            // alpha is continuous and has continuous first and second derivative, as smoothStepFunc has them
            if (_smoothAlgo < 0) {
                val *= exp(alpha);
            } else {
                val += alpha; 
            }
            ++i;
        }    
        if (val <= 0) val = 1e-9;
        _cacheTotal->SetBinContent(b, val);
    }
    double norm = _cacheTotal->Integral("width");
    if (norm > 0) _cacheTotal->Scale(1.0/norm);
    _sentry.reset();
}

void VerticalInterpHistPdf::setupCaches() const {
    int ndim = _coefList.getSize();
    _cacheTotal = std::unique_ptr<TH1>{(dynamic_cast<const RooRealVar &>(_x.arg())).createHistogram("total")};
    _cacheTotal->SetDirectory(0);
    _cacheSingle.resize(2*ndim+1);
    _cacheSingleGood.resize(2*ndim+1);
    for (int i = 0; i < 2*ndim+1; ++i) { 
        _cacheSingle[i] = nullptr;
        _cacheSingleGood[i] = 0; 
        syncComponent(i);  
    } 
    if (_sentry.empty()) _sentry.addVars(_coefList); 
    syncTotal();
}

//=============================================================================================
ClassImp(FastVerticalInterpHistPdfBase)
ClassImp(FastVerticalInterpHistPdf)
ClassImp(FastVerticalInterpHistPdf2D)
ClassImp(FastVerticalInterpHistPdf3D)


//_____________________________________________________________________________
FastVerticalInterpHistPdfBase::FastVerticalInterpHistPdfBase() 
{
  // Default constructor
}


//_____________________________________________________________________________
FastVerticalInterpHistPdfBase::FastVerticalInterpHistPdfBase(const char *name, const char *title, const RooArgSet &obs, const RooArgList& inFuncList, const RooArgList& inCoefList, Double_t smoothRegion, Int_t smoothAlgo) :
  RooAbsPdf(name,title),
  _funcList("funcList","List of pdfs",this),
  _coefList("coefList","List of coefficients",this), // we should get shapeDirty when coefficients change
  _smoothRegion(smoothRegion),
  _smoothAlgo(smoothAlgo),
  _init(false),
  _morphs(), _morphParams()
{ 
  TRACEME()

  if (inFuncList.getSize()!=2*inCoefList.getSize()+1) {
    coutE(InputArguments) << "VerticalInterpHistPdf::VerticalInterpHistPdf(" << GetName() 
			  << ") number of pdfs and coefficients inconsistent, must have Nfunc=1+2*Ncoef" << std::endl ;
    assert(0);
  }

  for (RooAbsArg *func : inFuncList) {
    RooAbsPdf *pdf = dynamic_cast<RooAbsPdf*>(func);
    if (!pdf) {
      coutE(InputArguments) << "ERROR: VerticalInterpHistPdf::VerticalInterpHistPdf(" << GetName() << ") function  " << func->GetName() << " is not of type RooAbsPdf" << std::endl;
      assert(0);
    }
    RooArgSet *params = pdf->getParameters(obs);
    if (params->getSize() > 0) {
      coutE(InputArguments) << "ERROR: VerticalInterpHistPdf::VerticalInterpHistPdf(" << GetName() << ") pdf  " << func->GetName() << " has some parameters." << std::endl;
      assert(0);
    }
    delete params;
    _funcList.add(*func) ;
  }

  for (RooAbsArg *coef : inCoefList) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      coutE(InputArguments) << "ERROR: VerticalInterpHistPdf::VerticalInterpHistPdf(" << GetName() << ") coefficient " << coef->GetName() << " is not of type RooAbsReal" << std::endl;
      assert(0);
    }
    _coefList.add(*coef) ;    
  }
}

//_____________________________________________________________________________
FastVerticalInterpHistPdfBase::FastVerticalInterpHistPdfBase(const FastVerticalInterpHistPdfBase& other, const char* name) :
  RooAbsPdf(other,name),
  _funcList("funcList",this,other._funcList),
  _coefList("coefList",this,other._coefList),
  _smoothRegion(other._smoothRegion),
  _smoothAlgo(other._smoothAlgo),
  _init(false),
  _morphs(other._morphs), _morphParams(other._morphParams)
{
  // Copy constructor
  _sentry.addVars(_coefList);
  _sentry.setValueDirty(); 
}


//_____________________________________________________________________________
FastVerticalInterpHistPdfBase::~FastVerticalInterpHistPdfBase() = default;


//_____________________________________________________________________________
Double_t FastVerticalInterpHistPdf::evaluate() const 
{
  TRACEME()
  if (_cache.size() == 0) setupCaches();
  if (!_sentry.good() || !_init) syncTotal();
  //return std::max<double>(1e-9, _cache.GetAt(_x));
  return _cache.GetAt(_x);
}

//_____________________________________________________________________________
Double_t FastVerticalInterpHistPdf2D::evaluate() const 
{
  TRACEME()
  if (_cache.size() == 0) setupCaches();

  if (!_sentry.good() || !_init) syncTotal();
  //return std::max<double>(1e-9, _cache.GetAt(_x, _y));
  return _cache.GetAt(_x, _y);
}
Double_t FastVerticalInterpHistPdf3D::evaluate() const 
{
  TRACEME()
  if (_cache.size() == 0) setupCaches();

  if (!_sentry.good() || !_init) syncTotal();
  //return std::max<double>(1e-9, _cache.GetAt(_x, _y));
  return _cache.GetAt(_x, _y, _z);
}



void FastVerticalInterpHistPdf::syncNominal() const {
    TRACEME()
    RooAbsPdf *pdf = dynamic_cast<RooAbsPdf *>(_funcList.at(0));
    const RooRealVar &x = dynamic_cast<const RooRealVar &>(_x.arg());
    std::unique_ptr<TH1> hist(pdf->createHistogram("",x));
    hist->SetDirectory(0); 
    _cacheNominal = FastHisto(*hist);
    _cacheNominal.Normalize();
    if (_smoothAlgo < 0) {
        _cacheNominalLog = _cacheNominal;
        _cacheNominalLog.Log();
    }
}

void FastVerticalInterpHistPdf2D::syncNominal() const {
    TRACEME()
    RooAbsPdf *pdf = dynamic_cast<RooAbsPdf *>(_funcList.at(0));
    const RooRealVar &x = dynamic_cast<const RooRealVar &>(_x.arg());
    const RooRealVar &y = dynamic_cast<const RooRealVar &>(_y.arg());
#ifdef PATCH_FOR_HZZ_TEMPLATES
    std::unique_ptr<TH1> hist(::safeCreateHist2D(pdf,x,y,_conditional));
#else
    const RooCmdArg &cond = _conditional ? RooFit::ConditionalObservables(RooArgSet(x)) : RooCmdArg::none();
    std::unique_ptr<TH1> hist(pdf->createHistogram("", x, RooFit::YVar(y), cond));
#endif
    hist->SetDirectory(0); 
    _cacheNominal = FastHisto2D(dynamic_cast<TH2F&>(*hist), _conditional);
    if (_conditional) _cacheNominal.NormalizeXSlices(); 
    else              _cacheNominal.Normalize(); 

    if (_smoothAlgo < 0) {
        _cacheNominalLog = _cacheNominal;
        _cacheNominalLog.Log();
    }
}


void FastVerticalInterpHistPdf3D::syncNominal() const {
    TRACEME()
    RooAbsPdf *pdf = dynamic_cast<RooAbsPdf *>(_funcList.at(0));
    const RooRealVar &x = dynamic_cast<const RooRealVar &>(_x.arg());
    const RooRealVar &y = dynamic_cast<const RooRealVar &>(_y.arg());
    const RooRealVar &z = dynamic_cast<const RooRealVar &>(_z.arg());
    const RooCmdArg &cond = _conditional ? RooFit::ConditionalObservables(RooArgSet(x)) : RooCmdArg::none();
    std::unique_ptr<TH1> hist(pdf->createHistogram("", x, RooFit::YVar(y), RooFit::ZVar(z),cond));
    hist->SetDirectory(0); 
    _cacheNominal = FastHisto3D(dynamic_cast<TH3F&>(*hist), _conditional);
    if (_conditional) _cacheNominal.NormalizeXSlices(); 
    else              _cacheNominal.Normalize(); 

    if (_smoothAlgo < 0) {
        _cacheNominalLog = _cacheNominal;
        _cacheNominalLog.Log();
    }
}

void FastVerticalInterpHistPdfBase::syncMorph(Morph &out, const FastTemplate &nominal, FastTemplate &lo, FastTemplate &hi) const {
    if (_smoothAlgo < 0)  {
        hi.LogRatio(nominal); lo.LogRatio(nominal);
        //printf("Log-ratios for dimension %d: \n", dim);  hi.Dump(); lo.Dump();
    } else {
        hi.Subtract(nominal); lo.Subtract(nominal);
        //printf("Differences for dimension %d: \n", dim);  hi.Dump(); lo.Dump();
    }
    FastTemplate::SumDiff(hi, lo, out.sum, out.diff);
    //printf("Sum and diff for dimension %d: \n", dim);  out.sum.Dump(); out.diff.Dump();
}

void FastVerticalInterpHistPdf::syncComponents(int dim) const {
    TRACEME()
    RooAbsPdf *pdfHi = dynamic_cast<RooAbsPdf *>(_funcList.at(2*dim+1));
    RooAbsPdf *pdfLo = dynamic_cast<RooAbsPdf *>(_funcList.at(2*dim+2));
    const RooRealVar &x = dynamic_cast<const RooRealVar &>(_x.arg());
    std::unique_ptr<TH1> histHi(pdfHi->createHistogram("",x)); histHi->SetDirectory(0); 
    std::unique_ptr<TH1> histLo(pdfLo->createHistogram("",x)); histLo->SetDirectory(0);

    FastHisto hi(*histHi), lo(*histLo); 
    //printf("Un-normalized templates for dimension %d: \n", dim);  hi.Dump(); lo.Dump();
    hi.Normalize(); lo.Normalize();
    //printf("Normalized templates for dimension %d: \n", dim);  hi.Dump(); lo.Dump();
    syncMorph(_morphs[dim], _cacheNominal, lo, hi);
}
void FastVerticalInterpHistPdf2D::syncComponents(int dim) const {
    RooAbsPdf *pdfHi = dynamic_cast<RooAbsPdf *>(_funcList.at(2*dim+1));
    RooAbsPdf *pdfLo = dynamic_cast<RooAbsPdf *>(_funcList.at(2*dim+2));
    const RooRealVar &x = dynamic_cast<const RooRealVar &>(_x.arg());
    const RooRealVar &y = dynamic_cast<const RooRealVar &>(_y.arg());
#ifdef PATCH_FOR_HZZ_TEMPLATES
    std::unique_ptr<TH1> histHi(::safeCreateHist2D(pdfHi,x,y,_conditional)); histHi->SetDirectory(0);
    std::unique_ptr<TH1> histLo(::safeCreateHist2D(pdfLo,x,y,_conditional)); histLo->SetDirectory(0);
#else
    const RooCmdArg &cond = _conditional ? RooFit::ConditionalObservables(RooArgSet(x)) : RooCmdArg::none();
    std::unique_ptr<TH1> histHi(pdfHi->createHistogram("", x, RooFit::YVar(y), cond)); histHi->SetDirectory(0); 
    std::unique_ptr<TH1> histLo(pdfLo->createHistogram("", x, RooFit::YVar(y), cond)); histLo->SetDirectory(0);
#endif

    FastHisto2D hi(dynamic_cast<TH2&>(*histHi), _conditional), lo(dynamic_cast<TH2&>(*histLo), _conditional); 
    //printf("Un-normalized templates for dimension %d: \n", dim);  hi.Dump(); lo.Dump();
    if (_conditional) {
        hi.NormalizeXSlices(); lo.NormalizeXSlices();
    } else {
        hi.Normalize(); lo.Normalize();
    }
    //printf("Normalized templates for dimension %d: \n", dim);  hi.Dump(); lo.Dump();
    syncMorph(_morphs[dim], _cacheNominal, lo, hi);
}
void FastVerticalInterpHistPdf3D::syncComponents(int dim) const {
    RooAbsPdf *pdfHi = dynamic_cast<RooAbsPdf *>(_funcList.at(2*dim+1));
    RooAbsPdf *pdfLo = dynamic_cast<RooAbsPdf *>(_funcList.at(2*dim+2));
    const RooRealVar &x = dynamic_cast<const RooRealVar &>(_x.arg());
    const RooRealVar &y = dynamic_cast<const RooRealVar &>(_y.arg());
    const RooRealVar &z = dynamic_cast<const RooRealVar &>(_z.arg());
    const RooCmdArg &cond = _conditional ? RooFit::ConditionalObservables(RooArgSet(x)) : RooCmdArg::none();
    std::unique_ptr<TH1> histHi(pdfHi->createHistogram("", x, RooFit::YVar(y),RooFit::ZVar(z), cond)); histHi->SetDirectory(0); 
    std::unique_ptr<TH1> histLo(pdfLo->createHistogram("", x, RooFit::YVar(y),RooFit::ZVar(z), cond)); histLo->SetDirectory(0);

    FastHisto3D hi(dynamic_cast<TH3&>(*histHi), _conditional), lo(dynamic_cast<TH3&>(*histLo), _conditional); 
    //printf("Un-normalized templates for dimension %d: \n", dim);  hi.Dump(); lo.Dump();
    if (_conditional) {
        hi.NormalizeXSlices(); lo.NormalizeXSlices();
    } else {
        hi.Normalize(); lo.Normalize();
    }
    //printf("Normalized templates for dimension %d: \n", dim);  hi.Dump(); lo.Dump();
    syncMorph(_morphs[dim], _cacheNominal, lo, hi);
}


void FastVerticalInterpHistPdfBase::syncTotal(FastTemplate &cache, const FastTemplate &cacheNominal, const FastTemplate &cacheNominalLog) const {
    TRACEME()
    /* === how the algorithm works, in theory ===
     * let  dhi = h_hi - h_nominal
     *      dlo = h_lo - h_nominal
     * and x be the morphing parameter
     * we define alpha = x * 0.5 * ((dhi-dlo) + (dhi+dlo)*smoothStepFunc(x));
     * which satisfies:
     *     alpha(0) = 0
     *     alpha(+1) = dhi 
     *     alpha(-1) = dlo
     *     alpha(x >= +1) = |x|*dhi
     *     alpha(x <= -1) = |x|*dlo
     *     alpha is continuous and has continuous first and second derivative, as smoothStepFunc has them
     * === and in practice ===
     * we already have computed the histogram for diff=(dhi-dlo) and sum=(dhi+dlo)
     * so we just do template += (0.5 * x) * (diff + smoothStepFunc(x) * sum)
     * ========================================== */

    // start from nominal
    cache.CopyValues(_smoothAlgo < 0 ? cacheNominalLog : cacheNominal);
    //printf("Cache initialized to nominal template: \n");  cacheNominal.Dump();

    // apply all morphs one by one
    for (int i = 0, ndim = _coefList.getSize(); i < ndim; ++i) {
        double x = _morphParams[i]->getVal();
        double a = 0.5*x, b = smoothStepFunc(x);
        cache.Meld(_morphs[i].diff, _morphs[i].sum, a, b);    
        //printf("Merged transformation for dimension %d, x = %+5.3f, step = %.3f: \n", i, x, b);  cache.Dump();
    }

    // if necessary go back to linear scale
    if (_smoothAlgo < 0) {
        cache.Exp();
        //printf("Done exponential tranformation\n");  cache.Dump();
    } else {
        cache.CropUnderflows();
    }
    
    // mark as done
    _sentry.reset();
    _init = true;
}

void FastVerticalInterpHistPdf::syncTotal() const {
    FastVerticalInterpHistPdfBase::syncTotal(_cache, _cacheNominal, _cacheNominalLog);

    // normalize the result
    _cache.Normalize(); 
    //printf("Normalized result\n");  _cache.Dump();
}

void FastVerticalInterpHistPdf2D::syncTotal() const {
    FastVerticalInterpHistPdfBase::syncTotal(_cache, _cacheNominal, _cacheNominalLog);
    // normalize the result
    if (_conditional) _cache.NormalizeXSlices(); 
    else              _cache.Normalize(); 
    //printf("Normalized result\n");  _cache.Dump();
}

void FastVerticalInterpHistPdf3D::syncTotal() const {
    FastVerticalInterpHistPdfBase::syncTotal(_cache, _cacheNominal, _cacheNominalLog);
    // normalize the result
    if (_conditional) _cache.NormalizeXSlices(); 
    else              _cache.Normalize(); 
    //printf("Normalized result\n");  _cache.Dump();
}


void FastVerticalInterpHistPdf::setupCaches() const {
    TRACEME()
    int ndim = _coefList.getSize();

    _morphs.resize(ndim);
    _morphParams.resize(ndim);
    syncNominal();
    //printf("Nominal template has been set up: \n");  _cacheNominal.Dump();
    int i = 0;
    for (RooAbsArg *a : _coefList) {
        _morphParams[i] = dynamic_cast<RooAbsReal *>(a);
        _morphs[i].sum.Resize(_cacheNominal.size());
        _morphs[i].diff.Resize(_cacheNominal.size());
        syncComponents(i);
        ++i;
    } 
    _cache = FastHisto(_cacheNominal);

    if (_sentry.empty()) _sentry.addVars(_coefList); 
    syncTotal();
}

void FastVerticalInterpHistPdf2D::setupCaches() const {
    TRACEME()
    int ndim = _coefList.getSize();

    _morphs.resize(ndim);
    _morphParams.resize(ndim);
    syncNominal();
    //printf("Nominal template has been set up: \n");  _cacheNominal.Dump();
    int i = 0;
    for (RooAbsArg *a : _coefList) {
        _morphParams[i] = dynamic_cast<RooAbsReal *>(a);
        _morphs[i].sum.Resize(_cacheNominal.size());
        _morphs[i].diff.Resize(_cacheNominal.size());
        syncComponents(i);
        ++i;
    } 
    _cache = FastHisto2D(_cacheNominal);

    if (_sentry.empty()) _sentry.addVars(_coefList); 
    syncTotal();
}
void FastVerticalInterpHistPdf3D::setupCaches() const {
    TRACEME()
    int ndim = _coefList.getSize();

    _morphs.resize(ndim);
    _morphParams.resize(ndim);
    syncNominal();
    //printf("Nominal template has been set up: \n");  _cacheNominal.Dump();
    int i = 0;
    for (RooAbsArg *a : _coefList) {
        _morphParams[i] = dynamic_cast<RooAbsReal *>(a);
        _morphs[i].sum.Resize(_cacheNominal.size());
        _morphs[i].diff.Resize(_cacheNominal.size());
        syncComponents(i);
        ++i;
    } 
    _cache = FastHisto3D(_cacheNominal);

    if (_sentry.empty()) _sentry.addVars(_coefList); 
    syncTotal();
}


FastVerticalInterpHistPdfV::FastVerticalInterpHistPdfV(const FastVerticalInterpHistPdf &hpdf, const RooAbsData &data, bool includeZeroWeights) :
    hpdf_(hpdf)
{
    // check init
    if (hpdf._cache.size() == 0) hpdf.setupCaches();
    if (!hpdf._sentry.good() || !hpdf._init) hpdf.syncTotal();
    // find bins
    std::vector<int> bins;
    RooArgSet obs(hpdf._x.arg());
    const RooRealVar &x = static_cast<const RooRealVar &>(*obs.first());
    bool aligned = true;
    for (int i = 0, n = data.numEntries(); i < n; ++i) {
        obs = *data.get(i);
        if (data.weight() == 0 && !includeZeroWeights) continue;
        int idx = hpdf._cache.FindBin(x.getVal());
        if (!bins.empty() && idx != bins.back() + 1) aligned = false;
        bins.push_back(idx);
    }
    if (aligned) {
        begin_ = bins.front();
        end_   = bins.back()+1;
        //std::cout << "Created FastVerticalInterpHistPdfV from " << hpdf.GetName() << ", aligned, " << (end_-begin_) << " bins." << std::endl;
    } else {
        nbins_ = bins.size();
        bins_.swap(bins);
        blocks_.clear();
        int start = bins_[0], istart = 0;
        for (int i = 1, n = bins_.size(); i < n; ++i) {
            if (bins_[i] != bins_[i-1]+1) { 
                blocks_.push_back(Block(istart,start,bins_[i-1]+1));
                start = bins_[i];
                istart = i;
            }
        }
        blocks_.push_back(Block(istart,start,bins_.back()+1));
        if (blocks_.size() < 4*bins_.size()) {
            //std::cout << "Created FastVerticalInterpHistPdfV from " << hpdf.GetName() << ", block-aligned, " << bins_.size() << " bins, " <<  blocks_.size() << " blocks." << std::endl;
            bins_.clear();
        } else {
            //std::cout << "Created FastVerticalInterpHistPdfV from " << hpdf.GetName() << ", non-aligned, " << bins_.size() << " bins, " <<  blocks_.size() << " blocks." << std::endl;
            blocks_.clear();
        }
    }
}

void FastVerticalInterpHistPdfV::fill(std::vector<Double_t> &out) const 
{
    if (!hpdf_._sentry.good()) hpdf_.syncTotal();
    if (begin_ != end_) {
        out.resize(end_-begin_);
        std::copy(& hpdf_._cache.GetBinContent(begin_), (&hpdf_._cache.GetBinContent(end_-1))+1, out.begin());
    } else if (!blocks_.empty()) {
        out.resize(nbins_);
        for (auto b : blocks_) std::copy(& hpdf_._cache.GetBinContent(b.begin), (&hpdf_._cache.GetBinContent(b.end-1))+1, out.begin()+b.index);
    } else {
        out.resize(bins_.size());
        for (int i = 0, n = bins_.size(); i < n; ++i) {
            // if ((int)hpdf_._cache.GetNbinsX()>bins_[i]) out[i] = hpdf_._cache.GetBinContent(bins_[i]);
            // else out[i] = 0;
            out[i] = hpdf_._cache.GetBinContent(bins_[i]);
        }
    }
}

//=============================================================================================
ClassImp(FastVerticalInterpHistPdf2Base)
ClassImp(FastVerticalInterpHistPdf2)
//ClassImp(FastVerticalInterpHistPdf2D)


//_____________________________________________________________________________
FastVerticalInterpHistPdf2Base::FastVerticalInterpHistPdf2Base() :
    _initBase(false)
{
  // Default constructor
}


//_____________________________________________________________________________
FastVerticalInterpHistPdf2Base::FastVerticalInterpHistPdf2Base(const char *name, const char *title, const RooArgSet &obs, const TList& inFuncList, const RooArgList& inCoefList, Double_t smoothRegion, Int_t smoothAlgo) :
  RooAbsPdf(name,title),
  _coefList("coefList","List of coefficients",this), // we should get shapeDirty when coefficients change
  _smoothRegion(smoothRegion),
  _smoothAlgo(smoothAlgo),
  _initBase(false),
  _morphs(), _morphParams()
{ 
  if (inFuncList.GetSize()!=2*inCoefList.getSize()+1) {
    coutE(InputArguments) << "VerticalInterpHistPdf::VerticalInterpHistPdf(" << GetName() 
			  << ") number of pdfs and coefficients inconsistent, must have Nfunc=1+2*Ncoef" 
                          << "while Nfunc= " << inFuncList.GetSize() << " and Ncoef= " << inCoefList.getSize() <<std::endl ;
    assert(0);
  }

  TIter funcIter(&inFuncList) ;
  TObject* func;
  while((func = (RooAbsArg*)funcIter.Next())) {
    TH1 *hist = dynamic_cast<TH1*>(func);
    RooAbsPdf *pdf = dynamic_cast<RooAbsPdf*>(func);
    if (!pdf && !hist) {
      coutE(InputArguments) << "ERROR: VerticalInterpHistPdf::VerticalInterpHistPdf(" << GetName() << ") function  " << func->GetName() << " is not of type TH1 or RooAbsPdf" << std::endl;
      assert(0);
    }
    if (pdf) {
        std::unique_ptr<RooArgSet> params{pdf->getParameters(obs)};
        if (params->getSize() > 0) {
          coutE(InputArguments) << "ERROR: VerticalInterpHistPdf::VerticalInterpHistPdf(" << GetName() << ") pdf  " << func->GetName() << " (" << func->ClassName()<<") has some parameters." << std::endl;
          obs.Print("");
          params->Print("");
          assert(0);
        }
    }
  }

  _coefList.add(inCoefList);
}

//_____________________________________________________________________________
FastVerticalInterpHistPdf2Base::FastVerticalInterpHistPdf2Base(const FastVerticalInterpHistPdf2Base& other, const char* name) :
  RooAbsPdf(other,name),
  _coefList("coefList", this, other._coefList),
  _smoothRegion(other._smoothRegion),
  _smoothAlgo(other._smoothAlgo),
  _initBase(false),
  _morphs(other._morphs), _morphParams(other._morphParams)
{
}

//_____________________________________________________________________________
FastVerticalInterpHistPdf2Base::FastVerticalInterpHistPdf2Base(const FastVerticalInterpHistPdfBase& other, const char* name) :
  RooAbsPdf(other,name),
  _coefList("coefList", this, other._coefList),
  _smoothRegion(other._smoothRegion),
  _smoothAlgo(other._smoothAlgo),
  _initBase(false),
  _morphs(), _morphParams()
{
  // Convert constructor
}


Bool_t FastVerticalInterpHistPdf2Base::importWorkspaceHook(RooWorkspace& ws) {
  _initBase = false;
  _morphParams.clear();
  _sentry.reset();
  return kFALSE;
}


//_____________________________________________________________________________
FastVerticalInterpHistPdf2Base::~FastVerticalInterpHistPdf2Base()
{
  // Destructor
}

void
FastVerticalInterpHistPdf2Base::initBase() const 
{
    if (_initBase) return;

    for (RooAbsArg *coef : _coefList) {
        const RooAbsReal *rrv = dynamic_cast<RooAbsReal*>(coef);
        if (!rrv) {
            coutE(InputArguments) << "ERROR: VerticalInterpHistPdf::VerticalInterpHistPdf(" << GetName() << ") coefficient " << coef->GetName() << " is not of type RooAbsReal" << std::endl;
            assert(0);
        }
        _morphParams.push_back(rrv);
    }

    _sentry.addVars(_coefList);
    _sentry.setValueDirty(); 
    _initBase = true;
}

FastVerticalInterpHistPdf2::FastVerticalInterpHistPdf2(const char *name, const char *title, const RooRealVar &x, const TList & funcList, const RooArgList& coefList, Double_t smoothRegion, Int_t smoothAlgo) :
    FastVerticalInterpHistPdf2Base(name,title,x,funcList,coefList,smoothRegion,smoothAlgo),
    _x("x","Independent variable",this,const_cast<RooRealVar&>(x)),
    _cache(), _cacheNominal(), _cacheNominalLog()
{
    initBase();
    initNominal(funcList.At(0));
    _morphs.resize(coefList.getSize());
    for (int i = 0, n = coefList.getSize(); i < n; ++i) {
        initComponent(i, funcList.At(2*i+1), funcList.At(2*i+2));
    }
}

namespace {

RooArgSet createRooArgSet(RooAbsArg const& arg1, RooAbsArg const& arg2) {
    RooArgSet out;
    out.add(arg1);
    out.add(arg2);
    return out;
}

} // namespace

FastVerticalInterpHistPdf2D2::FastVerticalInterpHistPdf2D2(const char *name, const char *title, const RooRealVar &x, const RooRealVar &y, bool conditional, const TList & funcList, const RooArgList& coefList, Double_t smoothRegion, Int_t smoothAlgo) :
    FastVerticalInterpHistPdf2Base(name,title,createRooArgSet(x, y),funcList,coefList,smoothRegion,smoothAlgo),
    _x("x","Independent variable",this,const_cast<RooRealVar&>(x)),
    _y("y","Independent variable",this,const_cast<RooRealVar&>(y)),
    _conditional(conditional),
    _cache(), _cacheNominal(), _cacheNominalLog()
{
    initBase();
    initNominal(funcList.At(0));
    _morphs.resize(coefList.getSize());
    for (int i = 0, n = coefList.getSize(); i < n; ++i) {
        initComponent(i, funcList.At(2*i+1), funcList.At(2*i+2));
    }
}



FastVerticalInterpHistPdf2::FastVerticalInterpHistPdf2(const FastVerticalInterpHistPdf& other, const char* name) :
    FastVerticalInterpHistPdf2Base(other,name),
    _x("x",this,other._x),
    _cache(), _cacheNominal(), _cacheNominalLog()
{
    initBase();
    other.getVal(_x.arg());
    _morphs = other._morphs;
    _cache = other._cache;
    _cacheNominal = other._cacheNominal;
    _cacheNominalLog = other._cacheNominalLog;
}


FastVerticalInterpHistPdf2D2::FastVerticalInterpHistPdf2D2(const FastVerticalInterpHistPdf2D& other, const char* name) :
    FastVerticalInterpHistPdf2Base(other,name),
    _x("x",this,other._x),
    _y("y",this,other._y),
    _conditional(other._conditional),
    _cache(), _cacheNominal(), _cacheNominalLog()
{
    initBase();
    RooArgSet normSet;
    normSet.add(_x.arg());
    normSet.add(_y.arg());
    other.getVal(normSet);
    _morphs = other._morphs;
    _cache = other._cache;
    _cacheNominal = other._cacheNominal;
    _cacheNominalLog = other._cacheNominalLog;
}



//_____________________________________________________________________________
Double_t FastVerticalInterpHistPdf2::evaluate() const 
{
  if (!_initBase) initBase();
  if (_cache.size() == 0) _cache = _cacheNominal; // _cache is not persisted
  if (!_sentry.good()) syncTotal();
  return _cache.GetAt(_x);
}
Double_t FastVerticalInterpHistPdf2D2::evaluate() const 
{
  if (!_initBase) initBase();
  if (_cache.size() == 0) _cache = _cacheNominal; // _cache is not persisted
  if (!_sentry.good()) syncTotal();
  return _cache.GetAt(_x,_y);
}


void FastVerticalInterpHistPdf2::setActiveBins(unsigned int bins) {
  assert(bins <= _cacheNominal.fullsize());
  if (_cache.size() == 0) _cache = _cacheNominal; // _cache is not persisted
  _cache.CropUnderflows(1e-9,false);        // set all bins, also the non-active ones
  _cacheNominal.CropUnderflows(1e-9,false); // set all bins, also the non-active ones
  _cache.SetActiveSize(bins);
  _cacheNominal.SetActiveSize(bins);
  _cacheNominalLog.SetActiveSize(bins);
  for (Morph & m : _morphs) {
    m.sum.SetActiveSize(bins);
    m.diff.SetActiveSize(bins);
  }
  //printf("Setting the number of active bins to be %d/%d for %s\n", bins, _cacheNominal.fullsize(), GetName());
}

void FastVerticalInterpHistPdf2::initNominal(TObject *templ) {
    TH1 *hist = dynamic_cast<TH1*>(templ);
    if (hist) {
        _cacheNominal = FastHisto(*hist);
    } else {
        RooAbsPdf *pdf = dynamic_cast<RooAbsPdf *>(templ);
        const RooRealVar &x = dynamic_cast<const RooRealVar &>(_x.arg());
        std::unique_ptr<TH1> hist{pdf->createHistogram("",x)};
        hist->SetDirectory(0); 
        _cacheNominal = FastHisto(*hist);
    }
    _cacheNominal.Normalize();
    //printf("Normalized nominal templated: \n");  _cacheNominal.Dump(); 
    if (_smoothAlgo < 0) {
        _cacheNominalLog = _cacheNominal;
        _cacheNominalLog.Log();
    }
    _cache = _cacheNominal;
}

void FastVerticalInterpHistPdf2D2::initNominal(TObject *templ) {
    if (TH2 *templHist = dynamic_cast<TH2*>(templ)) {
        _cacheNominal = FastHisto2D(*templHist, _conditional);
    } else {
        RooAbsPdf *pdf = dynamic_cast<RooAbsPdf *>(templ);
        const RooRealVar &x = dynamic_cast<const RooRealVar &>(_x.arg());
        const RooRealVar &y = dynamic_cast<const RooRealVar &>(_y.arg());
        const RooCmdArg &cond = _conditional ? RooFit::ConditionalObservables(RooArgSet(x)) : RooCmdArg::none();
        std::unique_ptr<TH1> hist{pdf->createHistogram("", x, RooFit::YVar(y), cond)};
        hist->SetDirectory(0); 
        _cacheNominal = FastHisto2D(dynamic_cast<TH2&>(*hist), _conditional);
    }
    if (_conditional) _cacheNominal.NormalizeXSlices(); 
    else              _cacheNominal.Normalize(); 
    if (_smoothAlgo < 0) {
        _cacheNominalLog = _cacheNominal;
        _cacheNominalLog.Log();
    }
    _cache = _cacheNominal;
}


void FastVerticalInterpHistPdf2::initComponent(int dim, TObject *thi, TObject *tlo) {
    FastHisto hi, lo; 
    TH1 *histHi = dynamic_cast<TH1*>(thi);
    TH1 *histLo = dynamic_cast<TH1*>(tlo);
    if (histHi && histLo) {
        hi = *histHi; lo = *histLo; 
    } else {
        RooAbsPdf *pdfHi = dynamic_cast<RooAbsPdf *>(thi);
        RooAbsPdf *pdfLo = dynamic_cast<RooAbsPdf *>(tlo);
        const RooRealVar &x = dynamic_cast<const RooRealVar &>(_x.arg());
        histHi = pdfHi->createHistogram("",x); histHi->SetDirectory(0); 
        histLo = pdfLo->createHistogram("",x); histLo->SetDirectory(0);
        hi = *histHi; lo = *histLo; 
        delete histHi; delete histLo; 
    }
    //printf("Un-normalized templates for dimension %d: \n", dim);  hi.Dump(); lo.Dump();
    hi.Normalize(); lo.Normalize();
    //printf("Normalized templates for dimension %d: \n", dim);  hi.Dump(); lo.Dump();
    initMorph(_morphs[dim], _cacheNominal, lo, hi);
}
void FastVerticalInterpHistPdf2D2::initComponent(int dim, TObject *thi, TObject *tlo) {
    FastHisto2D hi, lo; 
    TH2 *histHi = dynamic_cast<TH2*>(thi);
    TH2 *histLo = dynamic_cast<TH2*>(tlo);
    if (histHi && histLo) {
        hi = FastHisto2D(*histHi,_conditional); 
        lo = FastHisto2D(*histLo,_conditional); 
    } else {
        RooAbsPdf *pdfHi = dynamic_cast<RooAbsPdf *>(thi);
        RooAbsPdf *pdfLo = dynamic_cast<RooAbsPdf *>(tlo);
        const RooRealVar &x = dynamic_cast<const RooRealVar &>(_x.arg());
        const RooRealVar &y = dynamic_cast<const RooRealVar &>(_y.arg());
        const RooCmdArg &cond = _conditional ? RooFit::ConditionalObservables(RooArgSet(x)) : RooCmdArg::none();
        histHi =  dynamic_cast<TH2*>(pdfHi->createHistogram("", x, RooFit::YVar(y), cond)); histHi->SetDirectory(0);
        histLo =  dynamic_cast<TH2*>(pdfLo->createHistogram("", x, RooFit::YVar(y), cond)); histLo->SetDirectory(0);
        hi = FastHisto2D(*histHi,_conditional); 
        lo = FastHisto2D(*histLo,_conditional); 
        delete histHi; delete histLo; 
    }
    //printf("Un-normalized templates for dimension %d: \n", dim);  hi.Dump(); lo.Dump();
    if (_conditional) {
        hi.NormalizeXSlices(); lo.NormalizeXSlices();
    } else {
        hi.Normalize(); lo.Normalize();
    }
    //printf("Normalized templates for dimension %d: \n", dim);  hi.Dump(); lo.Dump();
    initMorph(_morphs[dim], _cacheNominal, lo, hi);
}



void FastVerticalInterpHistPdf2Base::initMorph(Morph &out, const FastTemplate &nominal, FastTemplate &lo, FastTemplate &hi) const {
    out.sum.Resize(hi.size());
    out.diff.Resize(hi.size());
    if (_smoothAlgo < 0)  {
        hi.LogRatio(nominal); lo.LogRatio(nominal);
        //printf("Log-ratios for dimension %d: \n", dim);  hi.Dump(); lo.Dump();
    } else {
        hi.Subtract(nominal); lo.Subtract(nominal);
        //printf("Differences for dimension %d: \n", dim);  hi.Dump(); lo.Dump();
    }
    FastTemplate::SumDiff(hi, lo, out.sum, out.diff);
    //printf("Sum and diff for dimension %d: \n", dim);  out.sum.Dump(); out.diff.Dump();
}




void FastVerticalInterpHistPdf2Base::syncTotal(FastTemplate &cache, const FastTemplate &cacheNominal, const FastTemplate &cacheNominalLog) const {
    /* === how the algorithm works, in theory ===
     * let  dhi = h_hi - h_nominal
     *      dlo = h_lo - h_nominal
     * and x be the morphing parameter
     * we define alpha = x * 0.5 * ((dhi-dlo) + (dhi+dlo)*smoothStepFunc(x));
     * which satisfies:
     *     alpha(0) = 0
     *     alpha(+1) = dhi 
     *     alpha(-1) = dlo
     *     alpha(x >= +1) = |x|*dhi
     *     alpha(x <= -1) = |x|*dlo
     *     alpha is continuous and has continuous first and second derivative, as smoothStepFunc has them
     * === and in practice ===
     * we already have computed the histogram for diff=(dhi-dlo) and sum=(dhi+dlo)
     * so we just do template += (0.5 * x) * (diff + smoothStepFunc(x) * sum)
     * ========================================== */

    // start from nominal
    cache.CopyValues(_smoothAlgo < 0 ? cacheNominalLog : cacheNominal);
    //printf("Cache initialized to nominal template: \n");  cacheNominal.Dump();

    // apply all morphs one by one
    for (int i = 0, ndim = _coefList.getSize(); i < ndim; ++i) {
        double x = _morphParams[i]->getVal();
        double a = 0.5*x, b = smoothStepFunc(x);
        cache.Meld(_morphs[i].diff, _morphs[i].sum, a, b);    
        //printf("Merged transformation for dimension %d, x = %+5.3f, step = %.3f: \n", i, x, b);  cache.Dump();
    }

    // if necessary go back to linear scale
    if (_smoothAlgo < 0) {
        cache.Exp();
        //printf("Done exponential tranformation\n");  cache.Dump();
    } else {
        cache.CropUnderflows();
    }
    
    // mark as done
    _sentry.reset();
}

void FastVerticalInterpHistPdf2::syncTotal() const {
    FastVerticalInterpHistPdf2Base::syncTotal(_cache, _cacheNominal, _cacheNominalLog);

    // normalize the result
    _cache.Normalize(); 
    //printf("Normalized result\n");  _cache.Dump();
}

void FastVerticalInterpHistPdf2D2::syncTotal() const {
    FastVerticalInterpHistPdf2Base::syncTotal(_cache, _cacheNominal, _cacheNominalLog);

    // normalize the result
    if (_conditional) _cache.NormalizeXSlices(); 
    else              _cache.Normalize(); 
    //printf("Normalized result\n");  _cache.Dump();
}

Int_t FastVerticalInterpHistPdf2D2::getMaxVal(const RooArgSet& vars) const {
    //static int ncalls = 0;
    //if (++ncalls < 100) {
    //    std::cout << "Called getMaxVal(" << GetName() << "), x  = " << _x.arg().GetName() << ", y = " << _y.arg().GetName() << ", conditional " << _conditional << std::endl;
    //    vars.Print("V");
    //}
    switch (vars.getSize()) {
        case 1:
            if (vars.contains(_x.arg())) return 1; 
            if (vars.contains(_y.arg())) return 2; 
            break;
        case 2:
            if (vars.contains(_x.arg()) && vars.contains(_y.arg())) {
                return 3;
            }
            break;
    }
    return 0;
}

Double_t FastVerticalInterpHistPdf2D2::maxVal(int code) const {
    if (!_initBase) initBase();
    if (_cache.size() == 0) _cache = _cacheNominal;
    if (!_sentry.good()) syncTotal();
    switch (code) {
        case 1:
            return _cache.GetMaxOnX(_y);
        case 2:
            return _cache.GetMaxOnY(_x);
        case 3:
            return _cache.GetMaxOnXY();
    }
    coutE(InputArguments) << "FastVerticalInterpHistPdf2D2::maxVal(" << GetName() 
			  << ") unsupported integration code " << code << "\n" << std::endl;
    assert(0);

    return 0.0;
}

FastVerticalInterpHistPdf2V::FastVerticalInterpHistPdf2V(const FastVerticalInterpHistPdf2 &hpdf, const RooAbsData &data, bool includeZeroWeights) :
    hpdf_(hpdf)
{
    // check init
    if (!hpdf._initBase) hpdf.initBase();
    if (hpdf._cache.size() == 0) hpdf._cache = hpdf._cacheNominal;
    if (!hpdf._sentry.good()) hpdf.syncTotal();
    // find bins
    std::vector<int> bins;
    RooArgSet obs(hpdf._x.arg());
    const RooRealVar &x = static_cast<const RooRealVar &>(*obs.first());
    bool aligned = true;
    for (int i = 0, n = data.numEntries(); i < n; ++i) {
        obs = *data.get(i);
        if (data.weight() == 0 && !includeZeroWeights) continue;
        int idx = hpdf._cache.FindBin(x.getVal());
        if (!bins.empty() && idx != bins.back() + 1) aligned = false;
        bins.push_back(idx);
    }
    if (bins.empty()) {
        // nothing to do.
    } else if (aligned) {
        begin_ = bins.front();
        end_   = bins.back()+1;
        //std::cout << "Created FastVerticalInterpHistPdf2V from " << hpdf.GetName() << ", aligned, " << (end_-begin_) << " bins." << std::endl;
    } else {
        nbins_ = bins.size();
        bins_.swap(bins);
        blocks_.clear();
        int start = bins_[0], istart = 0;
        for (int i = 1, n = bins_.size(); i < n; ++i) {
            if (bins_[i] != bins_[i-1]+1) { 
                blocks_.push_back(Block(istart,start,bins_[i-1]+1));
                start = bins_[i];
                istart = i;
            }
        }
        blocks_.push_back(Block(istart,start,bins_.back()+1));
        if (blocks_.size() < 4*bins_.size()) {
            //std::cout << "Created FastVerticalInterpHistPdf2V from " << hpdf.GetName() << ", block-aligned, " << bins_.size() << " bins, " <<  blocks_.size() << " blocks." << std::endl;
            bins_.clear();
        } else {
            //std::cout << "Created FastVerticalInterpHistPdf2V from " << hpdf.GetName() << ", non-aligned, " << bins_.size() << " bins, " <<  blocks_.size() << " blocks." << std::endl;
            blocks_.clear();
        }
    }
}

void FastVerticalInterpHistPdf2V::fill(std::vector<Double_t> &out) const 
{
    if (!hpdf_._sentry.good()) hpdf_.syncTotal();
    if (begin_ != end_) {
        out.resize(end_-begin_);
        std::copy(& hpdf_._cache.GetBinContent(begin_), (&hpdf_._cache.GetBinContent(end_-1))+1, out.begin());
    } else if (!blocks_.empty()) {
        out.resize(nbins_);
        for (auto b : blocks_) std::copy(& hpdf_._cache.GetBinContent(b.begin), (&hpdf_._cache.GetBinContent(b.end-1))+1, out.begin()+b.index);
    } else {
        out.resize(bins_.size());
        for (int i = 0, n = bins_.size(); i < n; ++i) {
          // if ((int)hpdf_._cache.GetNbinsX()>bins_[i]) out[i] = hpdf_._cache.GetBinContent(bins_[i]);
          //   else out[i] = 0;
            out[i] = hpdf_._cache.GetBinContent(bins_[i]);
        }
    }
}

