#include "../interface/AsymQuad.h"

#include <cmath>
#include <cassert>
#include <cstdio>

AsymQuad::AsymQuad() :
RooAbsReal(),
_funcList("funcList", "List of functions", this),
_coefList("coefList", "List of coefficients", this),
smoothRegion_(0),
smoothAlgo_(0)
{
}

AsymQuad::AsymQuad(const char *name, const char *title, const RooArgList& inFuncList, const RooArgList& inCoefList, Double_t smoothRegion, Int_t smoothAlgo) :
RooAbsReal(name, title),
_funcList("funcList", "List of functions", this),
_coefList("coefList", "List of coefficients", this),
smoothRegion_(smoothRegion),
smoothAlgo_(smoothAlgo)
{
  if (inFuncList.getSize()!=2*inCoefList.getSize()+1) {
    coutE(InputArguments) << "AsymQuad::AsymQuad(" << GetName()
      << ") number of functions and coefficients inconsistent, must have Nfunc=1+2*Ncoef" << std::endl;
    assert(0);
  }

  for (RooAbsArg* func : inFuncList) {
    if (!dynamic_cast<RooAbsReal*>(func)) {
      coutE(InputArguments) << "ERROR: AsymQuad::AsymQuad(" << GetName() << ") function  " << func->GetName() << " is not of type RooAbsReal" << std::endl;
      assert(0);
    }
    _funcList.add(*func);
  }

  for (RooAbsArg* coef : inCoefList) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      coutE(InputArguments) << "ERROR: AsymQuad::AsymQuad(" << GetName() << ") coefficient " << coef->GetName() << " is not of type RooAbsReal" << std::endl;
      assert(0);
    }
    _coefList.add(*coef);
  }
}

AsymQuad::AsymQuad(const AsymQuad& other, const char* name):
RooAbsReal(other, name),
_funcList("!funcList", this, other._funcList),
_coefList("!coefList", this, other._coefList),
smoothRegion_(other.smoothRegion_),
smoothAlgo_(other.smoothAlgo_)
{
}

AsymQuad::~AsymQuad() = default;

Double_t AsymQuad::evaluate() const {

  RooAbsReal* func = &(RooAbsReal&)_funcList[0];

  Double_t central = func->getVal();
  Double_t result = central;

  for (int iCoef = 0; iCoef < _coefList.getSize(); ++iCoef) {
    Double_t coefVal = static_cast<RooAbsReal&>(_coefList[iCoef]).getVal();
    RooAbsReal* funcUp = &(RooAbsReal&)_funcList[2 * iCoef + 1];
    RooAbsReal* funcDn = &(RooAbsReal&)_funcList[2 * iCoef + 2];
    result += interpolate(coefVal, central, funcUp->getVal(), funcDn->getVal());
  }

  return result;
}

Double_t AsymQuad::interpolate(Double_t theta_, Double_t valueCenter_, Double_t valueHigh_, Double_t valueLow_) const {
  if (smoothAlgo_<0) return 0;
  else{
    if (fabs(theta_)>=smoothRegion_) return theta_ * (theta_ > 0 ? valueHigh_ - valueCenter_ : valueCenter_ - valueLow_);
    if (smoothAlgo_ == 0) {
      // Quadratic interpolation null at zero and continuous at boundaries but not smooth at boundaries
      Double_t c_up  = +theta_ * (smoothRegion_ + theta_) / (2 * smoothRegion_);
      Double_t c_dn  = -theta_ * (smoothRegion_ - theta_) / (2 * smoothRegion_);
      Double_t c_cen = -theta_ * theta_ / smoothRegion_;
      return c_up * valueHigh_ + c_dn * valueLow_ + c_cen * valueCenter_;
    }
    else if (smoothAlgo_ == 1){
      // Quadratic interpolation that is everywhere differentiable but not null at zero
      Double_t c_up  = (smoothRegion_ + theta_) * (smoothRegion_ + theta_) / (4 * smoothRegion_);
      Double_t c_dn  = (smoothRegion_ - theta_) * (smoothRegion_ - theta_) / (4 * smoothRegion_);
      Double_t c_cen = -c_up - c_dn;
      return c_up * valueHigh_ + c_dn * valueLow_ + c_cen * valueCenter_;
    }
    else/* if (smoothAlgo_ == 2)*/{
      // Quadratic interpolation that is everywhere differentiable and null at zero
      Double_t cnorm = theta_/smoothRegion_;
      Double_t cnorm2 = pow(cnorm, 2);
      Double_t hi = valueHigh_ - valueCenter_;
      Double_t lo = valueLow_ - valueCenter_;
      Double_t sum = hi+lo;
      Double_t diff = hi-lo;
      Double_t a = theta_/2.; // cnorm*smoothRegion_
      Double_t b = 0.125 * cnorm * (cnorm2 * (3.*cnorm2 - 10.) + 15.);
      Double_t result = a*(diff + b*sum);
      return result;
    }
  }
}

ClassImp(AsymQuad)
