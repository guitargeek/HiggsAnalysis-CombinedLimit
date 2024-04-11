#include "../interface/VectorizedSimplePdfs.h"
#include "RooMath.h"
#include "vectorized.h"
#include <RooRealVar.h>
#include <stdexcept>
#include <memory>

VectorizedExponential::VectorizedExponential(const RooExponential &pdf, const RooAbsData &data, bool includeZeroWeights)
{
    RooArgSet obs(*data.get());
    std::unique_ptr<RooArgSet> params(pdf.getParameters(data));
    if (params->getSize() != 1) throw std::invalid_argument("Can't resolve which is the parameter of the exponential");

    x_ = dynamic_cast<const RooRealVar*>(obs.first());
    lambda_ = dynamic_cast<const RooAbsReal*>(params->first());

    xvals_.reserve(data.numEntries());
    for (unsigned int i = 0, n = data.numEntries(); i < n; ++i) {
        obs.assignValueOnly(*data.get(i), true);
        if (data.weight() || includeZeroWeights) xvals_.push_back(x_->getVal());        
    }
    work_.resize(xvals_.size());
}

void VectorizedExponential::fill(std::vector<Double_t> &out) const {
    Double_t xmax = x_->getMax(), xmin = x_->getMin(), lambda = lambda_->getVal();
    Double_t norm = (lambda != 0 ? (std::exp(lambda*xmax)-std::exp(lambda*xmin))/lambda  : xmax-xmin);
    out.resize(xvals_.size());
    vectorized::exponentials(xvals_.size(), lambda, norm, &xvals_[0], &out[0], &work_[0]);
}
