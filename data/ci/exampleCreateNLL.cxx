#include <iostream>
#include <limits>
#include <memory>
#include <string>

#include "TFile.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"

#include "../../interface/Combine.h"
#include "../../interface/CascadeMinimizer.h"

namespace {
void printUsage(const char *argv0) {
  std::cerr << "Usage: " << argv0
            << " <workspace.root> [workspaceName=w] [modelConfig=ModelConfig] [dataName=data_obs]\n";
}
} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    printUsage(argv[0]);
    return 1;
  }

  const std::string fileName = argv[1];
  const std::string workspaceName = argc > 2 ? argv[2] : "w";
  const std::string modelConfigName = argc > 3 ? argv[3] : "ModelConfig";
  const std::string dataName = argc > 4 ? argv[4] : "data_obs";

  std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: failed to open input file '" << fileName << "'.\n";
    return 2;
  }

  auto *workspace = dynamic_cast<RooWorkspace *>(file->Get(workspaceName.c_str()));
  if (!workspace) {
    std::cerr << "ERROR: workspace '" << workspaceName << "' not found in '" << fileName << "'.\n";
    file->ls();
    return 2;
  }

  auto *modelConfig =
      dynamic_cast<RooStats::ModelConfig *>(workspace->genobj(modelConfigName.c_str()));
  if (!modelConfig) {
    std::cerr << "ERROR: ModelConfig '" << modelConfigName << "' not found in workspace '"
              << workspaceName << "'.\n";
    return 2;
  }

  RooAbsPdf *pdf = modelConfig->GetPdf();
  if (!pdf) {
    std::cerr << "ERROR: ModelConfig '" << modelConfigName << "' does not define a pdf.\n";
    return 2;
  }

  RooAbsData *data = workspace->data(dataName.c_str());
  if (!data) {
    std::cerr << "ERROR: dataset '" << dataName << "' not found in workspace '" << workspaceName
              << "'.\n";
    return 2;
  }

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);

  const RooArgSet *nuisances = modelConfig->GetNuisanceParameters();
  RooArgSet constraintSet;
  const RooArgSet *constraintPtr = nullptr;
  if (nuisances && nuisances->getSize() > 0) {
    constraintSet.add(*nuisances);
    constraintPtr = &constraintSet;
  }

  auto nll = combineCreateNLL(*pdf, *data, constraintPtr, /*offset=*/false);
  if (!nll) {
    std::cerr << "ERROR: combineCreateNLL returned a null pointer.\n";
    return 3;
  }

  const RooArgSet *poi = modelConfig->GetParametersOfInterest();
  RooRealVar *primaryPoi = nullptr;
  if (poi) {
    for (RooAbsArg *arg : *poi) {
      if (auto *var = dynamic_cast<RooRealVar *>(arg)) {
        var->setConstant(false);
        if (primaryPoi == nullptr)
          primaryPoi = var;
      }
    }
  }

  std::cout << "Initial NLL value: " << nll->getVal() << '\n';

  CascadeMinimizer minim(*nll, CascadeMinimizer::Constrained, primaryPoi);
  bool minimOk = minim.minimize(/*verbose=*/0);
  if (!minimOk) {
    std::cerr << "ERROR: minimization failed.\n";
    return 4;
  }

  minim.hesse(/*verbose=*/0);
  auto fitResult = std::unique_ptr<RooFitResult>(minim.save());
  std::cout << "Global minimum NLL: " << nll->getVal() << '\n';

  if (fitResult) {
    std::cout << "Minimizer status: " << fitResult->status()
              << ", edm=" << fitResult->edm() << '\n';
  }

  if (poi && poi->getSize() > 0) {
    std::cout << "Best-fit POI values:\n";
    for (RooAbsArg *arg : *poi) {
      auto *poiVar = dynamic_cast<RooRealVar *>(arg);
      if (!poiVar)
        continue;

      const double val = poiVar->getVal();
      double errHi = std::numeric_limits<double>::quiet_NaN();
      double errLo = std::numeric_limits<double>::quiet_NaN();
      if (fitResult) {
        if (auto *fitVar = dynamic_cast<RooRealVar *>(fitResult->floatParsFinal().find(poiVar->GetName()))) {
          errHi = fitVar->getErrorHi();
          errLo = fitVar->getErrorLo();
        }
      }

      std::cout << "  " << poiVar->GetName() << " = " << val;
      if (fitResult) {
        std::cout << " +" << errHi << " / " << errLo;
      }
      std::cout << '\n';
    }
  } else {
    std::cout << "ModelConfig has no parameters of interest.\n";
  }

  std::cout << "Finished building and evaluating the NLL from workspace '" << workspaceName
            << "'.\n";
  return 0;
}
