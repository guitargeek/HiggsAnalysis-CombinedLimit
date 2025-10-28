#include <iostream>
#include <memory>
#include <string>

#include "TFile.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"

#include "../../interface/Combine.h"

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

  std::cout << "Initial NLL value: " << nll->getVal() << '\n';

  const RooArgSet *poi = modelConfig->GetParametersOfInterest();
  if (poi && poi->getSize() > 0) {
    for (RooAbsArg *arg : *poi) {
      auto *poiVar = dynamic_cast<RooRealVar *>(arg);
      if (!poiVar)
        continue;

      const double originalValue = poiVar->getVal();
      std::cout << "  " << poiVar->GetName() << " = " << originalValue << " --> NLL "
                << nll->getVal() << '\n';

      if (poiVar->hasMin()) {
        poiVar->setVal(poiVar->getMin());
        std::cout << "  " << poiVar->GetName() << " at lower bound (" << poiVar->getMin()
                  << ") --> NLL " << nll->getVal() << '\n';
      }
      if (poiVar->hasMax()) {
        poiVar->setVal(poiVar->getMax());
        std::cout << "  " << poiVar->GetName() << " at upper bound (" << poiVar->getMax()
                  << ") --> NLL " << nll->getVal() << '\n';
      }

      poiVar->setVal(originalValue);
    }
  } else {
    std::cout << "ModelConfig has no parameters of interest to scan.\n";
  }

  std::cout << "Finished building and evaluating the NLL from workspace '" << workspaceName
            << "'.\n";
  return 0;
}
