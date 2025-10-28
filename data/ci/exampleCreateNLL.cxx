#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

#include "TFile.h"

#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"

#include "../../interface/Combine.h"
#include "../../interface/CascadeMinimizer.h"
#include "../../interface/RooMultiPdfCombine.h"
#include "../../interface/utils.h"
#include "../../interface/ProfilingTools.h"

namespace {
void printUsage(const char *argv0) {
  std::cerr << "Usage: " << argv0
            << " <workspace.root> [workspaceName=w] [modelConfig=ModelConfig] [dataName=data_obs]\n";
}

bool loadSnapshotIfExists(RooWorkspace &ws, const char *name) {
  if (!name)
    return false;
  if (!ws.getSnapshot(name))
    return false;
  return ws.loadSnapshot(name);
}

void configureCascadeMinimizerState(RooWorkspace &ws, RooStats::ModelConfig &mc) {
  auto &cfg = CascadeMinimizerGlobalConfigs::O();

  cfg.parametersOfInterest = RooArgList();
  if (const RooArgSet *poi = mc.GetParametersOfInterest()) {
    for (RooAbsArg *arg : *poi) {
      auto *var = dynamic_cast<RooRealVar *>(arg);
      if (!var)
        continue;
      cfg.parametersOfInterest.add(*var);
    }
  }

  cfg.nuisanceParameters = RooArgList();
  if (const RooArgSet *nuis = mc.GetNuisanceParameters()) {
    for (RooAbsArg *arg : *nuis) {
      auto *var = dynamic_cast<RooRealVar *>(arg);
      if (!var)
        continue;
      if (!var->isConstant())
        cfg.nuisanceParameters.add(*var);
    }
  }

  cfg.allFloatingParameters = RooArgList();
  RooArgSet allVars(ws.allVars());
  for (RooAbsArg *arg : allVars) {
    auto *var = dynamic_cast<RooRealVar *>(arg);
    if (!var)
      continue;
    if (!var->isConstant())
      cfg.allFloatingParameters.add(*var);
  }

  cfg.pdfCategories = RooArgList();
  cfg.allRooMultiPdfParams = RooArgList();
  cfg.allRooMultiPdfs = RooArgList();

  std::unique_ptr<RooArgSet> discreteParameters{
      dynamic_cast<RooArgSet *>(ws.genobj("discreteParams"))};

  if (discreteParameters) {
    for (RooAbsArg *arg : *discreteParameters) {
      auto *cat = dynamic_cast<RooCategory *>(arg);
      if (!cat)
        continue;
      if (!cfg.pdfCategories.containsInstance(*cat))
        cfg.pdfCategories.add(*cat);
    }
  } else if (runtimedef::get("ADD_DISCRETE_FALLBACK")) {
    RooArgSet categories(ws.allCats());
    for (RooAbsArg *arg : categories) {
      auto *cat = dynamic_cast<RooCategory *>(arg);
      if (!cat)
        continue;
      if (std::string(cat->GetName()).find("pdfindex") == std::string::npos)
        continue;
      if (!cfg.pdfCategories.containsInstance(*cat))
        cfg.pdfCategories.add(*cat);
    }
  }

  if (cfg.pdfCategories.getSize() > 0) {
    RooArgSet clients;
    utils::getClients(cfg.pdfCategories, ws.allPdfs(), clients);
    for (RooAbsArg *arg : clients) {
      auto *multi = dynamic_cast<RooMultiPdf *>(arg);
      if (!multi)
        continue;
      if (!cfg.allRooMultiPdfs.containsInstance(*multi))
        cfg.allRooMultiPdfs.add(*multi);

      std::unique_ptr<RooArgSet> pdfParams{multi->getParameters((const RooArgSet *)nullptr)};
      for (RooAbsArg *a : *pdfParams) {
        auto *var = dynamic_cast<RooRealVar *>(a);
        if (!var || var->isConstant())
          continue;
        if (!cfg.allRooMultiPdfParams.containsInstance(*var))
          cfg.allRooMultiPdfParams.add(*var);
      }
    }
  }
}
}  // namespace

int main(int argc, char **argv) {
  CascadeMinimizer::initOptions();
  runtimedef::set("OPTIMIZE_BOUNDS", 1);
  runtimedef::set("ADDNLL_RECURSIVE", 1);
  runtimedef::set("ADDNLL_GAUSSNLL", 1);
  runtimedef::set("ADDNLL_HISTNLL", 1);
  runtimedef::set("ADDNLL_CBNLL", 1);
  runtimedef::set("TMCSO_AdaptivePseudoAsimov", 1);
  runtimedef::set("MINIMIZER_optimizeConst", 2);
  runtimedef::set("MINIMIZER_rooFitOffset", 1);
  runtimedef::set("ADDNLL_ROOREALSUM_FACTOR", 1);
  runtimedef::set("ADDNLL_ROOREALSUM_NONORM", 1);
  runtimedef::set("ADDNLL_ROOREALSUM_BASICINT", 1);
  runtimedef::set("ADDNLL_ROOREALSUM_KEEPZEROS", 1);
  runtimedef::set("ADDNLL_PRODNLL", 1);
  runtimedef::set("ADDNLL_HFNLL", 1);
  runtimedef::set("ADDNLL_HISTFUNCNLL", 1);
  runtimedef::set("ADDNLL_ROOREALSUM_CHEAPPROD", 1);

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

  utils::check_inf_parameters(workspace->allVars(), /*verbosity=*/0);
  loadSnapshotIfExists(*workspace, "clean");

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
  RooRealVar::enableSilentClipping();

  const RooArgSet *poi = modelConfig->GetParametersOfInterest();
  if (!poi || poi->getSize() == 0) {
    std::cerr << "ERROR: ModelConfig '" << modelConfigName << "' has no parameters of interest.\n";
    return 2;
  }
  for (RooAbsArg *arg : *poi) {
    auto *var = dynamic_cast<RooRealVar *>(arg);
    if (var)
      var->setConstant(false);
  }

  configureCascadeMinimizerState(*workspace, *modelConfig);

  const RooArgSet *nuisances = modelConfig->GetNuisanceParameters();
  RooArgSet constraintSet;
  const RooArgSet *constraintPtr = nullptr;
  if (nuisances && nuisances->getSize() > 0) {
    constraintSet.add(*nuisances);
    constraintPtr = &constraintSet;
  }

  auto nll = combineCreateNLL(*pdf, *data, constraintPtr, /*offset=*/true);
  if (!nll) {
    std::cerr << "ERROR: combineCreateNLL returned a null pointer.\n";
    return 3;
  }

  RooRealVar *primaryPoi = dynamic_cast<RooRealVar *>(poi->first());

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
  } else {
    std::cout << "WARNING: RooFitResult unavailable (minimizer did not provide one).\n";
  }

  std::cout << "Best-fit POI values:\n";
  for (RooAbsArg *arg : *poi) {
    auto *poiVar = dynamic_cast<RooRealVar *>(arg);
    if (!poiVar)
      continue;

    const double val = poiVar->getVal();
    double errHi = std::numeric_limits<double>::quiet_NaN();
    double errLo = std::numeric_limits<double>::quiet_NaN();

    if (fitResult) {
      if (auto *fitVar =
              dynamic_cast<RooRealVar *>(fitResult->floatParsFinal().find(poiVar->GetName()))) {
        errHi = fitVar->getAsymErrorHi();
        errLo = fitVar->getAsymErrorLo();
      }
    }

    std::cout << "  " << poiVar->GetName() << " = " << val;
    if (!std::isnan(errHi) && !std::isnan(errLo)) {
      std::cout << " +" << errHi << " / " << errLo;
    }
    std::cout << '\n';
  }

  std::cout << "Finished building and minimising the NLL from workspace '" << workspaceName
            << "'.\n";
  return 0;
}
