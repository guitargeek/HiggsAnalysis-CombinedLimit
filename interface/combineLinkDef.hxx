#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class std::vector<std::vector<Float_t>>+;
#pragma link C++ class std::vector<std::vector<std::vector<Float_t>>>+;
#pragma link C++ class std::vector<std::vector<std::vector<std::vector<Float_t>>>>+;
#pragma link C++ class std::vector<std::vector<std::vector<std::vector<std::vector<Float_t>>>>>+;
#pragma link C++ class std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<Float_t>>>>>>+;
#pragma link C++ class std::vector<std::vector<Float_t>>::iterator+;
#pragma link C++ class std::vector<std::vector<std::vector<Float_t>>>::iterator+;
#pragma link C++ class std::vector<std::vector<std::vector<std::vector<Float_t>>>>::iterator+;
#pragma link C++ class std::vector<std::vector<std::vector<std::vector<std::vector<Float_t>>>>>::iterator+;
#pragma link C++ class std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<Float_t>>>>>>::iterator+;
#pragma link C++ class std::vector<std::vector<Float_t>>::const_iterator+;
#pragma link C++ class std::vector<std::vector<std::vector<Float_t>>>::const_iterator+;
#pragma link C++ class std::vector<std::vector<std::vector<std::vector<Float_t>>>>::const_iterator+;
#pragma link C++ class std::vector<std::vector<std::vector<std::vector<std::vector<Float_t>>>>>::const_iterator+;
#pragma link C++ class std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<Float_t>>>>>>::const_iterator+;

#pragma link C++ class std::vector<std::vector<Double_t>>+;
#pragma link C++ class std::vector<std::vector<std::vector<Double_t>>>+;
#pragma link C++ class std::vector<std::vector<std::vector<std::vector<Double_t>>>>+;
#pragma link C++ class std::vector<std::vector<std::vector<std::vector<std::vector<Double_t>>>>>+;
#pragma link C++ class std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<Double_t>>>>>>+;
#pragma link C++ class std::vector<std::vector<Double_t>>::iterator+;
#pragma link C++ class std::vector<std::vector<std::vector<Double_t>>>::iterator+;
#pragma link C++ class std::vector<std::vector<std::vector<std::vector<Double_t>>>>::iterator+;
#pragma link C++ class std::vector<std::vector<std::vector<std::vector<std::vector<Double_t>>>>>::iterator+;
#pragma link C++ class std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<Double_t>>>>>>::iterator+;
#pragma link C++ class std::vector<std::vector<Double_t>>::const_iterator+;
#pragma link C++ class std::vector<std::vector<std::vector<Double_t>>>::const_iterator+;
#pragma link C++ class std::vector<std::vector<std::vector<std::vector<Double_t>>>>::const_iterator+;
#pragma link C++ class std::vector<std::vector<std::vector<std::vector<std::vector<Double_t>>>>>::const_iterator+;
#pragma link C++ class std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<Double_t>>>>>>::const_iterator+;

#pragma link C++ class std::map<int, std::vector<double>>+;
#pragma link C++ class std::map<int, std::pair<double,double>>+;
#pragma link C++ class std::map<std::vector<RooAbsReal*>, double>+;

#pragma link C++ class AsymPow+;

#pragma link C++ class FastTemplate+;
#pragma link C++ class FastHisto+;
#pragma link C++ class FastHisto2D+;
#pragma link C++ class FastHisto3D+;

#pragma link C++ class FastHistoAxis_d+;
#pragma link C++ class FastTemplate_d+;
#pragma link C++ class FastHisto_d+;
#pragma link C++ class FastHisto2D_d+;
#pragma link C++ class FastHisto3D_d+;

#pragma link C++ class FastHistoAxis_f+;
#pragma link C++ class FastTemplate_f+;
#pragma link C++ class FastHisto_f+;
#pragma link C++ class FastHisto2D_f+;
#pragma link C++ class FastHisto3D_f+;

#pragma link C++ class FastTemplateFunc_f+;
#pragma link C++ class FastHistoFunc_f+;
#pragma link C++ class FastHisto2DFunc_f+;
#pragma link C++ class FastHisto3DFunc_f+;

#pragma link C++ class FastTemplateFunc_d+;
#pragma link C++ class FastHistoFunc_d+;
#pragma link C++ class FastHisto2DFunc_d+;
#pragma link C++ class FastHisto3DFunc_d+;

#pragma link C++ class FastVerticalInterpHistPdf+;
#pragma link C++ class FastVerticalInterpHistPdf2+;
#pragma link C++ class FastVerticalInterpHistPdf2Base+;
#pragma link C++ class FastVerticalInterpHistPdf2D+;
#pragma link C++ class FastVerticalInterpHistPdf2D2+;
#pragma link C++ class FastVerticalInterpHistPdf3D+;
#pragma link C++ class FastVerticalInterpHistPdfBase+;
#pragma link C++ class FastVerticalInterpHistPdfBase::Morph+;
#pragma link C++ class HZZ4L_RooCTauPdf_1D+;
#pragma link C++ class HZZ4L_RooCTauPdf_1D_Expanded+;
#pragma link C++ class HZZ4L_RooCTauPdf_2D+;
#pragma link C++ class HZZ4L_RooSpinZeroPdf+;
#pragma link C++ class HZZ4L_RooSpinZeroPdf_1D+;
#pragma link C++ class HZZ4L_RooSpinZeroPdf_2D+;
#pragma link C++ class HZZ4L_RooSpinZeroPdf_phase+;
#pragma link C++ class VBFHZZ4L_RooSpinZeroPdf+;
#pragma link C++ class VBFHZZ4L_RooSpinZeroPdf_fast+;

#pragma link C++ class HZZ4L_RooSpinZeroPdf_1D_fast+;
#pragma link C++ class HZZ4L_RooSpinZeroPdf_2D_fast+;
#pragma link C++ class HZZ4L_RooSpinZeroPdf_phase_fast+;
#pragma link C++ class VVHZZ4L_RooSpinZeroPdf_1D_fast+;


#pragma link C++ class ProcessNormalization+;
#pragma link C++ class Roo2ExpPdf+;
#pragma link C++ class Roo4lMasses2D+;
#pragma link C++ class Roo4lMasses2D_Bkg+;
#pragma link C++ class Roo4lMasses2D_BkgGGZZ+;
#pragma link C++ class RooAlpha+;
#pragma link C++ class RooAlpha42ExpPdf+;
#pragma link C++ class RooAlpha4ErfPow2Pdf+;
#pragma link C++ class RooAlpha4ErfPowExpPdf+;
#pragma link C++ class RooAlpha4ErfPowPdf+;
#pragma link C++ class RooAlpha4ExpNPdf+;
#pragma link C++ class RooAlpha4ExpTailPdf+;
#pragma link C++ class RooAlpha4GausExpPdf+;
#pragma link C++ class RooAlphaExp+;
#pragma link C++ class RooAnaExpNPdf+;
#pragma link C++ class RooBWHighMassGGH+;
#pragma link C++ class RooBWRunPdf+;
#pragma link C++ class RooBernsteinFast<1>+;
#pragma link C++ class RooBernsteinFast<2>+;
#pragma link C++ class RooBernsteinFast<3>+;
#pragma link C++ class RooBernsteinFast<4>+;
#pragma link C++ class RooBernsteinFast<5>+;
#pragma link C++ class RooBernsteinFast<6>+;
#pragma link C++ class RooBernsteinFast<7>+;
#pragma link C++ class RooBetaFunc_v2+;
#pragma link C++ class RooCB+;
#pragma link C++ class RooCPSHighMassGGH+;
#pragma link C++ class RooCPSHighMassGGHNoInterf+;
#pragma link C++ class RooCPSHighMassVBF+;
#pragma link C++ class RooCPSHighMassVBFNoInterf+;
#pragma link C++ class RooChebyshevPDF+;
#pragma link C++ class RooDoubleCB+;
#pragma link C++ class RooDoubleCrystalBall+;
#pragma link C++ class RooErfExpPdf+;
#pragma link C++ class RooErfPdf+;
#pragma link C++ class RooErfPow2Pdf+;
#pragma link C++ class RooErfPowExpPdf+;
#pragma link C++ class RooErfPowPdf+;
#pragma link C++ class RooExpNPdf+;
#pragma link C++ class RooExpTailPdf+;
#pragma link C++ class RooFermi+;
#pragma link C++ class RooFourEMassRes+;
#pragma link C++ class RooFourEMassShapePdf2+;
#pragma link C++ class RooFourMuMassRes+;
#pragma link C++ class RooFourMuMassShapePdf2+;
#pragma link C++ class RooGausExpPdf+;
#pragma link C++ class RooGaussStepBernstein+;
#pragma link C++ class RooLevelledExp+;
#pragma link C++ class RooMultiPdf+;
#pragma link C++ class RooParamKeysPdf+;
#pragma link C++ class RooPow2Pdf+;
#pragma link C++ class RooPowPdf+;
#pragma link C++ class RooPowerExpPdf+;
#pragma link C++ class RooPowerFunction+;
#pragma link C++ class RooPowerLaw+;
#pragma link C++ class RooQCDPdf+;
#pragma link C++ class RooRelBW+;
#pragma link C++ class RooRelBW1+;
#pragma link C++ class RooRelBWHighMass+;
#pragma link C++ class RooRelBWUF+;
#pragma link C++ class RooRelBWUFParam+;
#pragma link C++ class RooRelBWUFParamWidth+;
#pragma link C++ class RooRelBWUF_SM4+;
#pragma link C++ class RooScaleHGamGamLOSM+;
#pragma link C++ class RooScaleHGamGamLOSMPlusX+;
#pragma link C++ class RooScaleHGluGluLOSM+;
#pragma link C++ class RooScaleHGluGluLOSMPlusX+;
#pragma link C++ class RooScaleLOSM+;
#pragma link C++ class RooSigPlusInt+;
#pragma link C++ class RooSimultaneousOpt+;
#pragma link C++ class RooSpline1D+;
#pragma link C++ class RooSplineND+;
#pragma link C++ class map<vector<int>,int>+;
#pragma link C++ class RooStarMomentMorph+;
#pragma link C++ class RooStats::HistFactory::RooBSpline+;
#pragma link C++ class RooStats::HistFactory::RooBSplineBases+;
#pragma link C++ class RooStepBernstein+;
#pragma link C++ class RooTH1DPdf+;
#pragma link C++ class RooTsallis+;
#pragma link C++ class RooTwoETwoMuMassRes+;
#pragma link C++ class RooTwoETwoMuMassShapePdf2+;
#pragma link C++ class RooUser1Pdf+;
#pragma link C++ class RooVBFZZPdf+;
#pragma link C++ class RooVBFZZPdf_v2+;
#pragma link C++ class RooaDoubleCBxBW+;
#pragma link C++ class RooggZZPdf+;
#pragma link C++ class RooggZZPdf_v2+;
#pragma link C++ class RooqqZZPdf+;
#pragma link C++ class RooqqZZPdf_v2+;
#pragma link C++ class SimpleCacheSentry+;
#pragma link C++ class SimpleGaussianConstraint+;
#pragma link C++ class SimplePoissonConstraint+;
#pragma link C++ class SimpleConstraintGroup+; // transient
#pragma link C++ class TH1Keys+;
#pragma link C++ class Triangle+;
#pragma link C++ class VerticalInterpHistPdf+;
#pragma link C++ class VerticalInterpPdf+;
#pragma link C++ class GaussExp+;
#pragma link C++ class RooParametricHist+;
#pragma link C++ class RooParametricHist2D+;
#pragma link C++ class RooParametricShapeBinPdf+;
#pragma link C++ class RooMorphingPdf+;
#pragma link C++ class CMSHistFunc+;
#pragma link C++ class CMSHistErrorPropagator+;
#pragma link C++ class CMSHistSum+;
#pragma link C++ class CMSHistFuncWrapper+;
#pragma link C++ class SimpleTaylorExpansion1D+;
#pragma link C++ class RooTaylorExpansion+;
#pragma link C++ class RooDoubleCBFast+;
#pragma link C++ class CombDataSetFactory+; // transient
#pragma link C++ class DebugProposal+; // transient
#pragma link C++ class cmsmath::SequentialMinimizer+; // transient
#pragma link C++ class TestProposal+; // transient
#pragma link C++ class RooCheapProduct+;
#pragma link C++ class CMSHggFormulaA1+;
#pragma link C++ class CMSHggFormulaA2+;
#pragma link C++ class CMSHggFormulaB1+;
#pragma link C++ class CMSHggFormulaB2+;
#pragma link C++ class CMSHggFormulaC1+;
#pragma link C++ class CMSHggFormulaD1+;
#pragma link C++ class CMSHggFormulaD2+;
#pragma link C++ class SimpleProdPdf+;
#pragma link C++ class CMSExternalMorph+;
#pragma link C++ class CMSInterferenceFunc+;
#pragma link C++ class RooEFTScalingFunction+;

#endif
