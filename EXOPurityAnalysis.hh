#ifndef EXOPurityAnalysis_hh_
#define EXOPurityAnalysis_hh_
#include "EXOAnalysisManager/EXOAnalysisModule.hh"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TMath.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

class EXOPurityAnalysis : public EXOAnalysisModule 
{
  public:
    EXOPurityAnalysis(EXOAnalysisManager* mgr) : EXOAnalysisModule(mgr) {}
    int Initialize();
    EventStatus ProcessEvent(EXOEventData *ED);
    int TalkTo(EXOTalkToManager *tm);
    int ShutDown();
    void SetOutputFile(std::string aval) {oName = aval;}
    void SetNBins(double aval) {NBINS = int(aval);}
    void SetDriftvelocity(double aval) {driftVelocity = aval;}
    void SetMinEntries(double aval) {minEntries = aval;}
    void SetSourceType(double aval) {SourceType = int(aval);}
    void SetSourcePosition(double aval) {SourcePosition = int(aval);}
    void SetMC(double aval) {MC = int(aval);}
    void DoAnalysis();
    void PrintResults();
    void GetEndPoint(TH1F *h, double *ep, double *perr, double *nerr);

  private:
    TH1F *hCRecon;
    TH1F *hCSumRecon;
    TH1F *hXRecon;
    TH1F *hYRecon;
    TH1F *hZRecon;
    TH2F *h2XYRecon;

    TH1F *hdtcl;
    TH1F *hPurAllCHPZ[100];
    TH1F *hPurAllCHNZ[100];

    TGraphAsymmErrors *grPurAllCHPZ;
    TGraphAsymmErrors *grPurAllCHNZ;

    double result_eLife_AllCHPZ;
    double result_eLife_err_AllCHPZ;
    double result_E0_AllCHPZ;
    double result_E0_err_AllCHPZ;
    double result_Chi2_AllCHPZ;

    double result_eLife_AllCHNZ;
    double result_eLife_err_AllCHNZ;
    double result_E0_AllCHNZ;
    double result_E0_err_AllCHNZ;
    double result_Chi2_AllCHNZ;

    int NBINS;
    double minEntries;
    double driftVelocity;
    int SourceType;
    int SourcePosition;
    int MC;

    TFile *oFile;
    std::string oName;

    int RunID;
    int StartTime;
    int EndTime;
    bool first;
    int MIN_CLUSTER_CUT;
    int MAX_CLUSTER_CUT;
    double R_CUT;
    int DT_MAX;
    int STP;

  DEFINE_EXO_ANALYSIS_MODULE( EXOPurityAnalysis )
}; 

#include "EXOUtilities/EXOPluginUtilities.hh"
EXO_DEFINE_PLUGIN_MODULE(EXOPurityAnalysis, "pam")
#endif
