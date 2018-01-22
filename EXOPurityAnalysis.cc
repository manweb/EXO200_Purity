#include "EXOPurityAnalysis.hh"
#include "EXOUtilities/EXOErrorLogger.hh"
#include "EXOUtilities/EXOTalkToManager.hh"
#include "EXOUtilities/EXOEventData.hh"
#include "EXOUtilities/EXOChannelMap.hh"
#include "EXOUtilities/EXODimensions.hh"
#include "EXOUtilities/SystemOfUnits.hh"
#include "EXOAnalysisManager/EXOAnalysisManager.hh"
#include "EXOCalibUtilities/EXOCalibManager.hh"
#include <iostream>

int EXOPurityAnalysis::Initialize()
{
  std::cout << "Initializing module pam..." << std::endl;

  MIN_CLUSTER_CUT = 0;
  MAX_CLUSTER_CUT = 100;
  R_CUT = 200;

  int nBins = 0;
  double hMax = 0.0;
  if (SourceType == 0) {nBins = 250; hMax = 2500;}
  else {nBins = 360; hMax = 3600;}

  double V_D = driftVelocity * 1000.0; // [mm/us]
  double dt_maxTMP = CATHODE_ANODE_x_DISTANCE / V_D;
  DT_MAX = TMath::Nint(dt_maxTMP / NBINS) * NBINS;

  hCRecon = new TH1F("hCRecon","Recon cluster energy",300,0,3000);
  hCSumRecon = new TH1F("hSumRecon","Sum recon energy",2000,0,4000);
  hXRecon = new TH1F("hXRecon","x recon cluster",200,-200,200);
  hYRecon = new TH1F("hYRecon","y recon cluster",200,-200,200);
  hZRecon = new TH1F("hZRecon","z recon cluster",200,-200,200);
  h2XYRecon = new TH2F("h2XYRecon","x-y recon cluster",200,-200,200,200,-200,200);
  hdtcl = new TH1F("hdtcl","Drift time",200,0,200000);

  STP = TMath::Nint(dt_maxTMP / NBINS);
  // initialize purity histograms
  for (int i = 0; i < NBINS; i++) {
     char hNamePZ[50];
     char hNameNZ[50];
     char hTitlePZ[50];
     char hTitleNZ[50];

     sprintf(hNamePZ,"h%iPurPZ",i);
     sprintf(hNameNZ,"h%iPurNZ",i);
     sprintf(hTitlePZ,"Purity histo (%i-%ius) PZ, all channels",i*STP,(i+1)*STP);
     sprintf(hTitleNZ,"Purity histo (%i-%ius) NZ, all channels",i*STP,(i+1)*STP);

     hPurAllCHPZ[i] = new TH1F(hNamePZ,hTitlePZ,nBins,0,hMax);
     hPurAllCHNZ[i] = new TH1F(hNameNZ,hTitleNZ,nBins,0,hMax);
  }

  std::cout << "EXOPurityAnalysis-> Number of z-bins = " << NBINS << std::endl;
  std::cout << "EXOPurityAnalysis-> Set drift velocity = " << driftVelocity << " mm/ns" << std::endl;
  std::cout << "EXOPurityAnalysis-> Set minimum entries = " << minEntries << std::endl;
  std::cout << "EXOPurityAnalysis-> Source type: ";
  if (SourceType == 0) {std::cout << "Co60";}
  if (SourceType == 1) {std::cout << "Th228";}
  std::cout << std::endl;
  std::cout << "EXOPurityAnalysis-> Source position: ";
  if (SourcePosition == 0) {std::cout << "Cathode";}
  if (SourcePosition == 1) {std::cout << "Anode";}
  std::cout << std::endl;
  std::cout << "EXOPurityAnalysis-> MC = " << MC << std::endl;

  first = true;

  return 0;
}

EXOAnalysisModule::EventStatus EXOPurityAnalysis::ProcessEvent(EXOEventData *ED)
{
  if (first) {
     StartTime = ED->fEventHeader.fTriggerSeconds;
     RunID = ED->fRunNumber;
     first = false;
  }
  EndTime = ED->fEventHeader.fTriggerSeconds;

// **** Fill purity histograms ********************************************************************

  int ncl = ED->GetNumChargeClusters();

  // cut on number of clusters
  if (ncl < MIN_CLUSTER_CUT) {return kOk;}
  if (ncl > MAX_CLUSTER_CUT) {return kOk;}

  double sumReconCharge = 0.0;
  for (int i = 0; i < ncl; i++) {
     double cX = ED->GetChargeCluster(i)->fX;
     double cY = ED->GetChargeCluster(i)->fY;
     double cZ = 0.0;
     double dt = 0.0;
     double eccl = 0.0;

     // *** This is used for simulated data **********
     if (MC == 1) {
        cZ = 0.0;
        if (ED->GetChargeCluster(i)->fCollectionTime < 1024000) {continue;}
        dt = ED->GetChargeCluster(i)->fCollectionTime / 1000.0 - 1024.0;
        eccl = ED->GetChargeCluster(i)->fCorrectedEnergy;
        //eccl = ED->GetChargeCluster(i)->fCorrectedEnergy * TMath::Exp(-dt/250.0);
        if (eccl < 57) {continue;}
     }
     // **********************************************

     // *** This is used for real data ***************
     else {
        cZ = ED->GetChargeCluster(i)->fZ;
        dt = ED->GetChargeCluster(i)->fDriftTime / 1000.0;
        eccl = ED->GetChargeCluster(i)->fCorrectedEnergy;
     }
     // **********************************************

     // keep only clusters in the Xenon
     if (cX < -200 || cX > 200 || cY < -200 || cY > 200 || cZ < -200 || cZ > 200) {continue;}

     // position cut in X-Y plane
     double R = TMath::Sqrt(cX*cX + cY*cY);
     if (R > R_CUT) {continue;}

     sumReconCharge += ED->GetChargeCluster(i)->fCorrectedEnergy;
     hCRecon->Fill(ED->GetChargeCluster(i)->fCorrectedEnergy);
     hXRecon->Fill(cX);
     hYRecon->Fill(cY);
     hZRecon->Fill(cZ);
     h2XYRecon->Fill(cX,cY);
     hdtcl->Fill(ED->GetChargeCluster(i)->fDriftTime);

     double PurIDTMP = floor(dt / STP);
     int PurID = int(PurIDTMP);

     if (PurID < 0 || PurID > NBINS-1) {continue;}
     //std::cout << "BinID = " << PurID << "  dt = " << dt << "  eccl = " << eccl << std::endl;

     if (ED->GetChargeCluster(i)->fDetectorHalf > 1) {continue;}
     if (ED->GetChargeCluster(i)->fDetectorHalf == 0) {hPurAllCHPZ[PurID]->Fill(eccl);}
     if (ED->GetChargeCluster(i)->fDetectorHalf == 1) {hPurAllCHNZ[PurID]->Fill(eccl);}
  }

  hCSumRecon->Fill(sumReconCharge);

// **** End filling histograms ********************************************************************

  return kOk;
}

int EXOPurityAnalysis::TalkTo(EXOTalkToManager *tm)
{
  tm->CreateCommand("/pam/file","name of output file",this,"output.root",&EXOPurityAnalysis::SetOutputFile);
  if ( oName == NULL ) {
    errorLogger->LogError("EXOPurityAnalysis","TalkTo","unable to create output file command",5);
  }

  tm->CreateCommand("/pam/NBins","Number of z-bins",this,40.0,&EXOPurityAnalysis::SetNBins);
  if ( NBINS == 0.0 ) {
    errorLogger->LogError("EXOPurityAnalysis","TalkTo","number of z-bins must be > 0",5);
  }
  if ( NBINS == 100.0 ) {
    errorLogger->LogError("EXOPurityAnalysis","TalkTo","number of z-bins should not exceed 100",5);
  }

  tm->CreateCommand("/pam/driftVelocity","Drift velocity in mm/ns",this,0.0018,&EXOPurityAnalysis::SetDriftvelocity);
  if ( driftVelocity == 0.0 ) {
    errorLogger->LogError("EXOPurityAnalysis","TalkTo","unable to set drift velocity",5);
  }

  tm->CreateCommand("/pam/minEntries","Minimum number of entries in histograms for individual channls",this,400.0,&EXOPurityAnalysis::SetMinEntries);
  if ( minEntries == 0.0 ) {
    errorLogger->LogError("EXOPurityAnalysis","TalkTo","unable to set mimimum number of histogram entries",5);
  }

  tm->CreateCommand("/pam/sourceType","Source type: 0 = Co60, 1 = Th228",this,0.0,&EXOPurityAnalysis::SetSourceType);

  tm->CreateCommand("/pam/sourcePosition","Position of the source: 0 = cathode, 1 = anode",this,0.0,&EXOPurityAnalysis::SetSourcePosition);

  tm->CreateCommand("/pam/MC","Set whether the data is MC or real data",this,0.0,&EXOPurityAnalysis::SetMC);

  std::cout << "at talk manager" << std::endl;

  return 0;
}

int EXOPurityAnalysis::ShutDown()
{
  DoAnalysis();

  PrintResults();

  return 0;
}

void EXOPurityAnalysis::DoAnalysis()
{
  //if (MC == 0) {
     int nBins = hPurAllCHPZ[0]->GetNbinsX();
     // zero out bins below threshold
     for (int i = 0; i < NBINS; i++) {
        for (int k = 0; k < nBins; k++) {
           double BinContPZ = hPurAllCHPZ[i]->GetBinContent(k+1);
           double BinContNZ = hPurAllCHNZ[i]->GetBinContent(k+1);

           if (BinContPZ < 2) {hPurAllCHPZ[i]->SetBinContent(k+1,0);}
           if (BinContNZ < 2) {hPurAllCHNZ[i]->SetBinContent(k+1,0);}
        }
     }
  //}

  //int stp = DT_MAX / 20;
  double x1[NBINS];
  double y1[NBINS];
  double exl1[NBINS];
  double exh1[NBINS];
  double eyl1[NBINS];
  double eyh1[NBINS];

  double x2[NBINS];
  double y2[NBINS];
  double exl2[NBINS];
  double exh2[NBINS];
  double eyl2[NBINS];
  double eyh2[NBINS];

  int nPointsPZ = 0;
  int nPointsNZ = 0;

  for (int i = 0; i < NBINS; i++) {
     if (hPurAllCHPZ[i]->GetEntries() < minEntries) {continue;}

     double EP1 = 0.0;
     double EP_PERR1 = 0.0;
     double EP_NERR1 = 0.0;

     GetEndPoint(hPurAllCHPZ[i], &EP1, &EP_PERR1, &EP_NERR1);

     //if (EP_PERR1 == 0) {EP_PERR1 = 100;}
     //if (EP_NERR1 == 0) {EP_NERR1 = 100;}

     x1[nPointsPZ] = (i+0.5)*double(STP);
     y1[nPointsPZ] = EP1;
     exl1[nPointsPZ] = double(STP) / 2.0;
     exh1[nPointsPZ] = double(STP) / 2.0;
     eyl1[nPointsPZ] = EP_NERR1;
     eyh1[nPointsPZ] = EP_PERR1;

     nPointsPZ++;
  }

  for (int i = 0; i < NBINS; i++) {
     if (hPurAllCHNZ[i]->GetEntries() < minEntries) {continue;}

     double EP2 = 0.0;
     double EP_PERR2 = 0.0;
     double EP_NERR2 = 0.0;

     GetEndPoint(hPurAllCHNZ[i], &EP2, &EP_PERR2, &EP_NERR2);

     //if (EP_PERR2 == 0) {EP_PERR2 = 100;}
     //if (EP_NERR2 == 0) {EP_NERR2 = 100;}

     x2[nPointsNZ] = (i+0.5)*double(STP);
     y2[nPointsNZ] = EP2;
     exl2[nPointsNZ] = double(STP) / 2.0;
     exh2[nPointsNZ] = double(STP) / 2.0;
     eyl2[nPointsNZ] = EP_NERR2;
     eyh2[nPointsNZ] = EP_PERR2;

     nPointsNZ++;
  }

  grPurAllCHPZ = new TGraphAsymmErrors(nPointsPZ,x1,y1,exl1,exh1,eyl1,eyh1);
  grPurAllCHNZ = new TGraphAsymmErrors(nPointsNZ,x2,y2,exl2,exh2,eyl2,eyh2);

  grPurAllCHPZ->SetTitle("Purity PZ");
  grPurAllCHNZ->SetTitle("Purity NZ");

  TF1 *fitPZ = new TF1("fitPZ","[0]*TMath::Exp(-x/[1])");
  TF1 *fitNZ = new TF1("fitNZ","[0]*TMath::Exp(-x/[1])");

  /*if (SourcePosition == 0) {
     fitPZ->SetRange(double(DT_MAX) / 2.0, DT_MAX - STP);
     fitNZ->SetRange(double(DT_MAX) / 2.0, DT_MAX - STP);
  }
  else {
     fitPZ->SetRange(STP,double(DT_MAX) / 2.0);
     fitNZ->SetRange(STP,double(DT_MAX) / 2.0);
  }*/

  fitPZ->SetRange(double(STP), double(DT_MAX - 2*STP));
  fitNZ->SetRange(double(STP), double(DT_MAX - 2*STP));

  fitPZ->SetParameters(3000,240);
  fitNZ->SetParameters(3000,240);

  fitPZ->SetParLimits(0,0,100000);
  fitPZ->SetParLimits(1,100,400);

  fitNZ->SetParLimits(0,0,100000);
  fitNZ->SetParLimits(1,100,400);

  fitPZ->SetLineWidth(1);
  fitPZ->SetLineColor(kRed);

  fitNZ->SetLineWidth(1);
  fitNZ->SetLineColor(kRed);

  grPurAllCHPZ->Fit(fitPZ,"rq");
  grPurAllCHNZ->Fit(fitNZ,"rq");

  result_eLife_AllCHPZ = fitPZ->GetParameter(1);
  result_eLife_err_AllCHPZ = fitPZ->GetParError(1);
  result_E0_AllCHPZ = fitPZ->GetParameter(0);
  result_E0_err_AllCHPZ = fitPZ->GetParError(0);

  result_eLife_AllCHNZ = fitNZ->GetParameter(1);
  result_eLife_err_AllCHNZ = fitNZ->GetParError(1);
  result_E0_AllCHNZ = fitNZ->GetParameter(0);
  result_E0_err_AllCHNZ = fitNZ->GetParError(0);

  result_Chi2_AllCHPZ = fitPZ->GetChisquare() / fitPZ->GetNDF();
  result_Chi2_AllCHNZ = fitNZ->GetChisquare() / fitNZ->GetNDF();

  return;
}

void EXOPurityAnalysis::PrintResults()
{
  std::cout << "EXOPurityAnalysis-> Opening the output file to write histograms..." << std::endl;
  oFile = new TFile(oName.c_str(),"RECREATE");
  std::cout << "EXOPurityAnalysis-> Writing histograms..." << std::endl;

  hCSumRecon->Write();
  hCRecon->Write();
  hXRecon->Write();
  hYRecon->Write();
  hZRecon->Write();
  h2XYRecon->Write();
  hdtcl->Write();

  for (int i = 0; i < NBINS; i++) {
     hPurAllCHPZ[i]->Write();
     hPurAllCHNZ[i]->Write();
  }

  TCanvas *c1 = new TCanvas("c1","Purity PZ (combined)");
  grPurAllCHPZ->Draw("A*");

  TCanvas *c2 = new TCanvas("c2","Purity NZ (combined)");
  grPurAllCHNZ->Draw("A*");

  c1->Write();
  c2->Write();

  std::cout << "EXOPurityAnalysis-> Closing file..." << std::endl;
  oFile->Close();

  std::cout << "EXOPurityAnalysis-> Printing results..." << std::endl;

  std::cout << "**** Result for combined channels ********************************************************************************************" << std::endl;

  std::cout << "Summary electron lifetime:" << std::endl;
  std::cout << "\tPZ: " << result_eLife_AllCHPZ << " +- " << result_eLife_err_AllCHPZ << std::endl;
  std::cout << "\tNZ: " << result_eLife_AllCHNZ << " +- " << result_eLife_err_AllCHNZ << std::endl;

  std::cout << RunID << "\t" << ((StartTime - 1304146800)+(EndTime-StartTime)/2.0)/3600.0/24.0 << "\t" << (result_eLife_AllCHPZ+result_eLife_AllCHNZ)/2.0 << "\t" << TMath::Sqrt(result_eLife_err_AllCHPZ*result_eLife_err_AllCHPZ + result_eLife_err_AllCHNZ*result_eLife_err_AllCHNZ)/2.0 << "\t" << SourceType << "\t" << SourcePosition << std::endl;

  if (SourceType == 0) {std::cout << "Co60 run";}
  else {std::cout << "Th228 run";}

  if (SourcePosition == 0) {std::cout << " at cathode" << std::endl;}
  else {std::cout << " at anode" << std::endl;}

  return;
}

void EXOPurityAnalysis::GetEndPoint(TH1F *h, double *ep, double *perr, double *nerr)
{
  int integral = int(h->Integral());

  int BinSum1 = 0;
  int BinSum3 = 0;

  int nBins = h->GetNbinsX();

  int i1 = 0;
  int i2 = 0;
  int i3 = 0;

  // get endpoint (int > 99.9%)
  for (i1 = 1; i1 <= nBins; i1++) {
     BinSum1 += int(h->GetBinContent(i1));
     if (double(BinSum1) / double(integral) >= 0.999) {break;}
  }

  int BinSum2 = BinSum1;

  // get error (int = 100%)
  for (i2 = i1; i2 <= nBins; i2++) {
     BinSum2 += int(h->GetBinContent(i2));
     if (double(BinSum2) / double(integral) >= 1.0) {break;}
  }

  // get error (int >= 99.8%)
  for (i3 = 1; i3 <= nBins; i3++) {
     BinSum3 += int(h->GetBinContent(i3));
     if (double(BinSum3) / double(integral) >= 0.998) {break;}
  }

  double BinWidth = h->GetBinWidth(1);
  double EndPoint = i1 * BinWidth;
  double PErr = (i2-i1) * BinWidth;
  double NErr = (i1 - i3) * BinWidth;

  TF1 *fit = new TF1("fit","[0]*TMath::Erfc((x-[1])/[2])",1300,3600);

  fit->SetParameters(40, EndPoint, 100);
  h->Fit("fit","rq");

  double par[3];
  fit->GetParameters(par);

  double *parErrors = fit->GetParErrors();

  //if (EndPointMethod == 0) {
     *ep = EndPoint;
     *perr = PErr;
     *nerr = NErr;
  //}
  /*if (EndPointMethod == 1) {
     *ep = par[1] + par[2];
     *perr = TMath::Sqrt(parErrors[1]*parErrors[1] + parErrors[2]*parErrors[2]);
     *nerr = TMath::Sqrt(parErrors[1]*parErrors[1] + parErrors[2]*parErrors[2]);
  }*/

  return;
}
