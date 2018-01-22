void Plot()
{
  TTree *t1 = new TTree;
  TTree *t2 = new TTree;
  TTree *t3 = new TTree;

  t1->ReadFile("StrongTh228/StrongTh228.dat");
  t2->ReadFile("StrongCo60/StrongCo60.dat");
  t3->ReadFile("WeakTh228/WeakTh228.dat");
  //t3->ReadFile("StrongTh228/StrongTh228_2.dat");

  double Time1;
  double Time2;
  double Time3;
  double ELife1;
  double ELife2;
  double ELife3;
  double ELifeErr1;
  double ELifeErr2;
  double ELifeErr3;
  int Type1;
  int Type2;
  int Type3;
  int Pos1;
  int Pos2;
  int Pos3;

  t1->SetBranchAddress("Time",&Time1);
  t1->SetBranchAddress("ELife",&ELife1);
  t1->SetBranchAddress("ELifeErr",&ELifeErr1);
  t1->SetBranchAddress("Type",&Type1);
  t1->SetBranchAddress("Pos",&Pos1);

  t2->SetBranchAddress("Time",&Time2);
  t2->SetBranchAddress("ELife",&ELife2);
  t2->SetBranchAddress("ELifeErr",&ELifeErr2);
  t2->SetBranchAddress("Type",&Type2);
  t2->SetBranchAddress("Pos",&Pos2);

  t3->SetBranchAddress("Time",&Time3);
  t3->SetBranchAddress("ELife",&ELife3);
  t3->SetBranchAddress("ELifeErr",&ELifeErr3);
  t3->SetBranchAddress("Type",&Type3);
  t3->SetBranchAddress("Pos",&Pos3);

  TMultiGraph *mg1 = new TMultiGraph();
  TMultiGraph *mg2 = new TMultiGraph();

  TGraphErrors *gr1C = new TGraphErrors();
  TGraphErrors *gr1A = new TGraphErrors();
  TGraphErrors *gr2C = new TGraphErrors();
  TGraphErrors *gr2A = new TGraphErrors();
  TGraphErrors *gr3C = new TGraphErrors();
  TGraphErrors *gr3A = new TGraphErrors();

  //gr1C->GetYaxis()->SetRangeUser(0,320);
  mg1->SetMinimum(0);
  mg1->SetMaximum(320);

  gr1C->SetMarkerStyle(26);
  gr1A->SetMarkerStyle(26);
  gr1C->SetMarkerSize(0.5);
  gr1A->SetMarkerSize(0.5);
  gr1C->SetMarkerColor(kBlack);
  gr1A->SetMarkerColor(kRed);

  gr2C->SetMarkerStyle(22);
  gr2A->SetMarkerStyle(22);
  gr2C->SetMarkerSize(0.8);
  gr2A->SetMarkerSize(0.8);
  gr2C->SetMarkerColor(kBlack);
  gr2A->SetMarkerColor(kRed);

  gr3C->SetMarkerStyle(4);
  gr3A->SetMarkerStyle(4);
  gr3C->SetMarkerSize(0.8);
  gr3A->SetMarkerSize(0.8);
  gr3C->SetMarkerColor(kBlack);
  gr3A->SetMarkerColor(kRed);

  mg1->Add(gr1C);
  mg2->Add(gr1A);
  mg1->Add(gr2C);
  mg2->Add(gr2A);
  mg1->Add(gr3C);

  TLegend *l = new TLegend(0.0,0.0,0.2,0.2);
  l->AddEntry(gr1C,"Strong Th228 (cathode)");
  l->AddEntry(gr1A,"Strong Th228 (anode)");
  l->AddEntry(gr2C,"Strong Co60 (cathode)");
  l->AddEntry(gr2A,"Strong Co60 (anode)");
  l->AddEntry(gr3C,"Weak Th228 (cathode)");
  l->AddEntry(gr3A,"Weak Th228 (anode)");

  int n1C = 0;
  int n1A = 0;
  int n2C = 0;
  int n2A = 0;
  int n3C = 0;
  int n3A = 0;
  for (int i = 0; i < t1->GetEntries(); i++) {
     t1->GetEntry(i);

     if (Pos1 == 0) {gr1C->SetPoint(n1C,Time1,ELife1); gr1C->SetPointError(n1C,0,ELifeErr1); n1C++;}
     else {gr1A->SetPoint(n1A,Time1,ELife1); gr1A->SetPointError(n1A,0,ELifeErr1); n1A++;}
  }

  for (int i = 0; i < t2->GetEntries(); i++) {
     t2->GetEntry(i);
  
     if (Pos2 == 0) {gr2C->SetPoint(n2C,Time2,ELife2); gr2C->SetPointError(n2C,0,ELifeErr2); n2C++;}
     else {gr2A->SetPoint(n2A,Time2,ELife2); gr2A->SetPointError(n2A,0,ELifeErr2); n2A++;}
  }

  for (int i = 0; i < t3->GetEntries(); i++) {
     t3->GetEntry(i);

     if (Pos3 == 0) {gr3C->SetPoint(n3C,Time3,ELife3); gr3C->SetPointError(n3C,0,ELifeErr3); n3C++;}
     else {gr3A->SetPoint(n3A,Time3,ELife3); gr3A->SetPointError(n3A,0,ELifeErr3); n3A++;}
  }

  TF1 *fit = new TF1("fit","pol4",0,80);
  fit->SetLineWidth(1);

  mg1->Fit("fit","r");

  TF1 *func = new TF1("func","pol4",0,80);
  func->SetParameters(-532.7,75.99,-2.568,0.0353,-0.0001679);
  func->SetLineColor(kBlue);
  func->SetLineWidth(1);
  func->SetLineStyle(2);

  TCanvas *c1 = new TCanvas();
  mg1->Draw("AP");
  //mg2->Draw("Psame");
  //gr1C->Draw("AP");
  gr1A->Draw("Psame");
  //gr2C->Draw("Psame");
  gr2A->Draw("Psame");
  gr3C->Draw("Psame");
  gr3A->Draw("Psame");
  l->Draw("same");
  func->Draw("same");

  return;
}
