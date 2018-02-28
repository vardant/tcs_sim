#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH2.h"
#include <TStyle.h>
#include <TGaxis.h>
#include <iostream>
#include "TPaveLabel.h"
#include <sstream>   //for istringstream

using namespace std;

#define NXCHAN 332
#define NYCHAN 150

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

void trkxy_hit_map(string dir, string phases, int maxfiles,
		  double e_thr, string comment) {

  TH2D* hTotHitPos = new TH2D("hh","Top TrackerY vs TrackerX",
			     NXCHAN+1,-0.5,NXCHAN+0.5,NYCHAN+1,-0.5,NYCHAN+0.5);
  hTotHitPos->GetXaxis()->SetTitle("TrackerX chanel");
  hTotHitPos->GetYaxis()->SetTitle("TrackerY chanel");
  TH2D *hTotHitNeg = (TH2D*)hTotHitPos->Clone("hTotHitNeg");
  hTotHitNeg->SetTitle("bottom TrackerY vs TrackerX");

  //Chain root files

  TChain chx("trackerx");
  TChain chy("trackery");
  istringstream iss;
  iss.str(phases);
  string phase;
  while (iss >> phase) {
    cout << "phase: " << phase << endl;
    for (int ifile=0; ifile<maxfiles; ifile++) {
      string rootfile = dir + "/" + phase + "/" + Form("%d.root",ifile);
      if (FILE *file = fopen(rootfile.c_str(), "r")) {
        fclose(file);
	chx.Add(rootfile.c_str(),0);
	chy.Add(rootfile.c_str(),0);
	cout << "  Add root file: " << rootfile << endl;
      }
      else
	cout << "  File: " << rootfile << " does not exist" << endl;
    }
  }

  // X.
  
  std::vector<int> *detx_vec = 0;
  std::vector<unsigned int> *chanx_vec = 0;
  std::vector<double> *edepx_vec = 0;

  TBranch *detx_br = 0;
  TBranch *chanx_br = 0;
  TBranch *edepx_br = 0;

  chx.SetBranchAddress("detcont",&detx_vec,&detx_br);
  chx.SetBranchAddress("chancont",&chanx_vec,&chanx_br);
  chx.SetBranchAddress("edepcont",&edepx_vec,&edepx_br);

  // Y.
  
  std::vector<int> *dety_vec = 0;
  std::vector<unsigned int> *chany_vec = 0;
  std::vector<double> *edepy_vec = 0;

  TBranch *dety_br = 0;
  TBranch *chany_br = 0;
  TBranch *edepy_br = 0;

  chy.SetBranchAddress("detcont",&dety_vec,&dety_br);
  chy.SetBranchAddress("chancont",&chany_vec,&chany_br);
  chy.SetBranchAddress("edepcont",&edepy_vec,&edepy_br);

  Long64_t nentries = chx.GetEntriesFast();
  cout << "nentries = " << nentries << endl;

  for (Int_t i = 0; i < nentries; i++) {
      
    Long64_t tentryx = chx.LoadTree(i);
    if (i < 10) cout << "tentryx = " << tentryx << "  i = " << i << endl;
    detx_br->GetEntry(tentryx);
    chanx_br->GetEntry(tentryx);
    edepx_br->GetEntry(tentryx);

    Long64_t tentryy = chy.LoadTree(i);
    dety_br->GetEntry(tentryy);
    chany_br->GetEntry(tentryy);
    edepy_br->GetEntry(tentryy);

    for (UInt_t j = 0; j < detx_vec->size(); ++j) {
 
      int detx = detx_vec->at(j);
      unsigned int chanx = chanx_vec->at(j);
      double edepx = edepx_vec->at(j);

      for (UInt_t k = 0; k < dety_vec->size(); ++k) {
	
	int dety = dety_vec->at(k);
	unsigned int chany = chany_vec->at(k);
	double edepy = edepy_vec->at(k);
      
	//cout << " hit " << j << ": " << det << " " << chan << " "
	//	     << edep << endl;

	if (edepx > e_thr && edepy > e_thr) {
	  if (detx > 0 && dety > 0)
	    hTotHitPos->Fill(chanx+1,chany+1,1.);
	  else
	    if (detx < 0 && dety < 0)
	      hTotHitNeg->Fill(chanx+1,chany+1,1.);   //mirror
	}
	
      }

    }

  }

    //      if (det_vec->size() !=0) getchar();

  // Since we passed the address of a local variable we need
  // to remove it.
  ////  ch.ResetBranchAddresses();

  // Scale histograms.

  const double CPSflux = 1.5e+12;      //photon/s

  hTotHitPos->Scale(CPSflux/double(nentries));
  hTotHitNeg->Scale(CPSflux/double(nentries));

   // Plot.

   gStyle->SetOptStat(0);

   TCanvas *ct = new TCanvas("trkxy_hit_map", comment.c_str(), 600, 720);

   TPaveLabel* title = new TPaveLabel(0.1,0.96,0.9,0.99,comment.c_str());
   title->SetBorderSize(0);
   title->SetFillColor(0);
   title->Draw();

   TPad* hPad = new TPad("histos","histos",0.01,0.05,0.95,0.95);
   hPad->Draw();
   hPad->cd();

   hPad->Divide(1,2);
   hPad->cd(1);
   hTotHitPos->Draw("cont");
   hPad->cd(2);
   hTotHitNeg->Draw("cont");

}
