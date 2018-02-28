#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH1F.h"
#include <TStyle.h>
#include <TGaxis.h>
#include <iostream>
#include "TPaveLabel.h"
#include <sstream>   //for istringstream

void ReverseXAxis (TH1 *h);

using namespace std;

#define NCHAN 150

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

void trky_hit_map(string hodo_tree, string dir, string phases, int maxfiles,
		  double e_thr, string comment) {

  //  string hTitle=Form("Total hit rate [Hz] (E >%4.1f MeV, ",e_thr) +
  //    comment + ")";

  TH1D* hTotHitPos = new TH1D("hh","top TrackerY", NCHAN+1,-0.5,NCHAN+0.5);
  hTotHitPos->GetXaxis()->SetTitle("channel");
  hTotHitPos->GetYaxis()->SetTitle("Rate [Hz]");
  TH1D *hTotHitNeg = (TH1D*)hTotHitPos->Clone("hTotHitNeg");
  hTotHitNeg->SetTitle("bottom TrackerY");

  Int_t nY = hTotHitPos->GetNbinsY();
  
  //Chain root files

  TChain ch(hodo_tree.c_str());   //a chain to process Tree hodo_tree
  istringstream iss;
  iss.str(phases);
  string phase;
  while (iss >> phase) {
    cout << "phase: " << phase << endl;
    for (int ifile=0; ifile<maxfiles; ifile++) {
      string rootfile = dir + "/" + phase + "/" + Form("%d.root",ifile);
      if (FILE *file = fopen(rootfile.c_str(), "r")) {
        fclose(file);
	ch.Add(rootfile.c_str(),0);
	cout << "  Add root file: " << rootfile << endl;
      }
      else
	cout << "  File: " << rootfile << " does not exist" << endl;
    }
  }

  std::vector<int> *det_vec = 0;
  std::vector<unsigned int> *chan_vec = 0;
  std::vector<double> *edep_vec = 0;

  TBranch *det_br = 0;
  TBranch *chan_br = 0;
  TBranch *edep_br = 0;

  ch.SetBranchAddress("detcont",&det_vec,&det_br);
  ch.SetBranchAddress("chancont",&chan_vec,&chan_br);
  ch.SetBranchAddress("edepcont",&edep_vec,&edep_br);

  Long64_t nentries = ch.GetEntriesFast();
  cout << "nentries = " << nentries << endl;

  for (Int_t i = 0; i < nentries; i++) {
      
    Long64_t tentry = ch.LoadTree(i);
    if (i < 10) cout << "tentry = " << tentry << "  i = " << i << endl;

    det_br->GetEntry(tentry);
    chan_br->GetEntry(tentry);
    edep_br->GetEntry(tentry);

    for (UInt_t j = 0; j < det_vec->size(); ++j) {
 
      int det = det_vec->at(j);
      unsigned int chan = chan_vec->at(j);
      double edep = edep_vec->at(j);
      //cout << " hit " << j << ": " << det << " " << chan << " "
      //	     << edep << endl;

      if (edep > e_thr) {
	if (det > 0)
	  //	  hTotHitPos->Fill(nY-1-chan-1,1.);
	  hTotHitPos->Fill(chan+1,1.);
	else
	  hTotHitNeg->Fill(chan+1,1.);   //mirror
      }

    }

    //      if (det_vec->size() !=0) getchar();

  }

  // Since we passed the address of a local variable we need
  // to remove it.
  ch.ResetBranchAddresses();

  // Scale histograms.

  const double CPSflux = 1.5e+12;      //photon/s

  hTotHitPos->Scale(CPSflux/double(nentries));
  hTotHitNeg->Scale(CPSflux/double(nentries));

   // Plot.

   gStyle->SetOptStat(0);

   TCanvas *ct = new TCanvas(("hit_map_" + hodo_tree).c_str(),
			     "hit map", 600, 720);

   TPaveLabel* title = new TPaveLabel(0.1,0.96,0.9,0.99,comment.c_str());
   title->SetBorderSize(0);
   title->SetFillColor(0);
   title->Draw();

   TPad* hPad = new TPad("histos","histos",0.01,0.05,0.95,0.95);
   hPad->Draw();
   hPad->cd();

   hPad->Divide(1,2);
   hPad->cd(1);
   hTotHitPos->Draw("hist");
   //   ct->SetRightMargin(0.15);
   //   ReverseXAxis(hTotHitPos);
   hPad->cd(2);
   hTotHitNeg->Draw("hist");
   //   ct->SetRightMargin(0.15);
   //   ReverseXAxis(hTotHitNeg);

}

//------------------------------------------------------------------------------

void ReverseXAxis (TH1 *h)
{
   // Remove the current axis
   h->GetXaxis()->SetLabelOffset(999);
   h->GetXaxis()->SetTickLength(0);

   // Redraw the new axis 
   gPad->Update();
   TGaxis *newaxis = new TGaxis(gPad->GetUxmax(), 
                                gPad->GetUymin(),
                                gPad->GetUxmin(),
                                gPad->GetUymin(),
                                h->GetXaxis()->GetXmin(),
                                h->GetXaxis()->GetXmax(),
                                510,"-");
   newaxis->SetLabelOffset(-0.03);
   newaxis->SetLabelFont(h->GetXaxis()->GetLabelFont());
   newaxis->SetLabelSize(h->GetXaxis()->GetLabelSize());
   newaxis->Draw();
}

void ReverseYAxis (TH1 *h)
{
   // Remove the current axis
   h->GetYaxis()->SetLabelOffset(999);
   h->GetYaxis()->SetTickLength(0);

   // Redraw the new axis 
   gPad->Update();
   TGaxis *newaxis = new TGaxis(gPad->GetUxmin(), 
                                gPad->GetUymax(),
                                gPad->GetUxmin()-0.001,
                                gPad->GetUymin(),
                                h->GetYaxis()->GetXmin(),
                                h->GetYaxis()->GetXmax(),
                                510,"+");
   newaxis->SetLabelOffset(-0.03);
   newaxis->Draw();
}
