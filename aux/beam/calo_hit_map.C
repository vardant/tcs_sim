#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH2F.h"
#include <TStyle.h>
#include <TGaxis.h>
#include <iostream>
#include "TPaveLabel.h"
#include <sstream>   //for istringstream

void ReverseXAxis (TH1 *h);

using namespace std;

#define NCOL 50
#define NROW 23

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

void calo_hit_map(string dir, string phases, int maxfiles, double e_thr,
	     string comment) {

  TH2D* hTotHitPos = new TH2D("hh","Top Calorimeter",
			      NCOL+1,-0.5,NCOL+0.5,NROW+1,-0.5,NROW+0.5);
  hTotHitPos->GetXaxis()->SetTitle("column #");
  hTotHitPos->GetYaxis()->SetTitle("row #");
  TH2D *hTotHitNeg = (TH2D*)hTotHitPos->Clone("hTotHitNeg");
  hTotHitNeg->SetTitle("Bottom Calorimeter");

  Int_t nX = hTotHitPos->GetNbinsX();
  
  //Chain root files

  TChain ch("calo");   //a chain to process Tree "calo"
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
  std::vector<unsigned int> *col_vec = 0;
  std::vector<unsigned int> *row_vec = 0;
  std::vector<double> *edep_vec = 0;

  TBranch *det_br = 0;
  TBranch *col_br = 0;
  TBranch *row_br = 0;
  TBranch *edep_br = 0;

  ch.SetBranchAddress("detcont",&det_vec,&det_br);
  ch.SetBranchAddress("colcont",&col_vec,&col_br);
  ch.SetBranchAddress("rowcont",&row_vec,&row_br);
  ch.SetBranchAddress("edepcont",&edep_vec,&edep_br);

  Long64_t nentries = ch.GetEntriesFast();
  cout << "nentries = " << nentries << endl;

  for (Int_t i = 0; i < nentries; i++) {
      
    Long64_t tentry = ch.LoadTree(i);
    if (i < 10) cout << "tentry = " << tentry << "  i = " << i << endl;

    det_br->GetEntry(tentry);
    col_br->GetEntry(tentry);
    row_br->GetEntry(tentry);
    edep_br->GetEntry(tentry);

    for (UInt_t j = 0; j < det_vec->size(); ++j) {
 
      int det = det_vec->at(j);
      unsigned int col = col_vec->at(j);
      unsigned int row = row_vec->at(j);
      double edep = edep_vec->at(j);
      //cout << " hit " << j << ": " << det << " " << col << " " << row << " "
      //	     << edep << endl;

      if (edep > e_thr) {
	if (det > 0)
	  ////	  hTotHitPos->Fill(nX-1-col-1,row+1,1.);
	  	  hTotHitPos->Fill(col+1,row+1,1.);
	else
	  ////	  hTotHitNeg->Fill(nX-1-col-1,row+1,1.);   //mirror
	  	  hTotHitNeg->Fill(col+1,row+1,1.);   //mirror
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

   TCanvas *ct = new TCanvas("calo_hit_map", "Total hit rate [Hz/uA]",
			     600, 720);

   TPaveLabel* title = new TPaveLabel(0.1,0.96,0.9,0.99,comment.c_str());
   title->SetBorderSize(0);
   title->SetFillColor(0);
   title->Draw();

   TPad* hPad = new TPad("histos","histos",0.01,0.05,0.95,0.95);
   hPad->Draw();
   hPad->cd();

   hPad->Divide(1,2);
   hPad->cd(1);
   hTotHitPos->Draw("colz");
   ////   ct->SetRightMargin(0.15);
   ////   ReverseXAxis(hTotHitPos);
   hPad->cd(2);
   hTotHitNeg->Draw("colz");
   ////   ct->SetRightMargin(0.15);
   ////   ReverseXAxis(hTotHitNeg);

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
