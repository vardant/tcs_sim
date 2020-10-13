#define tracker_cxx
#define calo_cxx
#define kin_cxx
#define beam_cxx

#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TPad.h>
#include <iostream>

#include "include/tracker.h"
#include "include/calo.h"
#include "include/kin.h"
#include "include/beam.h"

using namespace std;

#define ELEID   11
#define POSID  -11
#define PROID 2212

#define NQUARTER 4
#define NLAYER   3

//#define ESUM_MIN 5.e3 //MeV, threshold on sum of deposited energies in calo-s.
#define ESUM_MIN 0.

//Plot deposited in the calorimeter energies from the TCS e+, e- tracks
//for the trigger concept development.

//This is an adaptation of aceptance_tag2 code.

int n_track_points(int quarter, const bool track[NQUARTER][NLAYER]);

void tcs_edep(string tracked_file) {

  tracker Tracker(tracked_file.c_str());
  calo Calo(tracked_file.c_str());
  kin Kin(tracked_file.c_str());
  beam Beam(tracked_file.c_str());

  TH1D* hedep_ele = new TH1D("hedep_ele", "Edep(e-)", 110, 0., 11000.);
  TH1D* hedep_pos = new TH1D("hedep_pos", "Edep(e+)", 110, 0., 11000.);
  TH1D* hedep_sum = new TH1D("hedep_sum", "Edep(e+ + e-)", 110, 0., 11000.);
  TH2D* hedep2= new TH2D("hedep", "Edep e+ vs e-", 110, 0., 11000.,
			                           110, 0., 11000.);

  TH1D* he_ele_vrtx = new TH1D("he_ele_vrtx", "E(e-) at vertex", 110,0.,11000.);

  TH2D* hedep_evrtx = new TH2D("hedep_evrtx", "Edep vs Evertex for e-",
			       110, 0., 11000., 110, 0., 11000.);
  hedep_evrtx->GetXaxis()->SetTitle("E_{Vertex} [MeV]");
  hedep_evrtx->GetYaxis()->SetTitle("E_{DEP} [MeV]");

  TH1D* hedep_norm = new TH1D("hedep_norm", "E_{DEP}/E_{VERTEX} (e-)",
			      150,0.,1.5);

  //Loop over events.

  Long64_t nentries = Calo.GetEntriesFast();
  cout << "calo nentries = " << nentries << endl;

  int n_accepted = 0;
  int n_opposite = 0;
  int n_same_quad = 0;
  int n_neighbour = 0;

  //  for (int ientry=0; ientry<100000; ientry++) {
  for (int ientry=0; ientry<nentries; ientry++) {

    bool ele_track[NQUARTER][NLAYER];
    bool pos_track[NQUARTER][NLAYER];
    bool pro_track[NQUARTER][NLAYER];
    bool ele_calo[NQUARTER];
    bool pos_calo[NQUARTER];

    //Check tracks passed through GEM trackers and deposited energies
    //in the calorimeters.

    for (int q=0; q<NQUARTER; q++) {
      for (int l=0; l<NLAYER; l++) {
	ele_track[q][l] = Tracker.CheckTrack(ientry, ELEID, l, q);
	pos_track[q][l] = Tracker.CheckTrack(ientry, POSID, l, q);
	pro_track[q][l] = Tracker.CheckTrack(ientry, PROID, l, q);
      }
      ele_calo[q] = Calo.CheckEdep(ientry, ELEID, q) &&
	            Calo.CheckAccept(ientry, ELEID, q);
      pos_calo[q] = Calo.CheckEdep(ientry, POSID, q) &&
	            Calo.CheckAccept(ientry, POSID, q);
    }

    /*
    int n_fired = 0;
    for (int q=0; q<NQUARTER; q++) {
      for (int l=0; l<NLAYER; l++) {
	if (pos_track[q][l]) n_fired++;
      }
    }

    if (n_fired>0) {
      cout << "e+ tracks:" << endl;
      for (int q=0; q<NQUARTER; q++) {
	cout << " quarter " << q << ":";
	for (int l=0; l<NLAYER; l++) {
	  cout << " " << pos_track[q][l];
	}
	cout << "\n";
	cout << " calo     : " << pos_calo[q] << endl;
      }
      getchar();
    }
    */
    /*
    int n_fired = 0;
    for (int q=0; q<NQUARTER; q++) {
      for (int l=0; l<NLAYER; l++) {
	if (pro_track[q][l]) n_fired++;
      }
    }

    if (n_fired>0) {
      cout << "p tracks:" << endl;
      for (int q=0; q<NQUARTER; q++) {
	cout << " quarter " << q << ":";
	for (int l=0; l<NLAYER; l++) {
	  cout << " " << pro_track[q][l];
	}
	cout << "\n";
///	cout << " hodox    : " << pro_hodox[q] << endl;
///	cout << " hodoy    : " << pro_hodoy[q] << endl;
      }
      getchar();
    }
    */

    //Check if enough hits from GEM trackers for track reconstruction,
    //and signals from calorimeters.

    int q_ele = -1;
    double edep_ele = 0.;
    bool good_ele = false;
    for (int q=0; q<NQUARTER; q++) {
      if (n_track_points(q, ele_track) > 1 && ele_calo[q]) {
	edep_ele = Calo.Edep(ientry, ELEID, q);
	q_ele = q;
	good_ele = true;
	break;
      }
    }

    int q_pos = -1;
    double edep_pos = 0.;
    bool good_pos = false;
    for (int q=0; q<NQUARTER; q++) {
      if (n_track_points(q, pos_track) > 1 && pos_calo[q]) {
	edep_pos = Calo.Edep(ientry, POSID, q);
	q_pos = q;
	good_pos = true;
	break;
      }
    }

    bool good_pro = false;
    for (int q=0; q<NQUARTER; q++) {
      ///if (n_track_points(q, pro_track) > 1 && pro_hodox[q] && pro_hodoy[q]) {
      if (n_track_points(q, pro_track) > 1) {
	good_pro = true;
	break;
      }
    }

    //Request triple coincidence.
    bool good_esum = edep_ele + edep_pos > ESUM_MIN;
    bool acceptance_flag = good_ele && good_pos && good_pro && good_esum &&
      Kin.Cut(ientry);

    if (acceptance_flag) {

      if (q_ele<0 || q_pos<0)
	cout << "*** Wrong quadrant: q_ele = " << q_ele << "  q_pos = " << q_pos
	     << " ***" << endl;

      if (q_ele == q_pos)
	n_same_quad++;
      else if (q_ele+q_pos == 2 || q_ele+q_pos == 4)
	n_opposite++;
      else
	n_neighbour++;

      //      double w = Kin.GetWeight(ientry)/nentries; No x-sec for now
      double w=1.;
      hedep_ele->Fill(edep_ele,w);
      hedep_pos->Fill(edep_pos,w);
      hedep_sum->Fill(edep_pos+edep_ele,w);
      hedep2->Fill(edep_ele,edep_pos,w);
      //      cout << edep_ele << endl;

      Beam.GetEntry(ientry);
      he_ele_vrtx->Fill(Beam.e,w);      //e- energy
      hedep_evrtx->Fill(Beam.e,edep_ele,w);
      hedep_norm->Fill(edep_ele/Beam.e,w);

      n_accepted++;
    }

    //cout << good_ele << good_pos << good_pro << " " << acceptance_flag <<endl;
    //getchar();

  }   //entries

  cout << "n_same_quad = " << n_same_quad << "  n_neighbour = " << n_neighbour
       << "  n_opposite = " << n_opposite << "  n_accepted = " << n_accepted
       << endl;

  TCanvas* cedep_ele = new TCanvas();
  hedep_ele->Draw("HIST");
  cedep_ele->Print("cedep_ele.pdf");

  TCanvas* cedep_pos = new TCanvas();
  hedep_pos->Draw("HIST");
  cedep_pos->Print("cedep_pos.pdf");

  TCanvas* cedep_sum = new TCanvas();
  hedep_sum->Draw("HIST");
  cedep_sum->Print("cedep_sum.pdf");

  TCanvas* cedep2 = new TCanvas();
  hedep2->Draw();
  cedep2->Print("cedep2.pdf");

  TCanvas* ce_ele_vrtx = new TCanvas();
  he_ele_vrtx->Draw("HIST");
  ce_ele_vrtx->Print("ce_ele_vrtx.pdf");

  TCanvas* cedep_evrtx = new TCanvas();
  hedep_evrtx->Draw();
  cedep_evrtx->Print("cedep_evrtx.pdf");

  TCanvas* cedep_norm = new TCanvas();
  gPad->SetLogy();
  hedep_norm->Draw();
  cedep_norm->Print("cedep_norm.pdf");

}

//------------------------------------------------------------------------------

int n_track_points(int quarter, const bool track[NQUARTER][NLAYER]) {
  int n = 0;
  for (int l=0; l<NLAYER; l++) {
    if (track[quarter][l])
      n++;
  }
  return n;
}
