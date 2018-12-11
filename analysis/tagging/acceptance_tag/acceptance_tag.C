#define sim_tcs_cxx
#define tracker_cxx
#define hodox_cxx
#define hodoy_cxx
#define calo_cxx

#include "include/sim_tcs.h"
#include <iostream>
#include "include/tracker.h"
#include "include/hodox.h"
#include "include/hodoy.h"
#include "include/calo.h"

using namespace std;

#define ELEID   11
#define POSID  -11
#define PROID 2212

#define NQUARTER 4
#define NLAYER   3

//Acceptance tag TCS events. Tracks are accepted if they pass through
//at least 2 GEM trackers. e-. e+ are accepted if there is signifant
//energy deposition in the calorimeters. Recoil proton is accepted if
//there are significant energy depositions in the X, Y hodoscopes.

int n_track_points(int quarter, const bool track[NQUARTER][NLAYER]);

void acceptance_tag(string prefix) {

  //TCS simulation input file and class object.
  string primary_file = "RootFiles/"+prefix+".root";
  sim_tcs SimTCS(primary_file.c_str());

  //Tracked throgh setup events.
  string tracked_file = "RootFiles.tracked/"+prefix+"_tracked.root";
  tracker Tracker(tracked_file.c_str());
  calo Calo(tracked_file.c_str());
  hodox Hodox(tracked_file.c_str());
  hodoy Hodoy(tracked_file.c_str());

  //Output root file.
  string tagged_file  = "RootFiles.tagged/"+prefix+"_tagged.root";
  TFile* tag_file = new TFile(tagged_file.c_str(),"recreate");

  //Add branch of tag flag to the tree of simulated events.

  TTree *primTree = SimTCS.GetTree();
  TTree* taggedTree = primTree->CloneTree();

  int acceptance_flag = 0;
  TBranch *acceptance_br = taggedTree->Branch("acceptance_flag",
                                              &acceptance_flag,
                                              "acceptance_flag/I");

  //Loop over events.

  Long64_t nentries = SimTCS.GetEntriesFast();
  cout << "SimTCS nentries = " << nentries << endl;

  for (int ientry=0; ientry<nentries; ientry++) {

    bool ele_track[NQUARTER][NLAYER];
    bool pos_track[NQUARTER][NLAYER];
    bool pro_track[NQUARTER][NLAYER];
    bool ele_calo[NQUARTER];
    bool pos_calo[NQUARTER];
    bool pro_hodox[NQUARTER];
    bool pro_hodoy[NQUARTER];

    //Check tracks passed through GEM trackers and deposited energies
    //in the calorimeters and hodoscopes.

    for (int q=0; q<NQUARTER; q++) {
      for (int l=0; l<NLAYER; l++) {
	ele_track[q][l] = Tracker.CheckTrack(ientry, ELEID, l, q);
	pos_track[q][l] = Tracker.CheckTrack(ientry, POSID, l, q);
	pro_track[q][l] = Tracker.CheckTrack(ientry, PROID, l, q);
      }
      ele_calo[q] = Calo.CheckEdep(ientry, ELEID, q);
      pos_calo[q] = Calo.CheckEdep(ientry, POSID, q);
      pro_hodox[q] = Hodox.CheckEdep(ientry, PROID, q);
      pro_hodoy[q] = Hodoy.CheckEdep(ientry, PROID, q);
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
	cout << " hodox    : " << pro_hodox[q] << endl;
	cout << " hodoy    : " << pro_hodoy[q] << endl;
      }
      getchar();
    }
    */

    //Check if enough hits from GEM trackers for track reconstruction,
    //and signals from from calorimeters and hodoscopes.

    bool good_ele = false;
    for (int q=0; q<NQUARTER; q++) {
      if (n_track_points(q, ele_track) > 1 && ele_calo[q]) {
	good_ele = true;
	break;
      }
    }

    bool good_pos = false;
    for (int q=0; q<NQUARTER; q++) {
      if (n_track_points(q, pos_track) > 1 && pos_calo[q]) {
	good_pos = true;
	break;
      }
    }

    bool good_pro = false;
    for (int q=0; q<NQUARTER; q++) {
      if (n_track_points(q, pro_track) > 1 && pro_hodox[q] && pro_hodoy[q]) {
	good_pro = true;
	break;
      }
    }

    //Request triple coincidence.
    acceptance_flag = good_ele && good_pos && good_pro;
    acceptance_br->Fill();

    //cout << good_ele << good_pos << good_pro << " " << acceptance_flag <<endl;
    //getchar();

  }   //entries

  tag_file->cd();
  taggedTree->Write();

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
