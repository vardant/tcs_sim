#define deepgen_cxx
#include "deepgen.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// Create ascii file of TCS events in HEP format.

#define Mele 0.000511
#define Mpos 0.000511
#define Mpro 0.938272

#define Nele   11
#define Npos  -11
#define Npro 2212

#define PSF 269.5   //may be wrong.

void deepgen::Loop()
{
//   In a ROOT session, you can do:
//      root> .L deepgen.C
//      root> deepgen t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  if (fChain == 0) return;

  std::ofstream ofs("tcs_gen.data");
  std::ofstream ofs_kin("tcs_gen.kin_data");

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries = " << nentries << endl;
   
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
     
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    ofs << "3" << endl;
    ofs << "1 " << Nele << " 0 0 "
	<< ALV_minus_lab[1] << " "
	<< ALV_minus_lab[2] << " "
	<< ALV_minus_lab[3] << " "
	<< Mele << endl;
      ofs << "1 " << Npos << " 0 0 "
          << ALV_plus_lab[1] << " "
          << ALV_plus_lab[2] << " "
          << ALV_plus_lab[3] << " "
          << Mpos << endl;
      ofs << "1 " << Npro << " 0 0 "
          << ALV_Recoil_lab[1] << " "
          << ALV_Recoil_lab[2] << " "
          << ALV_Recoil_lab[3] << " "
          << Mpro << endl;

		     double s = 0.;
		     double tau = 0.;
		     double eta = 0.;
		     double xi = 0.;
		     if (Egamma > 0.) {
		       s = Mpro*Mpro + 2*Mpro*Egamma;
		       tau = Qp2/(s-Mpro*Mpro);
		       eta = tau/(tau - 2.);
		       xi = Qp2/(2.*(s-Mpro*Mpro)+tt-Qp2);
		     };
		     
      ofs_kin << Qp2 << " " << tt << " " << s << " " << xi << " " << tau << " "
	      << eta << " " << Phi_CMV << " " << Theta_CMV << " " << PSF << " "
	      << Flux_bmr << " " << cross_tot_unpol << " " << Egamma << endl;

      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}
