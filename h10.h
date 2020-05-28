//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb 26 11:22:22 2009 by ROOT version 5.18/00
// from TTree h10/All_out
// found on file: root_42011_01.pass2.root
//////////////////////////////////////////////////////////

#ifndef h10_h
#define h10_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TVirtualIndex.h>
#include <TFile.h>
#include <TCanvas.h>
#include <iostream>             // std::cout, std::endl
#include <fstream>              // std::ifstream
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TH2F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TPostScript.h>
#include <TPaveText.h>

using namespace std;

 class h10 : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t AvailableEvents, Nevta, Nevt_found;
   TFile *FileGam;
   TTree *TreeGam;
  
   Float_t eND, eNS;
   Float_t Npi0D, Npi0S;
   Float_t Rh;
   Int_t trig_dcfid_flag,target;
   Int_t trig_sect, sect;
   Float_t trig_phisect;
   Float_t trig_ql;
  
   TLorentzVector V4Beam, V4Elec, V4Q, V4Pi0, V4Pi0corr, VG1corr, VG2corr,  V4Prot;
   Int_t trig_sc_stat; 
   Float_t trig_mom, trig_theta,trig_theta_rad, trig_phi, trig_nphe, trig_ecin, trig_ecout, trig_ectot, trig_ece;
   Float_t trig_ecx, trig_ecy, trig_ecz, trig_ecu, trig_ecv, trig_ecw, trig_ect, trig_ecl;
   Float_t trig_sct, trig_scl, trig_dcECx, trig_dcECy, trig_dcECz;
   Float_t trig_cctheta, trig_ccphi;
   Float_t trig_W,trig_Q2, trig_y, trig_nu;   
   Float_t trig_vz, trig_vx, trig_vy;
   Float_t trig_vz_corr, trig_vx_corr, trig_vy_corr, trig_transv_pos; 
   Float_t trig_px,trig_py,trig_pz;
   Float_t trig_tl1x,  trig_tl1y, trig_tl1z,trig_tl1cx,  trig_tl1cy, trig_tl1cz;
   Float_t trig_phil1, trig_thetal1;
   Int_t   LorenzoFlag;
   Bool_t NpheCutSect;
   Int_t MirrorCode;
 
   Int_t nphot;
   Float_t gam_Eec[10], gam_bet[10], gam_E[10];
   Float_t  gam_Etrue_Corr[10],gam_Pxtrue_Corr[10], gam_Pytrue_Corr[10], gam_Pztrue_Corr[10];
   Float_t gam_Etrue[10], gam_Pxtrue[10], gam_Pytrue[10], gam_Pztrue[10];
   Float_t gam_ecx[10], gam_ecy[10], gam_ecz[10], gam_ecu[10], gam_ecv[10], gam_ecw[10], gam_ect[10], gam_ecpath[10];   
   Float_t gam_Px[10], gam_Py[10], gam_Pz[10];
   Float_t gam_theta[10], gam_phi[10];
   
   Int_t npi0, G1[45], G2[45];
   Float_t Pi0mass_EVNT[45], Pi0z_EVNT[45],Pi0energy_EVNT[45];
   Float_t opening_angle_evnt[45]; 
   Float_t Pi0mass_ECPB[45], Pi0z_ECPB[45],Pi0energy_ECPB[45];
   Float_t opening_angle_ecpb[45], opening_angle_Corr[45];
   Float_t pT2_val[45], pT2_val_Corr[45], phi_val_Corr[45]; 
   Float_t Pi0mass_Corr[45], Pi0z_Corr[45],Pi0energy_Corr[45];
   Float_t eG1Angle[45], eG2Angle[45];
   Float_t Pi0phi_Corr[45], Pi0theta_Corr[45], Pi0angle_gs[45];
   Float_t pionW[45];
  
  /* 
   Double_t PREDUMMY[200][200];
   Float_t Mx2[45], Mx2_fixed[45];  
   Double_t POSTDUMMY[200][200];
*/

  // Declaration of leaf types
   UChar_t         npart;
   UChar_t         evstat;
   UInt_t          evntid;
   Char_t          evntype;
   Int_t           evntclas;
   Float_t         q_l;
   Float_t         t_l;
   Float_t         tr_time;
   Float_t         rf_time;
   Int_t           l2bit;
   Int_t           l3bit;
   Int_t           hlsc;
   Int_t           intt;
   Int_t           gpart;
   Int_t           id[40];   //[gpart]
   Int_t           stat[40];   //[gpart]
   Int_t           dc[40];   //[gpart]
   Int_t           cc[40];   //[gpart]
   Int_t           sc[40];   //[gpart]
   Int_t           ec[40];   //[gpart]
   Int_t           lec[40];   //[gpart]
   Int_t           st[40];   //[gpart]
   Float_t         p[40];   //[gpart]
   Float_t         m[40];   //[gpart]
   Int_t           q[40];   //[gpart]
   Float_t         b[40];   //[gpart]
   Float_t         cx[40];   //[gpart]
   Float_t         cy[40];   //[gpart]
   Float_t         cz[40];   //[gpart]
   Float_t         vx[40];   //[gpart]
   Float_t         vy[40];   //[gpart]
   Float_t         vz[40];   //[gpart]
   Int_t           dc_part;
   Int_t           dc_sect[40];   //[dc_part]
   Int_t           dc_trk[40];   //[dc_part]
   Int_t           dc_stat[40];   //[dc_part]
   Int_t           tb_st[40];   //[dc_part]
   Float_t         dc_xsc[40];   //[dc_part]
   Float_t         dc_ysc[40];   //[dc_part]
   Float_t         dc_zsc[40];   //[dc_part]
   Float_t         dc_cxsc[40];   //[dc_part]
   Float_t         dc_cysc[40];   //[dc_part]
   Float_t         dc_czsc[40];   //[dc_part]
   Float_t         dc_vx[40];   //[dc_part]
   Float_t         dc_vy[40];   //[dc_part]
   Float_t         dc_vz[40];   //[dc_part]
   Float_t         dc_vr[40];   //[dc_part]
   Float_t         tl1_cx[40];   //[dc_part]
   Float_t         tl1_cy[40];   //[dc_part]
   Float_t         tl1_cz[40];   //[dc_part]
   Float_t         tl1_x[40];   //[dc_part]
   Float_t         tl1_y[40];   //[dc_part]
   Float_t         tl1_z[40];   //[dc_part]
   Float_t         tl1_r[40];   //[dc_part]
   Float_t         dc_c2[40];   //[dc_part]
   Int_t           ec_part;
   Int_t           ec_stat[40];   //[ec_part]
   Int_t           ec_sect[40];   //[ec_part]
   Int_t           ec_whol[40];   //[ec_part]
   Int_t           ec_inst[40];   //[ec_part]
   Int_t           ec_oust[40];   //[ec_part]
   Float_t         etot[40];   //[ec_part]
   Float_t         ec_ei[40];   //[ec_part]
   Float_t         ec_eo[40];   //[ec_part]
   Float_t         ec_t[40];   //[ec_part]
   Float_t         ec_r[40];   //[ec_part]
   Float_t         ech_x[40];   //[ec_part]
   Float_t         ech_y[40];   //[ec_part]
   Float_t         ech_z[40];   //[ec_part]
   Float_t         ec_m2[40];   //[ec_part]
   Float_t         ec_m3[40];   //[ec_part]
   Float_t         ec_m4[40];   //[ec_part]
   Float_t         ec_c2[40];   //[ec_part]
   Int_t           sc_part;
   Int_t           sc_sect[40];   //[sc_part]
   Int_t           sc_hit[40];   //[sc_part]
   Int_t           sc_pd[40];   //[sc_part]
   Int_t           sc_stat[40];   //[sc_part]
   Float_t         edep[40];   //[sc_part]
   Float_t         sc_t[40];   //[sc_part]
   Float_t         sc_r[40];   //[sc_part]
   Float_t         sc_c2[40];   //[sc_part]
   Int_t           cc_part;
   Int_t           cc_sect[40];   //[cc_part]
   Int_t           cc_hit[40];   //[cc_part]
   Int_t           cc_segm[40];   //[cc_part]
   Int_t           nphe[40];   //[cc_part]
   Float_t         cc_t[40];   //[cc_part]
   Float_t         cc_r[40];   //[cc_part]
   Float_t         cc_c2[40];   //[cc_part]
   Int_t           lac_part;
   Int_t           lec_sect[40];   //[lac_part]
   Int_t           lec_hit[40];   //[lac_part]
   Int_t           lec_stat[40];   //[lac_part]
   Float_t         lec_etot[40];   //[lac_part]
   Float_t         lec_ein[40];   //[lac_part]
   Float_t         lec_t[40];   //[lac_part]
   Float_t         lec_r[40];   //[lac_part]
   Float_t         lec_x[40];   //[lac_part]
   Float_t         lec_y[40];   //[lac_part]
   Float_t         lec_z[40];   //[lac_part]
   Float_t         lec_c2[40];   //[lac_part]
   Int_t           vidmvrt;
   Int_t           ntrmvrt;
   Float_t         xmvrt;
   Float_t         ymvrt;
   Float_t         zmvrt;
   Float_t         ch2mvrt;
   Float_t         cxxmvrt;
   Float_t         cxymvrt;
   Float_t         cxzmvrt;
   Float_t         cyymvrt;
   Float_t         cyzmvrt;
   Int_t           stamvrt;
   Int_t           mcnentr;
   UChar_t         mcnpart;
   Int_t           mcst[20];   //[mcnentr]
   Int_t           mcid[20];   //[mcnentr]
   Int_t           mcpid[20];   //[mcnentr]
   Float_t         mctheta[20];   //[mcnentr]
   Float_t         mcphi[20];   //[mcnentr]
   Float_t         mcp[20];   //[mcnentr]
   Float_t         mcm[20];   //[mcnentr]
   Float_t         mcvx[20];   //[mcnentr]
   Float_t         mcvy[20];   //[mcnentr]
   Float_t         mcvz[20];   //[mcnentr]
   Float_t         mctof[20];   //[mcnentr]
   Int_t           nprt;
   Int_t           pidpart[20];   //[nprt]
   Float_t         xpart[20];   //[nprt]
   Float_t         ypart[20];   //[nprt]
   Float_t         zpart[20];   //[nprt]
   Float_t         epart[20];   //[nprt]
   Float_t         pxpart[20];   //[nprt]
   Float_t         pypart[20];   //[nprt]
   Float_t         pzpart[20];   //[nprt]
   Float_t         qpart[20];   //[nprt]
   Int_t           Ipart10[20];   //[nprt]
   Float_t         Rpart11[20];   //[nprt]
   Float_t         Rpart12[20];   //[nprt]
   Int_t           Ipart13[20];   //[nprt]
   Int_t           nschit;
   Int_t           scsect[10];   //[nschit]
   Int_t           schid[10];   //[nschit]
   Int_t           scpid[10];   //[nschit]
   Float_t         sct[10];   //[nschit]
   Float_t         sce[10];   //[nschit]
   Float_t         scx[10];   //[nschit]
   Float_t         scy[10];   //[nschit]
   Float_t         scz[10];   //[nschit]

   // List of branches
   TBranch        *b_npart;   //!
   TBranch        *b_evstat;   //!
   TBranch        *b_evntid;   //!
   TBranch        *b_evntype;   //!
   TBranch        *b_evntclas;   //!
   TBranch        *b_q_l;   //!
   TBranch        *b_t_l;   //!
   TBranch        *b_tr_time;   //!
   TBranch        *b_rf_time;   //!
   TBranch        *b_l2bit;   //!
   TBranch        *b_l3bit;   //!
   TBranch        *b_hlsc;   //!
   TBranch        *b_intt;   //!
   TBranch        *b_gpart;   //!
   TBranch        *b_id;   //!
   TBranch        *b_stat;   //!
   TBranch        *b_dc;   //!
   TBranch        *b_cc;   //!
   TBranch        *b_sc;   //!
   TBranch        *b_ec;   //!
   TBranch        *b_lec;   //!
   TBranch        *b_st;   //!
   TBranch        *b_p;   //!
   TBranch        *b_m;   //!
   TBranch        *b_q;   //!
   TBranch        *b_b;   //!
   TBranch        *b_cx;   //!
   TBranch        *b_cy;   //!
   TBranch        *b_cz;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_dc_part;   //!
   TBranch        *b_dc_sect;   //!
   TBranch        *b_dc_trk;   //!
   TBranch        *b_dc_stat;   //!
   TBranch        *b_tb_st;   //!
   TBranch        *b_dc_xsc;   //!
   TBranch        *b_dc_ysc;   //!
   TBranch        *b_dc_zsc;   //!
   TBranch        *b_dc_cxsc;   //!
   TBranch        *b_dc_cysc;   //!
   TBranch        *b_dc_czsc;   //!
   TBranch        *b_dc_vx;   //!
   TBranch        *b_dc_vy;   //!
   TBranch        *b_dc_vz;   //!
   TBranch        *b_dc_vr;   //!
   TBranch        *b_tl1_cx;   //!
   TBranch        *b_tl1_cy;   //!
   TBranch        *b_tl1_cz;   //!
   TBranch        *b_tl1_x;   //!
   TBranch        *b_tl1_y;   //!
   TBranch        *b_tl1_z;   //!
   TBranch        *b_tl1_r;   //!
   TBranch        *b_dc_c2;   //!
   TBranch        *b_ec_part;   //!
   TBranch        *b_ec_stat;   //!
   TBranch        *b_ec_sect;   //!
   TBranch        *b_ec_whol;   //!
   TBranch        *b_ec_inst;   //!
   TBranch        *b_ec_oust;   //!
   TBranch        *b_etot;   //!
   TBranch        *b_ec_ei;   //!
   TBranch        *b_ec_eo;   //!
   TBranch        *b_ec_t;   //!
   TBranch        *b_ec_r;   //!
   TBranch        *b_ech_x;   //!
   TBranch        *b_ech_y;   //!
   TBranch        *b_ech_z;   //!
   TBranch        *b_ec_m2;   //!
   TBranch        *b_ec_m3;   //!
   TBranch        *b_ec_m4;   //!
   TBranch        *b_ec_c2;   //!
   TBranch        *b_sc_part;   //!
   TBranch        *b_sc_sect;   //!
   TBranch        *b_sc_hit;   //!
   TBranch        *b_sc_pd;   //!
   TBranch        *b_sc_stat;   //!
   TBranch        *b_edep;   //!
   TBranch        *b_sc_t;   //!
   TBranch        *b_sc_r;   //!
   TBranch        *b_sc_c2;   //!
   TBranch        *b_cc_part;   //!
   TBranch        *b_cc_sect;   //!
   TBranch        *b_cc_hit;   //!
   TBranch        *b_cc_segm;   //!
   TBranch        *b_nphe;   //!
   TBranch        *b_cc_t;   //!
   TBranch        *b_cc_r;   //!
   TBranch        *b_cc_c2;   //!
   TBranch        *b_lac_part;   //!
   TBranch        *b_lec_sect;   //!
   TBranch        *b_lec_hit;   //!
   TBranch        *b_lec_stat;   //!
   TBranch        *b_lec_etot;   //!
   TBranch        *b_lec_ein;   //!
   TBranch        *b_lec_t;   //!
   TBranch        *b_lec_r;   //!
   TBranch        *b_lec_x;   //!
   TBranch        *b_lec_y;   //!
   TBranch        *b_lec_z;   //!
   TBranch        *b_lec_c2;   //!
   TBranch        *b_vidmvrt;   //!
   TBranch        *b_ntrmvrt;   //!
   TBranch        *b_xmvrt;   //!
   TBranch        *b_ymvrt;   //!
   TBranch        *b_zmvrt;   //!
   TBranch        *b_ch2mvrt;   //!
   TBranch        *b_cxxmvrt;   //!
   TBranch        *b_cxymvrt;   //!
   TBranch        *b_cxzmvrt;   //!
   TBranch        *b_cyymvrt;   //!
   TBranch        *b_cyzmvrt;   //!
   TBranch        *b_stamvrt;   //!
   TBranch        *b_mcnentr;   //!
   TBranch        *b_mcnpart;   //!
   TBranch        *b_mcst;   //!
   TBranch        *b_mcid;   //!
   TBranch        *b_mcpid;   //!
   TBranch        *b_mctheta;   //!
   TBranch        *b_mcphi;   //!
   TBranch        *b_mcp;   //!
   TBranch        *b_mcm;   //!
   TBranch        *b_mcvx;   //!
   TBranch        *b_mcvy;   //!
   TBranch        *b_mcvz;   //!
   TBranch        *b_mctof;   //!
   TBranch        *b_nprt;   //!
   TBranch        *b_pidpart;   //!
   TBranch        *b_xpart;   //!
   TBranch        *b_ypart;   //!
   TBranch        *b_zpart;   //!
   TBranch        *b_epart;   //!
   TBranch        *b_pxpart;   //!
   TBranch        *b_pypart;   //!
   TBranch        *b_pzpart;   //!
   TBranch        *b_qpart;   //!
   TBranch        *b_Ipart10;   //!
   TBranch        *b_Rpart11;   //!
   TBranch        *b_Rpart12;   //!
   TBranch        *b_Ipart13;   //!
   TBranch        *b_nschit;   //!
   TBranch        *b_scsect;   //!
   TBranch        *b_schid;   //!
   TBranch        *b_scpid;   //!
   TBranch        *b_sct;   //!
   TBranch        *b_sce;   //!
   TBranch        *b_scx;   //!
   TBranch        *b_scy;   //!
   TBranch        *b_scz;   //!


   h10 (TTree * tree =0, int Ntoproc=0) {AvailableEvents=Ntoproc; }
   virtual ~h10() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   Bool_t GoodTrig(void);
   void MakeBranch(void);
   Bool_t IsAPhoton(int track);
   void MakePhoton(int part);
   // Bool_t IsAProton(int track);
  //  void MakeProton(int part);
   
   Bool_t SearchPi0s(void);
   void SortPhotons(void);
   void SwapPhotons(int p1, int p2);
   template <class T> void MySwap (T* a, T* b);

   TVector3 ec_xyz_uvw(TVector3 xyz);
   Bool_t SamplingCutC(void);
   Bool_t SamplingCutFe(void);
   Bool_t SamplingCutPb(void);
   Bool_t NpheSect(void);
   Bool_t TrigCCmatch(void); 
   
  TString fInFileLorenzo;
   TF1* fcfun[7][2][4];
   TF1* fcfidc_fun[2];
   Int_t CheckCut(void);
   void TFidCut(Char_t* FileName);
  
   void FillVlassovCoordinates(int track);
   void get_cc_special_coord(float *point, float *dir, float *cc_coor);
   int get_vcrpl(float *r0, float *dir, float *plane_par, float *dist, float *cross_point);

   float calculate_eG1_angle(void);       
   float calculate_eG2_angle(void);
   float calculate_pT2(void);
   float calculate_pT2_corr(void);
   float calculate_phi_corr(void);
   float GetCorrFactStep1(float E);


   ClassDef(h10,0);
};

#endif

#ifdef h10_cxx
void h10::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("npart", &npart, &b_npart);
   fChain->SetBranchAddress("evstat", &evstat, &b_evstat);
   fChain->SetBranchAddress("evntid", &evntid, &b_evntid);
   fChain->SetBranchAddress("evntype", &evntype, &b_evntype);
   fChain->SetBranchAddress("evntclas", &evntclas, &b_evntclas);
   fChain->SetBranchAddress("q_l", &q_l, &b_q_l);
   fChain->SetBranchAddress("t_l", &t_l, &b_t_l);
   fChain->SetBranchAddress("tr_time", &tr_time, &b_tr_time);
   fChain->SetBranchAddress("rf_time", &rf_time, &b_rf_time);
   fChain->SetBranchAddress("l2bit", &l2bit, &b_l2bit);
   fChain->SetBranchAddress("l3bit", &l3bit, &b_l3bit);
   fChain->SetBranchAddress("hlsc", &hlsc, &b_hlsc);
   fChain->SetBranchAddress("intt", &intt, &b_intt);
   fChain->SetBranchAddress("gpart", &gpart, &b_gpart);
   fChain->SetBranchAddress("id", id, &b_id);
   fChain->SetBranchAddress("stat", stat, &b_stat);
   fChain->SetBranchAddress("dc", dc, &b_dc);
   fChain->SetBranchAddress("cc", cc, &b_cc);
   fChain->SetBranchAddress("sc", sc, &b_sc);
   fChain->SetBranchAddress("ec", ec, &b_ec);
   fChain->SetBranchAddress("lec", lec, &b_lec);
   fChain->SetBranchAddress("st", st, &b_st);
   fChain->SetBranchAddress("p", p, &b_p);
   fChain->SetBranchAddress("m", m, &b_m);
   fChain->SetBranchAddress("q", q, &b_q);
   fChain->SetBranchAddress("b", b, &b_b);
   fChain->SetBranchAddress("cx", cx, &b_cx);
   fChain->SetBranchAddress("cy", cy, &b_cy);
   fChain->SetBranchAddress("cz", cz, &b_cz);
   fChain->SetBranchAddress("vx", vx, &b_vx);
   fChain->SetBranchAddress("vy", vy, &b_vy);
   fChain->SetBranchAddress("vz", vz, &b_vz);
   fChain->SetBranchAddress("dc_part", &dc_part, &b_dc_part);
   fChain->SetBranchAddress("dc_sect", dc_sect, &b_dc_sect);
   fChain->SetBranchAddress("dc_trk", dc_trk, &b_dc_trk);
   fChain->SetBranchAddress("dc_stat", dc_stat, &b_dc_stat);
   fChain->SetBranchAddress("tb_st", tb_st, &b_tb_st);
   fChain->SetBranchAddress("dc_xsc", dc_xsc, &b_dc_xsc);
   fChain->SetBranchAddress("dc_ysc", dc_ysc, &b_dc_ysc);
   fChain->SetBranchAddress("dc_zsc", dc_zsc, &b_dc_zsc);
   fChain->SetBranchAddress("dc_cxsc", dc_cxsc, &b_dc_cxsc);
   fChain->SetBranchAddress("dc_cysc", dc_cysc, &b_dc_cysc);
   fChain->SetBranchAddress("dc_czsc", dc_czsc, &b_dc_czsc);
   fChain->SetBranchAddress("dc_vx", dc_vx, &b_dc_vx);
   fChain->SetBranchAddress("dc_vy", dc_vy, &b_dc_vy);
   fChain->SetBranchAddress("dc_vz", dc_vz, &b_dc_vz);
   fChain->SetBranchAddress("dc_vr", dc_vr, &b_dc_vr);
   fChain->SetBranchAddress("tl1_cx", tl1_cx, &b_tl1_cx);
   fChain->SetBranchAddress("tl1_cy", tl1_cy, &b_tl1_cy);
   fChain->SetBranchAddress("tl1_cz", tl1_cz, &b_tl1_cz);
   fChain->SetBranchAddress("tl1_x", tl1_x, &b_tl1_x);
   fChain->SetBranchAddress("tl1_y", tl1_y, &b_tl1_y);
   fChain->SetBranchAddress("tl1_z", tl1_z, &b_tl1_z);
   fChain->SetBranchAddress("tl1_r", tl1_r, &b_tl1_r);
   fChain->SetBranchAddress("dc_c2", dc_c2, &b_dc_c2);
   fChain->SetBranchAddress("ec_part", &ec_part, &b_ec_part);
   fChain->SetBranchAddress("ec_stat", ec_stat, &b_ec_stat);
   fChain->SetBranchAddress("ec_sect", ec_sect, &b_ec_sect);
   fChain->SetBranchAddress("ec_whol", ec_whol, &b_ec_whol);
   fChain->SetBranchAddress("ec_inst", ec_inst, &b_ec_inst);
   fChain->SetBranchAddress("ec_oust", ec_oust, &b_ec_oust);
   fChain->SetBranchAddress("etot", etot, &b_etot);
   fChain->SetBranchAddress("ec_ei", ec_ei, &b_ec_ei);
   fChain->SetBranchAddress("ec_eo", ec_eo, &b_ec_eo);
   fChain->SetBranchAddress("ec_t", ec_t, &b_ec_t);
   fChain->SetBranchAddress("ec_r", ec_r, &b_ec_r);
   fChain->SetBranchAddress("ech_x", ech_x, &b_ech_x);
   fChain->SetBranchAddress("ech_y", ech_y, &b_ech_y);
   fChain->SetBranchAddress("ech_z", ech_z, &b_ech_z);
   fChain->SetBranchAddress("ec_m2", ec_m2, &b_ec_m2);
   fChain->SetBranchAddress("ec_m3", ec_m3, &b_ec_m3);
   fChain->SetBranchAddress("ec_m4", ec_m4, &b_ec_m4);
   fChain->SetBranchAddress("ec_c2", ec_c2, &b_ec_c2);
   fChain->SetBranchAddress("sc_part", &sc_part, &b_sc_part);
   fChain->SetBranchAddress("sc_sect", sc_sect, &b_sc_sect);
   fChain->SetBranchAddress("sc_hit", sc_hit, &b_sc_hit);
   fChain->SetBranchAddress("sc_pd", sc_pd, &b_sc_pd);
   fChain->SetBranchAddress("sc_stat", sc_stat, &b_sc_stat);
   fChain->SetBranchAddress("edep", edep, &b_edep);
   fChain->SetBranchAddress("sc_t", sc_t, &b_sc_t);
   fChain->SetBranchAddress("sc_r", sc_r, &b_sc_r);
   fChain->SetBranchAddress("sc_c2", sc_c2, &b_sc_c2);
   fChain->SetBranchAddress("cc_part", &cc_part, &b_cc_part);
   fChain->SetBranchAddress("cc_sect", cc_sect, &b_cc_sect);
   fChain->SetBranchAddress("cc_hit", cc_hit, &b_cc_hit);
   fChain->SetBranchAddress("cc_segm", cc_segm, &b_cc_segm);
   fChain->SetBranchAddress("nphe", nphe, &b_nphe);
   fChain->SetBranchAddress("cc_t", cc_t, &b_cc_t);
   fChain->SetBranchAddress("cc_r", cc_r, &b_cc_r);
   fChain->SetBranchAddress("cc_c2", cc_c2, &b_cc_c2);
   fChain->SetBranchAddress("lac_part", &lac_part, &b_lac_part);
   fChain->SetBranchAddress("lec_sect", lec_sect, &b_lec_sect);
   fChain->SetBranchAddress("lec_hit", lec_hit, &b_lec_hit);
   fChain->SetBranchAddress("lec_stat", lec_stat, &b_lec_stat);
   fChain->SetBranchAddress("lec_etot", lec_etot, &b_lec_etot);
   fChain->SetBranchAddress("lec_ein", lec_ein, &b_lec_ein);
   fChain->SetBranchAddress("lec_t", lec_t, &b_lec_t);
   fChain->SetBranchAddress("lec_r", lec_r, &b_lec_r);
   fChain->SetBranchAddress("lec_x", lec_x, &b_lec_x);
   fChain->SetBranchAddress("lec_y", lec_y, &b_lec_y);
   fChain->SetBranchAddress("lec_z", lec_z, &b_lec_z);
   fChain->SetBranchAddress("lec_c2", lec_c2, &b_lec_c2);
   fChain->SetBranchAddress("vidmvrt", &vidmvrt, &b_vidmvrt);
   fChain->SetBranchAddress("ntrmvrt", &ntrmvrt, &b_ntrmvrt);
   fChain->SetBranchAddress("xmvrt", &xmvrt, &b_xmvrt);
   fChain->SetBranchAddress("ymvrt", &ymvrt, &b_ymvrt);
   fChain->SetBranchAddress("zmvrt", &zmvrt, &b_zmvrt);
   fChain->SetBranchAddress("ch2mvrt", &ch2mvrt, &b_ch2mvrt);
   fChain->SetBranchAddress("cxxmvrt", &cxxmvrt, &b_cxxmvrt);
   fChain->SetBranchAddress("cxymvrt", &cxymvrt, &b_cxymvrt);
   fChain->SetBranchAddress("cxzmvrt", &cxzmvrt, &b_cxzmvrt);
   fChain->SetBranchAddress("cyymvrt", &cyymvrt, &b_cyymvrt);
   fChain->SetBranchAddress("cyzmvrt", &cyzmvrt, &b_cyzmvrt);
   fChain->SetBranchAddress("stamvrt", &stamvrt, &b_stamvrt);
   fChain->SetBranchAddress("mcnentr", &mcnentr, &b_mcnentr);
   fChain->SetBranchAddress("mcnpart", &mcnpart, &b_mcnpart);
   fChain->SetBranchAddress("mcst", mcst, &b_mcst);
   fChain->SetBranchAddress("mcid", mcid, &b_mcid);
   fChain->SetBranchAddress("mcpid", mcpid, &b_mcpid);
   fChain->SetBranchAddress("mctheta", mctheta, &b_mctheta);
   fChain->SetBranchAddress("mcphi", mcphi, &b_mcphi);
   fChain->SetBranchAddress("mcp", mcp, &b_mcp);
   fChain->SetBranchAddress("mcm", mcm, &b_mcm);
   fChain->SetBranchAddress("mcvx", mcvx, &b_mcvx);
   fChain->SetBranchAddress("mcvy", mcvy, &b_mcvy);
   fChain->SetBranchAddress("mcvz", mcvz, &b_mcvz);
   fChain->SetBranchAddress("mctof", mctof, &b_mctof);
   fChain->SetBranchAddress("nprt", &nprt, &b_nprt);
   fChain->SetBranchAddress("pidpart", pidpart, &b_pidpart);
   fChain->SetBranchAddress("xpart", xpart, &b_xpart);
   fChain->SetBranchAddress("ypart", ypart, &b_ypart);
   fChain->SetBranchAddress("zpart", zpart, &b_zpart);
   fChain->SetBranchAddress("epart", epart, &b_epart);
   fChain->SetBranchAddress("pxpart", pxpart, &b_pxpart);
   fChain->SetBranchAddress("pypart", pypart, &b_pypart);
   fChain->SetBranchAddress("pzpart", pzpart, &b_pzpart);
   fChain->SetBranchAddress("qpart", qpart, &b_qpart);
   fChain->SetBranchAddress("Ipart10", Ipart10, &b_Ipart10);
   fChain->SetBranchAddress("Rpart11", Rpart11, &b_Rpart11);
   fChain->SetBranchAddress("Rpart12", Rpart12, &b_Rpart12);
   fChain->SetBranchAddress("Ipart13", Ipart13, &b_Ipart13);
   fChain->SetBranchAddress("nschit", &nschit, &b_nschit);
   fChain->SetBranchAddress("scsect", scsect, &b_scsect);
   fChain->SetBranchAddress("schid", schid, &b_schid);
   fChain->SetBranchAddress("scpid", scpid, &b_scpid);
   fChain->SetBranchAddress("sct", sct, &b_sct);
   fChain->SetBranchAddress("sce", sce, &b_sce);
   fChain->SetBranchAddress("scx", scx, &b_scx);
   fChain->SetBranchAddress("scy", scy, &b_scy);
   fChain->SetBranchAddress("scz", scz, &b_scz);
}

Bool_t h10::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef h10_cxx
