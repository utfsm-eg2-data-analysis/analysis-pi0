#define h10_cxx
#include "h10.h"
#include <TH2.h>
#include <TStyle.h>



	/***********************************  BELOW CODE IS SET TO RUN Fe !!!! *****************************/

void h10::Begin(TTree *tree)
{
   TString option = GetOption();
   Init(tree);

   if(AvailableEvents==0)AvailableEvents = Int_t(tree->GetEntries());
   cout << "***********************************************************" << endl;
   cout << "********* WELCOME IN EC PHOTON RESOLUTION PROGRAM *********" << endl;
   cout << "***********************************************************" << endl << endl;
   cout << "Activating prefered styles... "; cout.flush();
   gROOT->SetStyle("Plain");
   gStyle->SetPalette(1);
   gStyle->SetOptFit(1111);
   gStyle->SetOptStat(1111111);
   gStyle->SetOptTitle(1);
   cout << "Done !" << endl;

  
   TFidCut("Lorenzo_5GeV.txt");
   FileGam = new TFile("outFe.root","recreate");
  
   TreeGam = new TTree("ec_gam","TreeGSIM");
   MakeBranch();
   cout << "Done !" << endl;
  eND = 0;
   eNS = 0;
   Nevta = 0;
   Nevt_found = 0;
   cout << "looping" << endl;
   for(int i=0;i<100;i++)cout << "*";
   cout << endl;

}

void h10::SlaveBegin(TTree *tree)
{

   TString option = GetOption();
    Init(tree);
}

Bool_t h10::Process(Long64_t entry)
{

   fChain->GetTree()->GetEntry(entry);
   TString option = GetOption();
   Nevta++;

   Int_t evtsPrint = Int_t(TMath::Floor(Nevta/(AvailableEvents/100)));
   evtsPrint = evtsPrint*Int_t(AvailableEvents/100);
   if(Nevta==evtsPrint){cout << "*";cout.flush();}

   npi0 =0;   
   nphot=0;
  
  if(GoodTrig()){
   V4Elec.SetXYZT(TMath::Cos(TMath::DegToRad()*trig_phi)*TMath::Sin(TMath::DegToRad()*trig_theta)*trig_mom,
		  TMath::Sin(TMath::DegToRad()*trig_phi)*TMath::Sin(TMath::DegToRad()*trig_theta)*trig_mom, 
		  TMath::Cos(TMath::DegToRad()*trig_theta)*trig_mom, 
		  trig_mom);
   V4Beam.SetXYZT(0.,0.,5.014,5.014);
   V4Q = V4Beam-V4Elec;
   
    for(int part=1;part<gpart;part++){
     if(nphot<10&&IsAPhoton(part))MakePhoton(part);
    }
    SearchPi0s();
    TreeGam->Fill();
  } //Good Trig
  
   return kTRUE;
}

void h10::SlaveTerminate()
{

}

void h10::Terminate()
{

   cout << endl;
   TreeGam->Write();
   FileGam->Close();
   cout << "Finished. Ciao, baby!" << endl;

}

////////////////////Function Declaration////////////////////////


//////////////////// Electron ////////////////////

Bool_t h10::GoodTrig(void){
 Bool_t trig = kFALSE;
 Int_t dci, cci, sci, eci, sti;
 sti = stat[0];
 dci = dc[0];
 cci = cc[0];
 eci = ec[0];
 sci = sc[0];
 
 
 if(sti>0 && eci>0 && cci>0 && dci>0 &&sci>0){
     eci--;cci--;dci--;sti--;sci--;
     trig_mom = p[0];
     trig_theta_rad = TMath::ACos((p[0]*cz[0])/p[0]);
     trig_Q2= (4*5.014*(trig_mom)*sin(trig_theta_rad/2)*sin(trig_theta_rad/2));
     trig_nu=(5.014-trig_mom);
     trig_W=sqrt(0.938*0.938+2*0.938*(trig_nu)-trig_Q2);
     trig_y=trig_nu/5.014;     
    
     if(q[0]==-1
	&&dc_stat[dci]>0
        &&ec_stat[eci]>0
	&&sc_stat[sci]>0 
	&&trig_Q2>1
	&&trig_W>2
	){
       trig_ql = q_l;
       trig_sect = dc_sect[dci];
       trig_theta = TMath::RadToDeg()*trig_theta_rad; 
       trig_phi = TMath::RadToDeg()*TMath::ATan2(cy[0]*trig_mom,cx[0]*trig_mom);
       trig_nphe = nphe[cci];
       trig_ecin = ec_ei[eci];
       trig_ecout = ec_eo[eci];
       trig_ectot = etot[eci];
       trig_ece = TMath::Max(trig_ectot, trig_ecin+trig_ecout);
       trig_vx = vx[0];
       trig_vy = vy[0];
       trig_vz = vz[0];
       trig_px = cx[0]*trig_mom;
       trig_py = cy[0]*trig_mom;
       trig_pz = cz[0]*trig_mom;
       trig_ecx = ech_x[eci];
       trig_ecy = ech_y[eci];
       trig_ecz = ech_z[eci];
       TVector3 xyzEC(trig_ecx,trig_ecy,trig_ecz);
       TVector3 uvwEC = ec_xyz_uvw(xyzEC);
       trig_ecu = uvwEC.X();
       trig_ecv = uvwEC.Y();
       trig_ecw = uvwEC.Z();
       trig_ect = ec_t[eci];    
       trig_ecl = ec_t[eci]; 
       trig_sct = sc_t[sci];
       trig_sc_stat = sc_stat[sci];
       trig_scl = sc_r[sci];
       MirrorCode = cc_segm[cci]%1000;

       trig_phisect = trig_phi - (trig_phi>32.)*(trig_sect-1)*60.+(trig_phi<-32.)*(7-trig_sect)*60.; // all sectors dist. on one  image;  trig_phi: -180 to 180
       trig_dcfid_flag = 0;
       trig_tl1x = tl1_x[dci];
       trig_tl1y = tl1_y[dci];
       trig_tl1z = tl1_z[dci];
       trig_tl1cx = tl1_cx[dci];
       trig_tl1cy = tl1_cy[dci];
       trig_tl1cz = tl1_cz[dci];
       trig_thetal1 = TMath::RadToDeg()*TMath::ACos((trig_tl1z-trig_vz_corr)/sqrt(trig_tl1x*trig_tl1x+trig_tl1y*trig_tl1y+((trig_tl1z-trig_vz_corr)*(trig_tl1z-trig_vz_corr))));
       trig_phil1 = TMath::RadToDeg()*TMath::ATan2(trig_tl1y,trig_tl1x);
       while(trig_phil1> 30)trig_phil1-=60.;
       while(trig_phil1<-30)trig_phil1+=60.;

      LorenzoFlag = CheckCut();
      FillVlassovCoordinates(0);

      Float_t phi;
      phi = trig_phi;
      phi += 30.;
      while(phi<0)phi+=360.;
      sect =  (int)TMath::Floor(phi/60.);//trig_phisect;
     
      //Vertex correction (x,y,z)->(x_corr,y_corr,z_corr). Correction preserves (x,y) in coordinates of Sector1 
      TVector3 RotatedVertPos(trig_vx,trig_vy,trig_vz);
      TVector3 RotatedVertDir(trig_px,trig_py,trig_pz);
      TVector3 TargetPos(0.043,-0.33,0);
      RotatedVertPos.RotateZ(-TMath::DegToRad()*60.*sect);
      RotatedVertDir.RotateZ(-TMath::DegToRad()*60.*sect);
      TargetPos.RotateZ(-TMath::DegToRad()*60.*sect);
      Float_t ShiftLength = (TargetPos.X()-RotatedVertPos.X())/RotatedVertDir.X();
      RotatedVertDir = ShiftLength*RotatedVertDir;
      RotatedVertPos = RotatedVertPos+RotatedVertDir;
      trig_vx_corr =  (RotatedVertPos-TargetPos).X();
      trig_vy_corr =  (RotatedVertPos-TargetPos).Y();
      trig_vz_corr =  RotatedVertPos.Z();
      TVector3 TransvShift = RotatedVertPos-TargetPos;
      TransvShift.SetZ(0);
      trig_transv_pos = TransvShift.Mag();
      if(TransvShift.Y()<0)trig_transv_pos *= -1.;
      target=0;

      if (trig_vy_corr>-2.4 && trig_vy_corr<2 &&  trig_vz_corr>-25.7 && trig_vz_corr<-24.0)target=2; //SOLID  GENERAL
      if (trig_vy_corr>-2.2 && trig_vy_corr<2. && trig_vz_corr>-31.8 &&  trig_vz_corr<-28.4)target=1;   //DEUTERIUM
     // if (trig_vy_corr>-2.2 && trig_vy_corr<2. &&  trig_vz_corr>-25.3 && trig_vz_corr<-24.1)target=2; //CARBON
      if (trig_vy_corr>-2.2 && trig_vy_corr<2 &&  trig_vz_corr>-25.65 && trig_vz_corr<-24.26)target=2; //Fe
     // if (trig_vy_corr>-2.2 && trig_vy_corr<2 &&  trig_vz_corr>-25.54 && trig_vz_corr<-24.36)target=2; //Pb
    

      //  NpheSect(); not nessesary when applying Trig CC matching procedure
         if(kTRUE 
           &&trig_mom>0.75 	        
           &&ec_ei[eci]>0.06
           &&ec_ei[eci]*ec_eo[eci]>0
           &&TrigCCmatch()
       //  &&NpheSect() 
           &&SamplingCutFe()
           &&TMath::Abs(trig_ect-trig_sct-0.7)<5*0.35
	   &&trig_ecu<400 
	   &&trig_ecu>40
	   &&trig_ecw<390
	   &&trig_ecv<360 
           &&LorenzoFlag==1
	   ){
	   if(target==1) eND++;
           if(target==2) eNS++;
           trig=kTRUE;
	 }

      } 
  }
  return trig;
   }


////////////////////// Photon ///////////////////////
 
Bool_t h10::IsAPhoton(int track){
 Bool_t result=kFALSE;
 Int_t eci;
 eci = ec[track];
 if(eci>0){
    eci--;
 if(q[track]==0){ 
   gam_Eec[nphot] = TMath::Max(etot[eci], ec_ei[eci]+ec_eo[eci]);
   gam_bet[nphot] = b[track];
   gam_E[nphot] = p[track];
   gam_Px[nphot] = p[track]*cx[track]; 
   gam_Py[nphot] = p[track]*cy[track];
   gam_Pz[nphot] = p[track]*cz[track];  
   gam_ecx[nphot] = ech_x[eci];
   gam_ecy[nphot] = ech_y[eci];
   gam_ecz[nphot] = ech_z[eci];
   TVector3 xyzEC(gam_ecx[nphot],gam_ecy[nphot],gam_ecz[nphot]);
   TVector3 uvwEC = ec_xyz_uvw(xyzEC);
   gam_ecu[nphot] = uvwEC.X();
   gam_ecv[nphot] = uvwEC.Y();
   gam_ecw[nphot] = uvwEC.Z();
   gam_ect[nphot] = ec_t[eci];
   gam_ecpath[nphot]=ec_r[eci];
   gam_phi[nphot]=TMath::RadToDeg()*TMath::ATan2(gam_ecy[nphot],gam_ecx[nphot]);
   float_t Rt = sqrt(gam_ecx[nphot]*gam_ecx[nphot]+gam_ecy[nphot]*gam_ecy[nphot]);
   float_t R = sqrt(gam_ecx[nphot]*gam_ecx[nphot]
		   +gam_ecy[nphot]*gam_ecy[nphot]
		   +(trig_vz_corr-gam_ecz[nphot])*(trig_vz_corr-gam_ecz[nphot]));
   gam_theta[nphot]=TMath::RadToDeg()*TMath::ASin(Rt/R);  

   gam_Etrue[nphot] = gam_Eec[nphot]/0.273;
   gam_Pxtrue[nphot] = gam_Etrue[nphot]*TMath::Cos(TMath::DegToRad()*gam_phi[nphot])*TMath::Sin(TMath::DegToRad()*gam_theta[nphot]);
   gam_Pytrue[nphot] = gam_Etrue[nphot]*TMath::Sin(TMath::DegToRad()*gam_phi[nphot])*TMath::Sin(TMath::DegToRad()*gam_theta[nphot]);
   gam_Pztrue[nphot] = gam_Etrue[nphot]*                                             TMath::Cos(TMath::DegToRad()*gam_theta[nphot]);

   gam_Etrue_Corr[nphot] = gam_Etrue[nphot]/GetCorrFactStep1(gam_Etrue[nphot]);
   gam_Pxtrue_Corr[nphot]= gam_Etrue_Corr[nphot]*TMath::Cos(TMath::DegToRad()*gam_phi[nphot])*TMath::Sin(TMath::DegToRad()*gam_theta[nphot]);
   gam_Pytrue_Corr[nphot]= gam_Etrue_Corr[nphot]*TMath::Sin(TMath::DegToRad()*gam_phi[nphot])*TMath::Sin(TMath::DegToRad()*gam_theta[nphot]);
   gam_Pztrue_Corr[nphot]= gam_Etrue_Corr[nphot]*                                             TMath::Cos(TMath::DegToRad()*gam_theta[nphot]);

 if(kTRUE
     &&gam_ecu[nphot]<410 && gam_ecu[nphot]>40
     &&gam_ecv[nphot]<370
     &&gam_ecw[nphot]<410
     &&(gam_ect[nphot]-tr_time-(gam_ecpath[nphot]/30))>-2.2&&(gam_ect[nphot]-tr_time-(gam_ecpath[nphot]/30))<1.3
     &&gam_Etrue[nphot]>0.1  //ATTENTION Ecut > 0.1 BUT later on for R I cut it to gam_Etrue[nphot]>0.3
   
   ){ 
    result=kTRUE;
   }
  }
 }
 return result;
}

//////////////////////////////////////////////////////////////////
void h10::MakePhoton(int part){
nphot++;return;  
}
//////////////////////////////////////////////////////////////////
void h10::SortPhotons(void){
 for(int gi=0;gi<nphot-1;gi++){
  Float_t MaxEnergyFound=gam_Etrue_Corr[gi];
  Int_t MaxEInd=-1;
  for(int li=gi+1;li<nphot;li++){
   if(gam_Etrue_Corr[li]>MaxEnergyFound){
    MaxEnergyFound=gam_Etrue_Corr[li];
    MaxEInd=li;
   }
  }
  if(MaxEnergyFound>gam_Etrue_Corr[gi]){
   //cout << "Swapping " << gam_Etrue[gi] << " and " << gam_Etrue[MaxEInd] << endl;
   SwapPhotons(gi,MaxEInd);
   //cout << "Swapped " << gam_Etrue[gi] << " and " << gam_Etrue[MaxEInd] << endl;
  }
 }
 return;
}

//////////////////////////////////////////////////////////////////////
void h10::SwapPhotons(int p1, int p2){
 MySwap(gam_Eec+p1,gam_Eec+p2);
 MySwap(gam_bet+p1,gam_bet+p2);
 MySwap(gam_E+p1,gam_E+p2);
 MySwap(gam_Px+p1,gam_Px+p2);
 MySwap(gam_Py+p1,gam_Py+p2);
 MySwap(gam_Pz+p1,gam_Pz+p2);
 MySwap(gam_Etrue+p1,gam_Etrue+p2);
 MySwap(gam_Pxtrue+p1,gam_Pxtrue+p2);
 MySwap(gam_Pytrue+p1,gam_Pytrue+p2);
 MySwap(gam_Pztrue+p1,gam_Pztrue+p2);
 MySwap(gam_Etrue_Corr+p1,gam_Etrue_Corr+p2);
 MySwap(gam_Pxtrue_Corr+p1,gam_Pxtrue_Corr+p2);
 MySwap(gam_Pytrue_Corr+p1,gam_Pytrue_Corr+p2);
 MySwap(gam_Pztrue_Corr+p1,gam_Pztrue_Corr+p2);
 MySwap(gam_Px+p1,gam_Px+p2);
 MySwap(gam_Py+p1,gam_Py+p2);
 MySwap(gam_Pz+p1,gam_Pz+p2);
 MySwap(gam_ecx+p1,gam_ecx+p2);
 MySwap(gam_ecy+p1,gam_ecy+p2);
 MySwap(gam_ecz+p1,gam_ecz+p2);
 MySwap(gam_ecu+p1,gam_ecu+p2);
 MySwap(gam_ecv+p1,gam_ecv+p2);
 MySwap(gam_ecw+p1,gam_ecw+p2);
 MySwap(gam_ect+p1,gam_ect+p2);
 MySwap(gam_ecpath+p1,gam_ecpath+p2);
 MySwap(gam_theta+p1,gam_theta+p2);
 MySwap(gam_phi+p1,gam_phi+p2);
}

/////////////////////////////////////////////////////////////////////
template <class T> void h10::MySwap (T* a, T* b) { 
 T tmp = *a;
 *a = *b;
 *b = tmp; 
}
//////////////////////////////////////////////////////////////////
Bool_t h10::SearchPi0s(void){
  Bool_t result=kFALSE;
 if(nphot>1){
  npi0=0;
  SortPhotons();
  for(int paire=0;paire<nphot-1;paire++){
   for(int paire2=paire+1;paire2<nphot&&npi0<45;paire2++){
    
     //Pi0 Mass from EVNT bank 
    TLorentzVector VG1(gam_Px[paire],gam_Py[paire],gam_Pz[paire],gam_E[paire]);
    TLorentzVector VG2(gam_Px[paire2],gam_Py[paire2],gam_Pz[paire2],gam_E[paire2]);
    TLorentzVector VPi = VG1 + VG2;
    opening_angle_evnt[npi0] =(VG1.Vect()).Angle(VG2.Vect());
    opening_angle_evnt[npi0] *= TMath::RadToDeg();

    // if(opening_angle>5.&&TMath::Abs(VPi.M()-0.135)<0.1&&VPi.E()>0.5){ //recheck opening angle; Recsis analyses of two neutral hits gives 9 if in one sector
    //   Float_t this_open_angle_cut = 2.*TMath::RadToDeg()*TMath::ATan(0.135/(gam_Etrue[G1[ip]]+gam_Etrue[G2[ip]])); 

     Pi0mass_EVNT[npi0]=VPi.M(); 
     Pi0z_EVNT[npi0]= VPi.E()/trig_nu; 
     Pi0energy_EVNT[npi0]=VPi.E();
  
     //Pi0mass using ECPB bank only
     TLorentzVector Vg1(gam_Pxtrue[paire],  gam_Pytrue[paire],  gam_Pztrue[paire],  gam_Etrue[paire]);
     TLorentzVector Vg2(gam_Pxtrue[paire2], gam_Pytrue[paire2], gam_Pztrue[paire2], gam_Etrue[paire2]);
     V4Pi0 = Vg1 + Vg2;  
     opening_angle_ecpb[npi0] =(Vg1.Vect()).Angle(Vg2.Vect());
     opening_angle_ecpb[npi0] *= TMath::RadToDeg();
     Pi0mass_ECPB[npi0]=V4Pi0.M(); 
     Pi0z_ECPB[npi0]= V4Pi0.E()/trig_nu; 
     Pi0energy_ECPB[npi0]=V4Pi0.E();
     pT2_val[npi0] = calculate_pT2();
  
    //Pi0mass after energy correction
    VG1corr.SetXYZT(gam_Pxtrue_Corr[paire],  gam_Pytrue_Corr[paire],  gam_Pztrue_Corr[paire],  gam_Etrue_Corr[paire]);
    VG2corr.SetXYZT(gam_Pxtrue_Corr[paire2],  gam_Pytrue_Corr[paire2],  gam_Pztrue_Corr[paire2],  gam_Etrue_Corr[paire2]);
   
    V4Pi0corr = VG1corr+VG2corr;
    Pi0mass_Corr[npi0]=V4Pi0corr.M(); 
    Pi0z_Corr[npi0]= V4Pi0corr.E()/trig_nu; 
    Pi0energy_Corr[npi0]= V4Pi0corr.E();
    Pi0theta_Corr[npi0] = V4Pi0corr.Theta();
    Pi0theta_Corr[npi0]*= TMath::RadToDeg();
    Pi0phi_Corr[npi0] = V4Pi0corr.Phi();
    Pi0phi_Corr[npi0] *= TMath::RadToDeg();
    Pi0angle_gs[npi0] = (V4Pi0corr.Vect()).Angle(V4Q.Vect());
    Pi0angle_gs[npi0] *= TMath::RadToDeg();

    opening_angle_Corr[npi0] =(VG1corr.Vect()).Angle(VG2corr.Vect());
    opening_angle_Corr[npi0] *= TMath::RadToDeg();
    pT2_val_Corr[npi0] = calculate_pT2_corr();
    phi_val_Corr[npi0] = calculate_phi_corr();   
   
    eG1Angle[npi0] = calculate_eG1_angle();
    eG2Angle[npi0] = calculate_eG2_angle();

     G1[npi0] = paire;
     G2[npi0] = paire2;
     npi0++;
   }
   if(npi0>0)result=kTRUE;
  }
 }
 return result;
}

//////////////////////////////////////////////////////////////////////////

void h10::MakeBranch(void){
  if(TreeGam!=NULL){ 
  TreeGam->Branch("trig_ql",&trig_ql,"trig_ql/F");
  TreeGam->Branch("tr_time",&tr_time,"tr_time/F");
  TreeGam->Branch("evntid",&evntid,"evntid/i");
  TreeGam->Branch("trig_sect",&trig_sect,"trig_sect/I"); //sector where e detected
  TreeGam->Branch("sect",&sect,"sect/I"); //sector where e detected
  TreeGam->Branch("trig_mom",&trig_mom,"trig_mom/F");   //e mom
  TreeGam->Branch("trig_theta",&trig_theta,"trig_theta/F"); 
  TreeGam->Branch("trig_phi",&trig_phi,"trig_phi/F");
  TreeGam->Branch("trig_phisect",&trig_phisect,"trig_phisect/F"); //phi_e - phi_sect, phi_sect = 60*(sect-1) +\- 360
  TreeGam->Branch("trig_nphe",&trig_nphe,"trig_nphe/F");
  TreeGam->Branch("trig_ecin",&trig_ecin,"trig_ecin/F");
  TreeGam->Branch("trig_ecout",&trig_ecout,"trig_ecout/F");
  TreeGam->Branch("trig_ectot",&trig_ectot,"trig_ectot/F");
  TreeGam->Branch("trig_ece",&trig_ece,"trig_ece/F"); // e energy
  TreeGam->Branch("trig_vz",&trig_vz,"trig_vz/F");
  TreeGam->Branch("trig_vx",&trig_vx,"trig_vx/F");
  TreeGam->Branch("trig_vy",&trig_vy,"trig_vy/F");
  TreeGam->Branch("trig_vx_corr",&trig_vx_corr,"trig_vx_corr/F");
  TreeGam->Branch("trig_vy_corr",&trig_vy_corr,"trig_vy_corr/F");
  TreeGam->Branch("trig_vz_corr",&trig_vz_corr,"trig_vz_corr/F");
  TreeGam->Branch("trig_px",&trig_px,"trig_px/F");
  TreeGam->Branch("trig_py",&trig_py,"trig_py/F");
  TreeGam->Branch("trig_pz",&trig_pz,"trig_pz/F");
  TreeGam->Branch("trig_ecx",&trig_ecx,"trig_ecx/F");
  TreeGam->Branch("trig_ecy",&trig_ecy,"trig_ecy/F");
  TreeGam->Branch("trig_ecz",&trig_ecz,"trig_ecz/F");
  TreeGam->Branch("trig_ecu",&trig_ecu,"trig_ecu/F");
  TreeGam->Branch("trig_ecv",&trig_ecv,"trig_ecv/F");
  TreeGam->Branch("trig_ecw",&trig_ecw,"trig_ecw/F");
  TreeGam->Branch("trig_ect",&trig_ect,"trig_ect/F");//ec time
  TreeGam->Branch("trig_ecl",&trig_ecl,"trig_ecl/F");
  TreeGam->Branch("trig_sct",&trig_sct,"trig_sct/F");// e tof
  TreeGam->Branch("trig_scl",&trig_scl,"trig_scl/F");
  TreeGam->Branch("trig_sc_stat",&trig_sc_stat,"trig_sc_stat/I");
  TreeGam->Branch("trig_cctheta",&trig_cctheta,"trig_cctheta/F"); 
  TreeGam->Branch("trig_ccphi",&trig_ccphi,"trig_ccphi/F");
  TreeGam->Branch("trig_tl1x",&trig_tl1x,"trig_tl1x/F");
  TreeGam->Branch("trig_tl1y",&trig_tl1y,"trig_tl1y/F");
  TreeGam->Branch("trig_tl1z",&trig_tl1z,"trig_tl1z/F");
  TreeGam->Branch("trig_tl1cx",&trig_tl1cx,"trig_tl1cx/F");
  TreeGam->Branch("trig_tl1cy",&trig_tl1cy,"trig_tl1cy/F");
  TreeGam->Branch("trig_tl1cz",&trig_tl1cz,"trig_tl1cz/F");
  TreeGam->Branch("trig_phil1",&trig_phil1,"trig_phil1/F");
  TreeGam->Branch("trig_thetal1",&trig_thetal1,"trig_thetal1/F");
  TreeGam->Branch("trig_W",&trig_W,"trig_W/F");
  TreeGam->Branch("trig_Q2",&trig_Q2,"trig_Q2/F");
  TreeGam->Branch("trig_nu",&trig_nu,"trig_nu/F");
  TreeGam->Branch("trig_y",&trig_y,"trig_y/F");
  TreeGam->Branch("target",&target,"target/I");  
  TreeGam->Branch("trig_transv_pos",&trig_transv_pos,"trig_transv_pos/F");
  TreeGam->Branch("LorenzoFlag",&LorenzoFlag,"LorenzoFlag/I");  
  TreeGam->Branch("NpheCutSect",&NpheCutSect,"NpheCutSect/O");  
  TreeGam->Branch("trig_dcfid_flag",&trig_dcfid_flag,"trig_dcfid_flag/I");
  TreeGam->Branch("MirrorCode",&MirrorCode,"MirrorCode/I");
    
  TreeGam->Branch("nphot",&nphot,"nphot/I");
  TreeGam->Branch("gam_Eec",gam_Eec,"gam_Eec[nphot]/F");//photon energy in ec
  TreeGam->Branch("gam_bet",gam_bet,"gam_bet[nphot]/F");//beta
  TreeGam->Branch("gam_Pxtrue",gam_Pxtrue,"gam_Pxtrue[nphot]/F");
  TreeGam->Branch("gam_Pytrue",gam_Pytrue,"gam_Pytrue[nphot]/F");
  TreeGam->Branch("gam_Pztrue",gam_Pztrue,"gam_Pztrue[nphot]/F");
  TreeGam->Branch("gam_E",gam_E,"gam_E[nphot]/F");
  TreeGam->Branch("gam_Px",gam_Px,"gam_Px[nphot]/F");
  TreeGam->Branch("gam_Py",gam_Py,"gam_Py[nphot]/F");
  TreeGam->Branch("gam_Pz",gam_Pz,"gam_Pz[nphot]/F");
  TreeGam->Branch("gam_Etrue",gam_Etrue,"gam_Etrue[nphot]/F");
  TreeGam->Branch("gam_Pxtrue",gam_Pxtrue,"gam_Pxtrue[nphot]/F");
  TreeGam->Branch("gam_Pytrue",gam_Pytrue,"gam_Pytrue[nphot]/F");
  TreeGam->Branch("gam_Pztrue",gam_Pztrue,"gam_Pztrue[nphot]/F");
  TreeGam->Branch("gam_Etrue_Corr",gam_Etrue_Corr,"gam_Etrue_Corr[nphot]/F");
  TreeGam->Branch("gam_Pxtrue_Corr",gam_Pxtrue_Corr,"gam_Pxtrue_Corr[nphot]/F");
  TreeGam->Branch("gam_Pytrue_Corr",gam_Pytrue_Corr,"gam_Pytrue_Corr[nphot]/F");
  TreeGam->Branch("gam_Pztrue_Corr",gam_Pztrue_Corr,"gam_Pztrue_Corr[nphot]/F");
  TreeGam->Branch("gam_ecx",gam_ecx,"gam_ecx[nphot]/F");
  TreeGam->Branch("gam_ecy",gam_ecy,"gam_ecy[nphot]/F");
  TreeGam->Branch("gam_ecz",gam_ecz,"gam_ecz[nphot]/F");
  TreeGam->Branch("gam_ecu",gam_ecu,"gam_ecu[nphot]/F");
  TreeGam->Branch("gam_ecv",gam_ecv,"gam_ecv[nphot]/F");
  TreeGam->Branch("gam_ecw",gam_ecw,"gam_ecw[nphot]/F");
  TreeGam->Branch("gam_ect",gam_ect,"gam_ect[nphot]/F");
  TreeGam->Branch("gam_ecpath",gam_ecpath,"gam_ecpath[nphot]/F");
  TreeGam->Branch("gam_phi",gam_phi,"gam_phi[nphot]/F");
  TreeGam->Branch("gam_theta",gam_theta,"gam_theta[nphot]/F");
  
  TreeGam->Branch("npi0",&npi0,"npi0/I");
  TreeGam->Branch("Pi0z_EVNT", Pi0z_EVNT, "Pi0z_EVNT[npi0]/F");
  TreeGam->Branch("Pi0mass_EVNT",Pi0mass_EVNT,"Pi0mass_EVNT[npi0]/F");
  TreeGam->Branch("opening_angle_evnt",opening_angle_evnt,"opening_angle_evnt[npi0]/F");
  TreeGam->Branch("Pi0energy_EVNT", Pi0energy_EVNT, "Pi0energy_EVNT[npi0]/F");  
  TreeGam->Branch("Pi0z_ECPB", Pi0z_ECPB, "Pi0z_ECPB[npi0]/F");
  TreeGam->Branch("Pi0mass_ECPB",Pi0mass_ECPB,"Pi0mass_ECPB[npi0]/F");
  TreeGam->Branch("opening_angle_ecpb",opening_angle_ecpb,"opening_angle_ecpb[npi0]/F");
  TreeGam->Branch("Pi0energy_ECPB", Pi0energy_ECPB, "Pi0energy_ECPB[npi0]/F");    
  TreeGam->Branch("Pi0mass_Corr",Pi0mass_Corr,"Pi0mass_Corr[npi0]/F");
  TreeGam->Branch("Pi0z_Corr", Pi0z_Corr, "Pi0z_Corr[npi0]/F");
  TreeGam->Branch("Pi0energy_Corr", Pi0energy_Corr, "Pi0energy_Corr[npi0]/F");
  TreeGam->Branch("Pi0pheta_Corr", Pi0theta_Corr, "Pi0theta_Corr[npi0]/F"); 
  TreeGam->Branch("Pi0phi_Corr", Pi0phi_Corr, "Pi0phi_Corr[npi0]/F");
  TreeGam->Branch("Pi0angle_gs", Pi0angle_gs, "Pi0angle_gs[npi0]/F");
  TreeGam->Branch("pT2_val_Corr",pT2_val_Corr,"pT2_val_Corr[npi0]/F");
  TreeGam->Branch("phi_val_Corr",phi_val_Corr,"phi_val_Corr[npi0]/F");
  TreeGam->Branch("pT1_val",pT2_val,"pT2_val[npi0]/F");
  TreeGam->Branch("opening_angle_Corr",opening_angle_Corr,"opening_angle_Corr[npi0]/F");
  TreeGam->Branch("G1",G1,"G1[npi0]/I");
  TreeGam->Branch("G2",G2,"G2[npi0]/I");
  TreeGam->Branch("eG1Angle", eG1Angle,"eG1Angle[npi0]/F");
  TreeGam->Branch("eG2Angle", eG2Angle,"eG2Angle[npi0]/F");
  //TreeGam->Branch("",,"");
  // TreeGam->Branch(",", ,"[npi0]/F");
 
}
 return;
}

//////////////////AUXILARY FNC//////////////////////
float h10::calculate_eG1_angle(){
float_t result = 0;
TVector3 elec =  V4Elec.Vect();
TVector3 phot1 = VG1corr.Vect();
 result = TMath::RadToDeg()*(elec.Angle(phot1));
 return result;
}

float h10::calculate_eG2_angle(){
float_t result = 0;
TVector3 elec =  V4Elec.Vect();
TVector3 phot2 = VG2corr.Vect();
 result = TMath::RadToDeg()*(elec.Angle(phot2));
 return result;
}

float h10::calculate_pT2_corr(){
  TVector3 qmom = V4Q.Vect();
  float result =0;
  TVector3 pimom = V4Pi0corr.Vect();
  TVector3 pimomlong;
  Float_t dot_prod = pimom.Dot(qmom);
  Float_t qnorm = qmom.Mag();
  if(qnorm>0){
    pimomlong = qmom;
    pimomlong.SetMag(dot_prod/qnorm);
    result = (pimom-pimomlong).Mag2();
  }
  return result;
}
float h10::calculate_phi_corr(){
 float result = -99.;
 //cout << "Go again " << V4Elec.E() << " theta=" << TMath::RadToDeg()*V4Elec.Theta() << " phi=" << TMath::RadToDeg()*V4Elec.Phi() << endl;
 TVector3 qmom = V4Q.Vect();
 TVector3 pimom = V4Pi0corr.Vect();
 TVector3 scat = V4Elec.Vect();
 TVector3 Lepto = scat.Cross(qmom);
 TVector3 Hadro = qmom.Cross(pimom);
 result = TMath::RadToDeg()*Lepto.Angle(Hadro);
 if(pimom.Dot(Lepto)<0)result=360.-result;
 return result;
}

float h10::calculate_pT2(){
  TVector3 qmom = V4Q.Vect();
  float result =0;
  TVector3 pimom = V4Pi0.Vect();
  TVector3 pimomlong;
  Float_t dot_prod = pimom.Dot(qmom);
  Float_t qnorm = qmom.Mag();
  if(qnorm>0){
    pimomlong = qmom;
    pimomlong.SetMag(dot_prod/qnorm);
    result = (pimom-pimomlong).Mag2();
  }
  return result;
}





//////////////////////////Electron Functions/////////////////////////////////


//********************************UVW****************************************

TVector3 h10::ec_xyz_uvw(TVector3 xyz){
  // Converts x,y,z EC hit in CLAS coordinate system
  // into u,v,w distances of the EC hit.

  Float_t ex=0.;
  Float_t wy=0.;
  Float_t zd=0.;
  Float_t yu=0.;
  Float_t ve=0.;
  Float_t wu=0.;
  Float_t xi=0.; 
  Float_t yi=0.; 
  Float_t zi=0.;
  Float_t ec_phy = 0.;  
  Float_t phy = 0.;
  Float_t rot[3][3];

  // Parameters
  Float_t ec_the = 0.4363323;
  Float_t ylow = -182.974;
  Float_t yhi = 189.956;
  Float_t tgrho = 1.95325; 
  Float_t sinrho = 0.8901256; 
  Float_t cosrho = 0.455715;

  // Variables
  ex = xyz[0];
  wy = xyz[1];
  zd = xyz[2];
  
  phy = TMath::ATan2(wy,ex)*57.29578;
  if(phy<0.){phy = phy + 360;}
  phy = phy+30.;
  if(phy>360.){phy = phy-360.;}

  ec_phy = ((Int_t) (phy/60.))*1.0471975;

  rot[0][0] = TMath::Cos(ec_the)*TMath::Cos(ec_phy);
  rot[0][1] = -TMath::Sin(ec_phy);
  rot[0][2] = TMath::Sin(ec_the)*TMath::Cos(ec_phy);
  rot[1][0] = TMath::Cos(ec_the)*TMath::Sin(ec_phy);
  rot[1][1] = TMath::Cos(ec_phy);
  rot[1][2] = TMath::Sin(ec_the)*TMath::Sin(ec_phy);
  rot[2][0] = -TMath::Sin(ec_the);
  rot[2][1] = 0.;
  rot[2][2] = TMath::Cos(ec_the);

  yi = ex*rot[0][0]+wy*rot[1][0]+zd*rot[2][0];
  xi = ex*rot[0][1]+wy*rot[1][1]+zd*rot[2][1];
  zi = ex*rot[0][2]+wy*rot[1][2]+zd*rot[2][2];
  zi = zi-510.32 ;

  yu = (yi-ylow)/sinrho;
  ve = (yhi-ylow)/tgrho - xi + (yhi-yi)/tgrho;
  wu = ((yhi-ylow)/tgrho + xi + (yhi-yi)/tgrho)/2./cosrho;

  TVector3 result3(yu,ve,wu);

  return result3;
}


//****************************Sampling fraction cut**************************

Bool_t h10::SamplingCutC(void){
  Bool_t result=kFALSE;

 if(trig_sect==1)result = TMath::Abs(trig_ece/trig_mom -( 0.252164 + 0.0122263 *trig_mom -0.000793937 *trig_mom*trig_mom))                                                                                     
                            <2.5*TMath::Sqrt( 9.55113e-03*9.55113e-03 + 3.40672e-02 *3.40672e-02  /trig_mom);
 else if(trig_sect==2)result  = TMath::Abs(trig_ece/trig_mom -( 0.278574 +0.0187482*trig_mom  -0.00238217*trig_mom*trig_mom))
			    <2.5*TMath::Sqrt( 1.39889e-02*1.39889e-02 +  3.74682e-02* 3.74682e-02 /trig_mom);
 else if(trig_sect==3)result = TMath::Abs(trig_ece/trig_mom -(0.262079 +   0.0230685*trig_mom  -0.00354741*trig_mom*trig_mom))
			      <2.5*TMath::Sqrt( 9.32762e-03*9.32762e-03 + 2.90046e-02*2.90046e-02 /trig_mom);
 else if(trig_sect==4)result = TMath::Abs(trig_ece/trig_mom -( 0.251108 +  0.0201568*trig_mom  -0.00332367*trig_mom*trig_mom))
				<2.5*TMath::Sqrt( 8.21055e-03*8.21055e-03 + 2.98893e-02*2.98893e-02 /trig_mom);
 else if(trig_sect==5)result = TMath::Abs(trig_ece/trig_mom -( 0.263396 +  0.00955238*trig_mom  -0.00102038*trig_mom*trig_mom))
				  <2.5*TMath::Sqrt( 2.25684e-02*2.25684e-02 + 3.06508e-02*3.06508e-02 /trig_mom);
 else if(trig_sect==6)result = TMath::Abs(trig_ece/trig_mom -(0.255245  + 0.0232659 *trig_mom  -0.00304798*trig_mom*trig_mom))
				    <2.5*TMath::Sqrt( 1.17254e-02*1.17254e-02 + 3.64221e-02*3.64221e-02 /trig_mom);
 return result;

}
 

Bool_t h10::SamplingCutFe(void){
  Bool_t result=kFALSE;
 if (trig_sect==1)result = TMath::Abs(trig_ece/trig_mom -(0.222404  +0.0222688 *trig_mom  -0.0024153*trig_mom*trig_mom))   
            		    <2.5*TMath::Sqrt( 9.23027e-03*9.23027e-03 + 2.98343e-02*2.98343e-02 /trig_mom);
 else if(trig_sect==2)result = TMath::Abs(trig_ece/trig_mom -(0.234623 +0.0194985*trig_mom  -0.00208357*trig_mom*trig_mom))
			<2.5*TMath::Sqrt( 8.66367e-03*8.66367e-03 + 3.08858e-02*3.08858e-02/trig_mom);
 else if(trig_sect==3)result =TMath::Abs(trig_ece/trig_mom -(0.252287 +0.024248*trig_mom  -0.00338846*trig_mom*trig_mom))
			<2.5*TMath::Sqrt( 1.07826e-02*1.07826e-02 + 2.63854e-02*2.63854e-02/trig_mom);
 else if(trig_sect==4)result =TMath::Abs(trig_ece/trig_mom -( 0.250946 +0.0208409*trig_mom -0.00326824*trig_mom*trig_mom))
			<2.5*TMath::Sqrt( 7.22581e-03*7.22581e-03 + 2.98809e-02*2.98809e-02 /trig_mom);
 else if(trig_sect==5)result = TMath::Abs(trig_ece/trig_mom -(0.271956 +0.0118487*trig_mom  -0.00187084*trig_mom*trig_mom))
			<2.5*TMath::Sqrt( 1.84073e-02*1.84073e-02 +  3.48029e-02* 3.48029e-02/trig_mom);
 else if(trig_sect==6)result = TMath::Abs(trig_ece/trig_mom -(0.252613 + 0.022819*trig_mom  -0.00311242*trig_mom*trig_mom))
			<2.5*TMath::Sqrt(4.11461e-03 *4.11461e-03 + 3.55081e-02*3.55081e-02 /trig_mom);
 return result;

}
 

Bool_t h10::SamplingCutPb(void){
  Bool_t result=kFALSE;

  if (trig_sect==1)result = TMath::Abs(trig_ece/trig_mom -(0.253431+0.0138251*trig_mom -0.0014016*trig_mom*trig_mom))
		     <2.5*TMath::Sqrt( 7.67408e-03*7.67408e-03 + 3.54391e-02*3.54391e-02 /trig_mom);
  else if(trig_sect==2)result = TMath::Abs(trig_ece/trig_mom -(0.249059 +  0.0147784 *trig_mom  -0.00148693*trig_mom*trig_mom))
                            <2.5*TMath::Sqrt( 7.52798e-03*7.52798e-03 + 3.38371e-02*3.38371e-02 /trig_mom);
  else if(trig_sect==3)result = TMath::Abs(trig_ece/trig_mom -(0.254573 +   0.022589*trig_mom  -0.00305686*trig_mom*trig_mom))
			    <2.5*TMath::Sqrt( 8.13241e-03*8.13241e-03 + 2.77300e-02*2.77300e-02 /trig_mom);
  else if(trig_sect==4)result = TMath::Abs(trig_ece/trig_mom -(0.255589 +  0.0190419 *trig_mom   -0.00305263*trig_mom*trig_mom))
			    <2.5*TMath::Sqrt( 7.20303e-03*7.20303e-03 + 3.03627e-02*3.03627e-02/trig_mom);
  else if(trig_sect==5)result = TMath::Abs(trig_ece/trig_mom -(0.276739+   0.0111585*trig_mom  -0.00175784*trig_mom*trig_mom))
			    <2.5*TMath::Sqrt( 1.80841e-02*1.80841e-02 + 3.53020e-02 *3.53020e-02  /trig_mom);
  else if(trig_sect==6)result = TMath::Abs(trig_ece/trig_mom -(0.262587 +   0.0191659*trig_mom  -0.0026264*trig_mom*trig_mom))
			    <2.5*TMath::Sqrt( 1.99220e-03*1.99220e-03 + 3.76172e-02*3.76172e-02 /trig_mom);
  return result;

}

//*****************************Nphe cut**************************************

Bool_t h10::NpheSect(void){
 Bool_t result=kFALSE;
 if(trig_sect==1||trig_sect==2)result = trig_nphe>25.; 
 else if(trig_sect==3)result = trig_nphe>26.;
 else if(trig_sect==4)result = trig_nphe>21.;
 else if(trig_sect==5||trig_sect==6)result = trig_nphe>28.; 
 else cout << "NpheSect could not find the sector ! " << trig_sect << endl;
 NpheCutSect = result;
 return result;
}


//*****************************CC match**************************************

Bool_t h10::TrigCCmatch(void){
 Bool_t result = 
       TMath::Abs(
       trig_cctheta-(7.306737e+00+1.383339e-01*MirrorCode+3.672767e-04*MirrorCode*MirrorCode)
                 )<2.5;
 return result;
}


// ******************************Fiducial cuts********************************        


Double_t fitf5(Double_t *y, Double_t *par2)
{
  // par2[4] is 1 for the upper part (p1) and -1 for the lower (p0)
  // par2[5] is the sector number

  
  Double_t arg = 0., fitval = 0.;
  arg = y[0];
  if(arg <= par2[2] || arg >= par2[3]) fitval = 60.*(par2[5] - 1); 
  else fitval = 60.*(par2[5] - 1) + par2[4]*par2[0]*(1. - 1./((arg - par2[2])/(par2[1]) + 1.));
 
  return fitval;
}



Int_t h10::CheckCut(void){ 

  Double_t moment = trig_mom;  //momentum
  Double_t phi = trig_phi;     //degrees
  Double_t theta = trig_theta; //degrees
  Double_t  theta_low=0.001;


  fcfidc_fun[0] = new TF1("func2",fitf5,12,60,6);
  fcfidc_fun[1] = new TF1("func3",fitf5,12,60,6);

  if (phi < -30) {
    phi+=360;
  }
  if (phi > 330) {
    phi-=360;
  }
  
 Int_t sector =  (Int_t)(phi + 90) /60;


  for( Int_t side = 0 ; side < 2 ; side++) {
    for(Int_t paraf = 0 ; paraf < 3 ; paraf++) {
      if(side==0) {
	fcfidc_fun[side]->SetParameter(paraf,fcfun[sector][side][paraf]->Eval(moment));
	if(paraf==2) {
         theta_low=fcfun[sector][side][paraf]->Eval(moment);
	}
      }
      else if(side==1 && paraf==2) {
	fcfidc_fun[side]->SetParameter(paraf,theta_low);
      }
      else {
	fcfidc_fun[side]->SetParameter(paraf,fcfun[sector][side][paraf]->Eval(moment));
      }
    }

    fcfidc_fun[side]->SetParameter(3,54); // theta max
    fcfidc_fun[side]->SetParameter(4, side*2 -1);//fixed: upper-lower switch
    fcfidc_fun[side]->SetParameter(5, sector);//fixed: sector switch

    
  }


  if (theta > theta_low && theta < 54 && phi > fcfidc_fun[0]->Eval(theta) && phi < fcfidc_fun[1]->Eval(theta) ) {
    return 1; // Fiducial Cut passed
  }
  else {
    return 0; // Fiducial Cut not passed
  }
}


void h10::TFidCut(Char_t* FileName){
  cout << "Executing constructor" << endl;
 
  Double_t fidc_par[7][2][4][6] ;
  Int_t sector;
  Int_t side;
  Int_t paraf;
  Char_t tmp[10];
  Char_t Title[40];
  float p0, p1, p2, p3, p4, p5;

  cout << "Copying name " << endl;
  fInFileLorenzo = TString( FileName );


  cout << "Opening file " << fInFileLorenzo.Data() <<endl;
  ifstream inputfile;
  inputfile.open(fInFileLorenzo.Data(), ifstream::in);
  if( !inputfile ) {
    cerr << "Error opening input stream" << endl;
    return;
  }              

  cout << "Reading files " << fInFileLorenzo.Data() << endl; 

 
  
  
  inputfile >> tmp  >> tmp >> tmp  >> tmp >> tmp  >> tmp >> tmp >> tmp  >> tmp >> tmp ;
  while( !inputfile.eof() ){
 
    inputfile >> sector >> side >> paraf >> p0 >> p1 >> p2 >> p3 >> p4 >> p5;

    
    fidc_par[sector][side][paraf][0] = p0;
    fidc_par[sector][side][paraf][1] = p1 ;
    fidc_par[sector][side][paraf][2] = p2;
    fidc_par[sector][side][paraf][3] = p3;
    fidc_par[sector][side][paraf][4] = p4;
    fidc_par[sector][side][paraf][5] = p5;

  }
  inputfile.close();
  
  for(sector = 1 ; sector < 7 ; sector++) {
    for(side = 0 ; side < 2 ; side ++) {
      for(paraf = 0 ; paraf < 3 ; paraf++) {
	sprintf(Title,"fidc_fun_se%i_si%i_p_%i",sector,side,paraf);
	if (paraf==0) {
	  fcfun[sector][side][paraf] = new TF1(Title,"[0]+[1]*exp([2]*(x-[3]))",0.2,6);
	}
	else if (paraf==1) {
	  fcfun[sector][side][paraf] = new TF1(Title,"[0]+[1]*x*exp([2]*(x-[3])**2)",0.2,6);
	}
	else if (paraf==2) {
	  fcfun[sector][side][paraf] = new TF1(Title,"[0]+[1]/x**2+[2]*x+[3]/x+[4]*exp([5]*x)",0.2,6);
	}
	else {
	  fcfun[sector][side][paraf] = new TF1(Title,"[0]+[1]*x+[2]/x+[3]*x**2",0.2,6);
	}
	
	
	fcfun[sector][side][paraf]->SetParameter(0,fidc_par[sector][side][paraf][0]);
	fcfun[sector][side][paraf]->SetParameter(1,fidc_par[sector][side][paraf][1]);
	fcfun[sector][side][paraf]->SetParameter(2,fidc_par[sector][side][paraf][2]);
	fcfun[sector][side][paraf]->SetParameter(3,fidc_par[sector][side][paraf][3]);
	if (paraf==2) {
          fcfun[sector][side][paraf]->SetParameter(4,fidc_par[sector][side][paraf][4]);
          fcfun[sector][side][paraf]->SetParameter(5,fidc_par[sector][side][paraf][5]);
	}
      }
    }
  }
}




//**********************CC Special Coordinates****************************

void h10:: FillVlassovCoordinates(int track){

  float point[3], dir[3], cc_coor[2];
  point[0] = dc_xsc[dc[track]-1];
  point[1] = dc_ysc[dc[track]-1];
  point[2] = dc_zsc[dc[track]-1];
    dir[0] = dc_cxsc[dc[track]-1];
    dir[1] = dc_cysc[dc[track]-1];
    dir[2] = dc_czsc[dc[track]-1];
  get_cc_special_coord(point,dir,cc_coor);
  trig_cctheta = cc_coor[0];
  trig_ccphi = cc_coor[1];


}

void h10::get_cc_special_coord(float *point, float *dir, float *cc_coor){
/*
C----------------------------------------------------------------------
C-  This macro was extracted from $CLAS_PACK/cc/cc_stuff.c
C-
C-  Input : point[3] and dir[3]
C-  -------
C-
C-          point[3] - coordinates of the point on the track
C-                     somewhere after CC (where B=0):
C-                     (it could be EC or SC  matching point)
C-          dir[3]   - direction vector at that point.
C-
C-  *** IMPORTANT ! *** Point and direction must be in Sector RS
C-
C-  Example : one can call with in particular (for track num 0, trigger) :
C-  ---------
C-
C-          point[0] = dc_xsc[dc[0]-1];
C-          point[1] = dc_ysc[dc[0]-1];
C-          point[2] = dc_zsc[dc[0]-1];
C-            dir[0] = dc_cxsc[dc[0]-1];
C-            dir[1] = dc_cysc[dc[0]-1];
C-            dir[2] = dc_czsc[dc[0]-1];
C-
C-  Ouput : cc_coor[2]
C   -------
C-          cc_coor[0] = theta
C-          cc_coor[1] = phy
C-          projective angles from the target to the track cross-point with "special plane".
C-
C----------------------------------------------------------------------
*/
  int nt;
  float x[3], dist, r, s, ccTheta, ccPhi;
  float cc_pln[3] = { - 0.0007840784063, 0., - 0.001681461571 }; //was static in original code
  nt = get_vcrpl(point, dir, cc_pln, &dist, x);
  ccTheta = -999.;
  ccPhi = -999.;
  if(nt){
   r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
   s = sqrt(x[0]*x[0] + x[1]*x[1]);
   ccTheta = TMath::RadToDeg()*TMath::ACos(x[2]/r);
   ccPhi = TMath::RadToDeg()*TMath::ATan2(x[1]/s,x[0]/s);
  }
  cc_coor[0] = ccTheta;
  cc_coor[1] = ccPhi;
 }//source for void get_cc_special_coord(float *, float *, float *)
                                                                                                                         
int h10::get_vcrpl(float *r0, float *dir, float *plane_par, float *dist, float *cross_point){
/*
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : crossing of the stright line(R0,d)
C-                         with a plane
C-
C-   Inputs  :   r0(3) - initial point of line
C-               dir(3) - vector direction: r = R0 + s*D
C-               plane_par(3) - array of plane parameters:
C-      plane_par(1)*x + plane_par(2)*y + plane_par(3)*z + 1 = 0
C-
C-   Outputs :   cc_vcrpl =  0 - no cross with the plane.
C-                           1 - cross in positive direction
C-                          -1 - cross in negative direction
C-               dist    =  Distance to the cross point
C-               cross_point(3) =  Cross point coordinates.
C-
C-   Created    23-NOV-1998   Alexander V. Vlassov
C-   Modified
C-
C----------------------------------------------------------------------
                                                                                                                         
*/
                                                                                                                         
  double a,b,t,c,d[3];
  const double un = 1.0000000000;
  const float vsmall = 0.000001;
  int i, ires;
  a = b = c = 0.;
  *dist = 0.;
  for(i=0;i<3;i++)
    c += un*dir[i]*dir[i];
  c = sqrt(c);
  if(c <= vsmall)
    {
        ires = 0;
        return(ires);
    }
  for(i=0;i<3;i++)
    d[i] = un*dir[i]/c;
  for(i=0;i<3;i++)
    { a += un*plane_par[i]*d[i]; b += un*plane_par[i]*r0[i]; }
  b += un;
  if(fabs(b) <= vsmall)
    {
    for(i=0;i<3;i++)
      cross_point[i] = r0[i];
    ires = 1;
    }
  else
    {
    if(fabs(a) <= vsmall)
      {
      for(i=0;i<3;i++)
        cross_point[i] = 0.;
      ires = 0;
      }
    else
      {
      t = -b/a;
      for(i=0;i<3;i++)
        cross_point[i] = t*dir[i] + r0[i];
      ires = 1;
      *dist = t;
      if(t < 0.)
        { *dist = -t; ires = -1;}
      }
        }
  return(ires);
}//source for int get_vcrpl(float *, float *, float *, float *, float *)


//*************************Photon Energy Correction****************//


Float_t h10::GetCorrFactStep1(Float_t E){

//return  1.129-0.05793/E-1.0773e-12/(E*E);   //corretion for C, Pb 
return  1.116-0.09213/E+0.01007/(E*E);         //correction for Fe

}
 

