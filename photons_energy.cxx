#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TF1.h>
#include "h10.h"

using namespace std;
int main(int argc, char *argv[])
{
TChain *chain;
chain = new TChain("h10","Analysis chain");
chain->Add("root_41146_01.pass2.root");
/*
gSystem->Exec("rm -rf list.dat;touch list.dat;");
gSystem->Exec("ls -lh *.root | awk '{print $9}' >> list.dat;");
ifstream list("list.dat");
TString filepath;
do{
 list >> filepath;
 if(filepath.Length()>0){
  TFile* FirstFile = new TFile(Form("%s",filepath.Data()));
  TList *FileContent;
  FileContent = FirstFile->GetListOfKeys();
  if(FileContent->Contains("h10")){
   cout << "Adding >" << filepath.Data() << "<" << endl;
   chain->Add(filepath.Data());
  }
 }
}while(!list.eof());
list.close();
*/
int Ntoproc;
Ntoproc = chain->GetEntries();
cout << "Number of events in chain : " << Ntoproc;
if(argc>1)Ntoproc=atoi(argv[1]);
cout << " Will process : " << Ntoproc << endl;
h10 *selec = new h10(0,Ntoproc);
chain->Process(selec,"",Ntoproc);
cout << "Back in main" << endl;

return 0;
}

