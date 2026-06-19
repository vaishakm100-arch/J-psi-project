// pp collisions at 13 TeV

#include "Pythia8/Pythia.h"
// You need to include this to get access to the HIInfo object for
// HeavyIons.
#include "Pythia8/HeavyIons.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TRandom.h"
#include "TTree.h"
#include "TVirtualPad.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
// #include "Pythia8Plugins/ColourReconnectionHooks.h"

using namespace Pythia8;
using namespace std;

int main(int argc, char *argv[])
{

    TRandom *rnd = new TRandom();
    rnd->SetSeed();
    int seed = rnd->Uniform(2, 900000000);

    Pythia pythia;
    pythia.readString("Random:setSeed=on");
    pythia.readString(Form("Random:seed=%d", seed));

    cout << "Seed: " << seed << endl;

    // pythia.readString("PDF:lepton = off");

    pythia.readString("Beams:idA = 1000822080");
    pythia.readString("Beams:idB = 1000822080");
    pythia.readString("Beams:eCM = 5360.0");
    pythia.readString("HeavyIon:mode=1");
    pythia.readString("PhaseSpace:pTHatMinDiverge = 0.5");
    //pythia.readString("Tune:pp = 5");
    pythia.readString("Tune:pp = 5");

    pythia.readString("HardQCD:all=on");
    //pythia.readString("HardQCD:hardccbar=on"); 
    //pythia.readString("SoftQCD:all=on");
    pythia.readString("Charmonium:all=on");
    pythia.readString("Bottomonium:all=on");
    pythia.readString("PartonLevel:MPI = on");
    //pythia.readString("MultipartonInteractions:bProfile = 3");
    //pythia.readString("MultipartonInteractions:processLevel = 3");
    pythia.readString("ColourReconnection:reconnect=on");
    pythia.readString("ColourReconnection:mode = 2");
    pythia.readString("HadronLevel:decay = on");
    //pythia.readString("-313:mayDecay = true");
    //pythia.particleData.readString("421:0:products = -321 211");
    //pythia.readString("421:onMode = off");    
    //pythia.readString("421:onIfAny = -321 211");
    //pythia.readString();
    //Initialization.
  
    pythia.init();

    TFile *file = new TFile(Form("/home/vaishak/pythia8317/Pythia Data/Pb_Jpsitagged%d.root", time(0)), "recreate");
    
    //ofstream logfile;
   
    //logfile.open(Form("simulation_log_%d.txt", time(0)));
    //logfile << "Simulation started" << endl;
    //logfile << "Total events: " << nevents << endl;
    //logfile << "Seed: " << seed << endl;

    const int n = 10000;
    int ntrack = 0, tNch = 0, PID[n] = {0}, dPID[n] = {0}, charge[n] = {0}, Nchmid = 0, tNchV0 = 0;
    int tag[n] = {0}; //tagging 0 as default
    double px[n] = {0}, py[n] = {0}, pz[n] = {0}, pT[n] = {0}, mass[n] = {0};
    double eta[n] = {0}, rap[n] = {0}, phi[n] = {0};
    double pTlead = 0;
    double spherocitymid(0);
    double spherocityfwd(0);
    int AccEvtNch_v0 = 0, tnMPI = 0;
    TTree *theTree = new TTree("t", "Pythia8Tree");
    cout << "Setting the branchs" << endl;
    theTree->Branch("ntrack", &ntrack,"ntrack/I"); // ntrack is the size of the tree (fill only
                                 // final chargedparticles)
    theTree->Branch("Nchmid", &Nchmid, "Nchmid/I");
    theTree->Branch("tNchV0", &tNchV0, "tNchV0/I");
    theTree->Branch("tnMPI", &tnMPI, "tnMPI/I");
    theTree->Branch("PID", PID, "PID[ntrack]/I");
    theTree->Branch("dPID", dPID, "dPID[ntrack]/I");
    theTree->Branch("tag",tag, "tag[ntrack]/I"); //tagging
    theTree->Branch("charge", charge, "charge[ntrack]/I");
    theTree->Branch("px", px, "px[ntrack]/D");
    theTree->Branch("py", py, "py[ntrack]/D");
    theTree->Branch("pz", pz, "pz[ntrack]/D");
    theTree->Branch("mass", mass, "mass[ntrack]/D");
    theTree->Branch("eta", eta, "eta[ntrack]/D");
    theTree->Branch("pT", pT, "pT[ntrack]/D");
    theTree->Branch("pTlead", &pTlead, "pTlead/I");
    theTree->Branch("rap", rap, "rap[ntrack]/D");
    theTree->Branch("phi", phi, "phi[ntrack]/D");

    Double_t pixel = 100.0;
    Int_t nchcut = 0;
    Double_t dpy, dpx, spy, spx;
    Int_t j;
    Double_t crossp = 0;
    Double_t count = 0;
    int trks = 0;

    int bHadrons[]={511,521,531,541,5122,5132,5232,5332};  //B hadron array
    int nB=8;

    int cHadrons[]={411,421,431,441,4122,4132,4232,4332};  //c Hadron array
    int nC=8;

    //Helper lambda to check if ID belongs to a froup
    auto isInList=[](int id,int* list, int size){
        for (int k=0;k<size;k++)
           if(abs(id)==list[k]) return true;
        return false;
    };



// Simulate the events
     const int nevents = 10000;
    for (int iEvent = 0; iEvent < nevents; ++iEvent)
    {
        if (!pythia.next())
            continue;
        int cntnch = 0;
        int Nparticles = 0;
        
        // nmuon=0;
        AccEvtNch_v0 = 0;
        
        if (iEvent % 1000 == 0) {
            cout << "\n--- Progress: " << iEvent << "/" << nevents
                 << " (" << (100.0*iEvent/nevents) << "%) ---\n" << endl;
            //logfile << "--- Progress: " << iEvent << "/" << nevents
                    //<< " (" << (100.0*iEvent/nevents) << "%) ---" << endl;
            //logfile.flush();
        }
        if (iEvent % 10000 == 0 && iEvent > 0) {
            file->cd();
            theTree->AutoSave("SaveSelf");
            //logfile << "Checkpoint saved at event " << iEvent << endl;
            //logfile.flush();
        }
        for (Int_t i = 0; i < pythia.event.size(); i++)
        {
            if((pythia.event[i].id()==13 || pythia.event[i].id()==-13)){

                //Getting the parent of this muon
                int motherIndex = pythia.event[i].mother1();
                int motherID = pythia.event[motherIndex].id();

                //Only store if mother is J/psi
                
                    PID[Nparticles] = pythia.event[i].id();
                    dPID[Nparticles] = motherID; 
                     //cout<<"PID = "<<PID[Nparticles]<<endl;
                charge[Nparticles] = pythia.event[i].charge();
                px[Nparticles] = pythia.event[i].px();
                py[Nparticles] = pythia.event[i].py();
                pz[Nparticles] = pythia.event[i].pz();
                mass[Nparticles] = pythia.event[i].m();
                eta[Nparticles] = pythia.event[i].eta();
                pT[Nparticles] = pythia.event[i].pT();                
                rap[Nparticles] = pythia.event[i].y();
                phi[Nparticles] = pythia.event[i].phi();
                
                //Find J/psi mother
                int jpsiMotherID=0;
                int jpsiMotherIndex=0;
                int currentIndex=i;

                //Find the J/psi in the chain first
                for(int gen=0; gen<5; gen++){
                    int momIndex=pythia.event[currentIndex].mother1();
                    if(momIndex <= 0) break;
                    int momID=abs(pythia.event[momIndex].id());

                    if(momID==443|| momID==100443){
                        currentIndex=momIndex;
                        continue;
                    }
                    jpsiMotherID=momID;
                    jpsiMotherIndex=momIndex;
                    cout << "Muon mother chain: J/psi <- parent(" << jpsiMotherID<<")" <<endl;
                    break;

                    }
                        
                    }
                    
                }
                //For checking whether J/psi parent is bottom -directly non-prompt
                if(isInList(jpsiMotherID, bHadrons, nB)){
                    tag[Nparticles]=2;
                    cout << "Tag=2 (non-prompt,direct B hadron)" <<endl;
                } 
                //if parent is charm,tracing back further
                else if(isInList(jpsiMotherID, cHadrons, nC)){
                    cout<< " Charm parent found,tracing back..." <<endl;
                    currentIndex=jpsiMotherIndex;
                    bool foundB=false;

                    for(int gen=0; gen<5; gen++){
                        int momIndex=pythia.event[currentIndex].mother1();
                        if(momIndex <=0) break;
                        int momID=abs(pythia.event[momIndex].id());

                        cout <<" gen" << gen+1 << "parent ID:" <<momID<<endl;

                        if(isInList(momID, bHadrons, nB)){
                            foundB = true;
                            cout << "TAG=2 (non-prompt, B found at gen" <<gen+1 <<")" << endl;
                            break;
                        }
                        currentIndex=momIndex;
                
                    }
                    if(foundB)
                        tag[Nparticles]=2;
                    else
                        tag[Nparticles]=1;            
                }
                //direct J/psi, prompt
                else{
                    tag[Nparticles]=1;
                    cout << " TAG=1 (prompt)" <<endl;
                }

                Nparticles++;
            }
    
        
                
                
               
            
        }
        
        tnMPI = pythia.info.nMPI();

        if (Nparticles > 0) 
        {
            ntrack = Nparticles;
            theTree->Fill();
            cout<<"Event: " << iEvent+1 << "/" << nevents
                << " | Muons found: " << Nparticles
                << " | pT = " << pT[0] << " GeV"
                << " | eta = " << eta[0]
                << " | rap = " << rap[0]
                << " | nMPI = " << tnMPI
                <<endl;
           // logfile << "Event: " << iEvent+1 << "/" << nevents
                    //<< " | Muons found: " << Nparticles
                   // << " | pT = " << pT[0] << "GeV"
                   // << " | eta = " << eta[0]
                   // << " | rap = " << rap[0]
                    //<< " | nMPI = " << tnMPI
                    //<<endl;
        }
    }
    pythia.stat();
    file->cd();
    theTree->Write();
    file->Close();
    
    //logfile << "simulation completed" << endl;
    //logfile.close();
    return 0;
}
