#include "Riostream.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TClonesArray.h"
#include "TClasTool.h"
#include "TIdentificator.h"
#include "TMath.h"
#include "TString.h"
#include "TMath.h"
#include "massConst.h"
using namespace std;

#define MAX_ELECTRONS 1
#define MAX_PIPLUS 1
#define MAX_PIMINUS 1
#define MAX_PHOTONS 2

#define PDG_PHOTON 22
#define PDG_POSITRON -11
#define PDG_ELECTRON 11
#define PDG_PIPLUS 211
#define PDG_PIMINUS -211
#define PDG_KPLUS 321
#define PDG_KMINUS -321
#define PDG_NEUTRON 2112
#define PDG_PROTON 2212

#define GEANT3_PHOTON 1
#define GEANT3_POSITRON 2
#define GEANT3_ELECTRON 3
#define GEANT3_PIPLUS 8
#define GEANT3_PIMINUS 9
#define GEANT3_KPLUS 11
#define GEANT3_KMINUS 12
#define GEANT3_NEUTRON 13
#define GEANT3_PROTON 14

#define CUT_Q2 1.0
#define CUT_W 2.0
#define CUT_NU 0.85

#define EBEAM 5.015  // e- beam energy in GeV

//declarations of functions
void PrintAnalysisTime(float tStart, float tStop);
void PrintUsage(char *processName);
int GetPID(string partName, int kind);

typedef struct{
    Float_t EvtNum, ElecVertTarg, Q2, Nu, Xb, W;
    Float_t Xcorr, Ycorr, Zcorr;
    Int_t nElec, nPip, nPim, nGam, nProton, nNeutron, nKp, nKm, nPositron;
} KINEVAR;

int main(int argc, char **argv)
{
    extern char *optarg;
    int c;
    extern int optind;
    
    int i, j, k;
    int nRows, kind, tempPid, savePid;
    int candCtr = 0;
    int ctr_nElec = 0;
    int dEvents = 1000; // increment of events for processing print statement
    int MaxEvents = 0; // max. number of events to process
    int nfiles = 0; // number of processed files
    
    TString catPid;
    
    bool bBatchMode = false;    // quiet mode
    bool simul_key = false;  // simulation flag (true = simulation, false = data)
    bool mflag = true; // cut flag for GetCategorization(k,tt,mflag)
    int cat_key = 0; // PID categorization 0 = EVNT (default), 1 = Full
    int tgt_key = 1;  // intitialize target flag 1 = Carbon (default), 2 = Iron, 3 = Lead
    string target; // solid target name

    char *inFile;
    string outFile = "PipPimPi0.root";
    
    bool partFound = false;
    bool topology = false;

    TVector3 *vert;
    TVector3 *ECxyz = new TVector3(0.0,0.0,0.0);
    TVector3 *ECuvw;
    
    float timeStart = clock(); // start time
    
    TClasTool *input = new TClasTool();
    input->InitDSTReader("ROOTDSTR");
  
    TIdentificator *t = new TIdentificator(input);
    
    for (i = 0; i < argc; ++i) cerr << argv[i] << " "; cerr << endl;
    while ((c = getopt(argc,argv, "o:M:D:c:t:Sih")) != -1 ) {
        switch (c) {
            case 'o': outFile = optarg; break;
            case 'M': MaxEvents = atoi(optarg); break;
            case 'D': dEvents = atoi(optarg); break;
            case 'c': cat_key = atoi(optarg); break;
            case 't': tgt_key = atoi(optarg); break;
            case 'S': simul_key = true; break;
            case 'i': bBatchMode = true; break;
            case 'h':
                PrintUsage(argv[0]);
                exit(0);
                break;
                
            default:
                cerr << "Unrecognized argument: " << optarg << endl;
                PrintUsage(argv[0]);
                exit(0);
                break;
        }
    }
    
    TFile *output; // ROOT output file
    
    // check target selection
    switch(tgt_key){
        case 1: target = "C"; break;
        case 2: target = "Fe"; break;
        case 3: target = "Pb"; break;
        default: cout<<"Unknown target "<<target<<endl; exit(0); break;
    }
    cout<<"Analyzing " << target << " target data"<<endl;
  
    TClonesArray *Sector = new TClonesArray("Int_t");
    TClonesArray *Charge = new TClonesArray("Double_t");
    TClonesArray *Pid = new TClonesArray("Double_t");
    TClonesArray *Beta = new TClonesArray("Double_t");
    TClonesArray *Px = new TClonesArray("Double_t");
    TClonesArray *Py = new TClonesArray("Double_t");
    TClonesArray *Pz = new TClonesArray("Double_t");
    TClonesArray *Mom = new TClonesArray("Double_t");
    TClonesArray *Mass2 = new TClonesArray("Double_t");
    TClonesArray *X = new TClonesArray("Double_t");
    TClonesArray *Y = new TClonesArray("Double_t");
    TClonesArray *Z = new TClonesArray("Double_t");
    TClonesArray *ECx = new TClonesArray("Double_t");
    TClonesArray *ECy = new TClonesArray("Double_t");
    TClonesArray *ECz = new TClonesArray("Double_t");
    TClonesArray *ECu = new TClonesArray("Double_t");
    TClonesArray *ECv = new TClonesArray("Double_t");
    TClonesArray *ECw = new TClonesArray("Double_t");
    TClonesArray *ECtot = new TClonesArray("Double_t");
    TClonesArray *ECin = new TClonesArray("Double_t");
    TClonesArray *ECout = new TClonesArray("Double_t");
    TClonesArray *ECtime = new TClonesArray("Double_t");
    TClonesArray *ECpath = new TClonesArray("Double_t");
    TClonesArray *EChit_M2 = new TClonesArray("Double_t");
    TClonesArray *EChit_M3 = new TClonesArray("Double_t");
    TClonesArray *EChit_M4 = new TClonesArray("Double_t");
    TClonesArray *Chi2EC = new TClonesArray("Double_t");
    TClonesArray *SCpath = new TClonesArray("Double_t");
    TClonesArray *SCtime = new TClonesArray("Double_t");
    TClonesArray *CCnphe = new TClonesArray("Double_t");
    TClonesArray *T = new TClonesArray("Double_t");
    TClonesArray *Xf = new TClonesArray("Double_t");
    TClonesArray *Mx2 = new TClonesArray("Double_t");
    TClonesArray *Pt = new TClonesArray("Double_t");
    TClonesArray *Zh = new TClonesArray("Double_t");
    TClonesArray *ThetaPQ = new TClonesArray("Double_t");
    TClonesArray *PhiPQ = new TClonesArray("Double_t");
    TClonesArray *TimeCorr4 = new TClonesArray("Double_t");
    
    string kineList = "EvtNum/F:ElecVertTarg/F:Q2/F:Nu/F:Xb/F:W:Xcorr/F:Ycorr/F:Zcorr/F:nElec/I:nPip/I:nPim/I:nGam/I:nProton/I:nNeutron/I:nKp/I:nKm/I:nPositron/I";
    KINEVAR myKine;
    
    TTree *dataTree = new TTree("Data","Experimental Data Tree");
    dataTree->Branch("Kinematics",&myKine,kineList.c_str());

    Int_t iTCA = 15;
    dataTree->Branch("Sector",&Sector,iTCA);
    dataTree->Branch("Charge",&Charge,iTCA);
    dataTree->Branch("Pid",&Pid,iTCA);
    dataTree->Branch("Beta",&Beta,iTCA);
    dataTree->Branch("Px",&Px,iTCA);
    dataTree->Branch("Py",&Py,iTCA);
    dataTree->Branch("Pz",&Pz,iTCA);
    dataTree->Branch("Mom",&Mom,iTCA);
    dataTree->Branch("Mass2",&Mass2,iTCA);
    dataTree->Branch("X",&X,iTCA);
    dataTree->Branch("Y",&Y,iTCA);
    dataTree->Branch("Z",&Z,iTCA);
    dataTree->Branch("ECx",&ECx,iTCA);
    dataTree->Branch("ECy",&ECy,iTCA);
    dataTree->Branch("ECz",&ECz,iTCA);
    dataTree->Branch("ECu",&ECu,iTCA);
    dataTree->Branch("ECv",&ECv,iTCA);
    dataTree->Branch("ECw",&ECw,iTCA);
    dataTree->Branch("ECtot",&ECtot,iTCA);
    dataTree->Branch("ECin",&ECin,iTCA);
    dataTree->Branch("ECout",&ECout,iTCA);
    dataTree->Branch("ECtime",&ECtime,iTCA);
    dataTree->Branch("ECpath",&ECpath,iTCA);
    dataTree->Branch("EChit_M2",&EChit_M2,iTCA);
    dataTree->Branch("EChit_M3",&EChit_M3),
    dataTree->Branch("EChit_M4",&EChit_M4,iTCA);
    dataTree->Branch("Chi2EC",&Chi2EC,iTCA);
    dataTree->Branch("SCpath",&SCpath,iTCA);
    dataTree->Branch("SCtime",&SCtime,iTCA);
    dataTree->Branch("CCnphe",&CCnphe,iTCA);
    dataTree->Branch("T",&T,iTCA);
    dataTree->Branch("Xf",&Xf,iTCA);
    dataTree->Branch("Mx2",&Mx2,iTCA);
    dataTree->Branch("Pt",&Pt,iTCA);
    dataTree->Branch("Zh",&Zh,iTCA);
    dataTree->Branch("ThetaPQ",&ThetaPQ,iTCA);
    dataTree->Branch("PhiPQ",&PhiPQ,iTCA);
    dataTree->Branch("TimeCorr4",&TimeCorr4,iTCA);

    output = new TFile(outFile.c_str(), "RECREATE", "Experimental Data");
    
    for (i = optind; i < argc; ++i) {
        inFile = argv[i]; // process all arguments on command line.
        if (*inFile != '-') { // we have a file to process
            cout << "Analyzing file " << inFile << endl; // let user know which file is being processed
            
            input->Add(inFile); // read file into ClasTool object
            
            nfiles++; // increment file counter
        }
    }

    Long_t nEntries = (Long_t) input->GetEntries(); // get total number of events
    
    cout<<"Analyzing "<<nEntries<<" from "<<nfiles<< " files."<<endl; // print out stats
  
    input->Next();
  
    k = 0; // event counter
    
    if(MaxEvents == 0) MaxEvents = nEntries; // if user does not set max. number of events, set to nEntries
    
    while (k < MaxEvents) {
    	if (!bBatchMode && ((k % dEvents) == 0)){
    		cerr << k << "\r";
    	}

        if(simul_key){
            kind = 1;
            nRows = input->GetNRows("GSIM");
        }else{
            kind = 0;
            nRows = input->GetNRows("EVNT");
        }
        
        memset(&myKine,0,sizeof(myKine)); // init kinematics struct to zeros
        Sector.Clear();
        Charge.Clear();
        Pid.Clear();
        Beta.Clear();
        Px.Clear();
        Py.Clear();
        Pz.Clear();
        Mom.Clear();
        Mass2.Clear();
        X.Clear();
        Y.Clear();
        Z.Clear();
        ECx.Clear();
        ECy.Clear();
        ECz.Clear();
        ECu.Clear();
        ECv.Clear();
        ECw.Clear();
        ECtot.Clear();
        ECin.Clear();
        ECout.Clear();
        ECtime.Clear();
        ECpath.Clear();
        EChit_M2.Clear();
        EChit_M3.Clear();
        EChit_M4.Clear();
        Chi2EC.Clear();
        SCpath.Clear();
        SCtime.Clear();
        CCnphe.Clear();
        T.Clear();
        Xf.Clear();
        Mx2.Clear();
        Pt.Clear();
        Zh.Clear();
        ThetaPQ.Clear();
        PhiPQ.Clear();
        TimeCorr4.Clear();
        
        if(nRows>0){
	    	topology = false; // init. the event topology cut
            ip = 0; // found particle index
            
            for (j = 0; j < nRows; j++) {

                partFound = false; // init the found particle flag to false
                
                // select the PID selection scheme
                if(simul_key){
                    catPid = t -> GetCategorizationGSIM(j);
                }else{
                    switch(cat_key){
                        case 0: catPid = t -> GetCategorizationEVNT(j); break;
                        case 1: catPid = t -> GetCategorization(j,target.c_str(),mflag); break;
                        default: cout<<"Incorrect PID categorization.  Try again."<<endl; exit(0); break;
                    }
                }
                tempPid = t -> Id(j,kind);
                if(tempPid == GetPID("Proton",kind)) myKine.nProton++;
                if(tempPid == GetPID("Neutron",kind)) myKine.nNeutron++;
                if(tempPid == GetPID("KPlus",kind)) myKine.nKp++;
                if(tempPid == GetPID("KMinus",kind)) myKine.nKm++;
                if(tempPid == GetPID("Positron",kind)) myKine.nPositron++;
                
//                cout<<"Particle "<< tempPid <<"\t"<<catPid<<endl;
//                if(tempPid == GetPID("Electron",kind)){
                if(catPid.EqualTo("electron")){
                    myKine.nElec++;
                    ctr_nElec++;
                    partFound = true; // init the found particle flag to true
                    savePid = GetPID("Electron",kind); // set the correct particle id
                }
                if(catPid.EqualTo("high energy pion +") || catPid.EqualTo("low energy pion +") || catPid.EqualTo("pi+")){
//                if(tempPid == GetPID("PiPlus",kind)){
                    myKine.nPip++;
                    partFound = true; // init the found particle flag to true
                    savePid = GetPID("PiPlus",kind); // set the correct particle id
                }
                if(catPid.EqualTo("pi-")){
//                if(tempPid == GetPID("PiMinus",kind)){
                    myKine.nPim++;
                    partFound = true; // init the found particle flag to true
                    savePid = GetPID("PiMinus",kind); // set the correct particle id
                }
                if(catPid.EqualTo("gamma")){
//                if(tempPid == GetPID("Photon",kind)){
                    myKine.nGam++;
                    partFound = true; // init the found particle flag to true
                    savePid = GetPID("Photon",kind); // set the correct particle id
                }

        		if (partFound) {
                    new ((*Sector)[ip]) t->Sector(j,kind);
                    new ((*Charge)[ip]) t->Charge(j,kind);
                    new ((*Beta)[ip]) t->Betta(j,kind);
//                    new ((*Pid)[ip]) t->Id(j,kind);
                    new ((*Pid)[ip]) savePid;
                    new ((*Mom)[ip]) t->Momentum(j,kind);
                    new ((*Px)[ip]) t->Px(j,kind);
                    new ((*Py)[ip]) t->Py(j,kind);
                    new ((*Pz)[ip]) t->Pz(j,kind);
                    new ((*X)[ip]) t->X(j,kind);
                    new ((*Y)[ip]) t->Y(j,kind);
                    new ((*Z)[ip]) t->Z(j,kind);
                    new ((*Mass2)[ip]) t->Mass2(j,kind);
                    new ((*ThetaPQ)[ip]) t -> ThetaPQ(j,kind);
                    new ((*PhiPQ)[ip]) t -> PhiPQ(j,kind);
                    new ((*Zh)[ip]) t -> Zh(j,kind);
                    new ((*Pt)[ip]) TMath::Sqrt(t -> Pt2(j,kind));
                    new ((*Mx2)[ip]) t -> Mx2(j,kind);
                    new ((*Xf)[ip]) t -> Xf(j,kind);
                    new ((*T)[ip]) t -> T(j,kind);
                         
                    if(simul_key == false){
                        new ((*ECtot)[ip]) TMath::Max(t->Etot(j),t->Ein(j)+t->Eout(j));
                        new ((*ECin)[ip]) t->Ein(j);
                        new ((*ECout)[ip]) t->Eout(j);
                        new ((*ECx)[ip]) t->XEC(j);
                        new ((*ECy)[ip]) t->YEC(j);
                        new ((*ECz)[ip]) t->ZEC(j);
                        ECxyz->SetXYZ(t->XEC(j),t->YEC(j),t->ZEC(j));
                        ECuvw = t->XYZToUVW(ECxyz);
                        new ((*ECu)[ip]) ECuvw->X();
                        new ((*ECv)[ip]) ECuvw->Y();
                        new ((*ECw)[ip]) ECuvw->Z();
                        new ((*ECtime)[ip]) t->TimeEC(j);
                        new ((*ECpath)[ip]) t->PathEC(j);
                        new ((*EChit_M2)[ip]) t->EChit_Moment2(j);
                        new ((*EChit_M3)[ip]) t->EChit_Moment3(j);
                        new ((*EChit_M4)[ip]) t->EChit_Moment4(j);
                        new ((*Chi2EC)[ip]) t->Chi2EC(j);
                        new ((*SCtime)[ip]) t->TimeSC(j);
                        new ((*SCpath)[ip]) t->PathSC(j);
                        new ((*CCnphe)[ip]) t->Nphe(j);

                        if(t->Id(j,kind) == GetPID("Electron",kind)) new ((*TimeCorr4)[ip]) t -> TimeCorr4(0.000511,j);
                        if(t->Id(j,kind) == GetPID("PiPlus",kind)) new ((*TimeCorr4)[ip]) t -> TimeCorr4(kMassPi_plus,j);
                        if(t->Id(j,kind) == GetPID("PiMinus",kind)) new ((*TimeCorr4)[ip]) t -> TimeCorr4(kMassPi_min,j);
                        if(t->Id(j,kind) == GetPID("Photon",kind)) new ((*TimeCorr4)[ip]) t -> TimeCorr4(0.0,j);
                    }
                    ip++; // increment the found particle index
                }
            }
            
            topology = (myKine.nElec>=MAX_ELECTRONS && myKine.nPip>=MAX_PIPLUS && myKine.nPim>=MAX_PIMINUS && myKine.nGam>=MAX_PHOTONS); // check event topology
                
            if(topology && t->Q2(kind) > CUT_Q2 && t->W(kind) > CUT_W && t->Nu(kind)/EBEAM < CUT_NU) {
                candCtr++;
                myKine.EvtNum = t -> NEvent();
                myKine.ElecVertTarg = t -> ElecVertTarg(kind);
                myKine.Q2 = t -> Q2(kind);
                myKine.Nu = t -> Nu(kind);
                myKine.Xb = t -> Xb(kind);
                myKine.W = t -> W(kind);
                
                if(simul_key){
                    myKine.Xcorr = t->X(0, kind);
                    myKine.Ycorr = t->Y(0, kind);
                    myKine.Zcorr = t->Z(0, kind);
                }else{
                    vert = t->GetCorrectedVert();
                    myKine.Xcorr = vert->X();
                    myKine.Ycorr = vert->Y();
                    myKine.Zcorr = vert->Z();
                }
                    
                dataTree->Fill();
            }
    	} 
    	k++; // increment event counter
        input->Next();
    }
//    dataTree->Print();
//    dataTree->Scan("Kinematics.EvtNum:Electron.Pid:PiPlus.Pid:PiMinus.Pid:Photon1.Pid:Photon2.Pid:Photon2.Beta");

    dataTree->Write();
    cout<<"Candidate data events = "<<candCtr<<endl;
    cout<<"#(e-) = "<<ctr_nElec<<endl;
                   
    output->Write();
    output->Close();
    
    float timeStop = clock();
    PrintAnalysisTime(timeStart,timeStop);

    // delete the TClonesArray objects
    ~Sector;
    ~Charge;
    ~Pid;
    ~Beta;
    ~Px;
    ~Py;
    ~Pz;
    ~Mom;
    ~Mass2;
    ~X;
    ~Y;
    ~Z;
    ~ECx;
    ~ECy;
    ~ECz;
    ~ECu;
    ~ECv;
    ~ECw;
    ~ECtot;
    ~ECin;
    ~ECout;
    ~ECtime;
    ~ECpath;
    ~EChit_M2;
    ~EChit_M3;
    ~EChit_M4;
    ~Chi2EC;
    ~SCpath;
    ~SCtime;
    ~CCnphe;
    ~T;
    ~Xf;
    ~Mx2;
    ~Pt;
    ~Zh;
    ~ThetaPQ;
    ~PhiPQ;
    ~TimeCorr4;
    
    return 0;
}

void PrintUsage(char *processName)
{
    cerr << processName << " <options> <filename>\n";
    cerr << "\toptions are:\n";
    cerr << "\t-o<filename>\tROOT output file (def. = PipPimPi0.root).\n";
    cerr << "\t-M#\t\tprocess maximum # of events.\n";
    cerr << "\t-D#\t\tinform user when # of events have been processed (def. = 1000).\n";
    cerr << "\t-c#\t\tType in categoization # scheme of 0=EVNT or 1=Full (def. = 0).\n";
    cerr << "\t-t#\t\tTarget # of 1=C, 2=Fe, or 3=Pb (def. = 1).\n";
    cerr << "\t-S\t\tAnalyze simulation.\n";
    cerr << "\t-i\t\tquiet mode (no counter).\n";
    cerr << "\t-h\t\tprint the above" << endl;
}


void PrintAnalysisTime(float tStart, float tStop){
    //time to complete function
    float minutes = 0;
    float seconds = 0;
    minutes = (tStop - tStart)/1000000;
    minutes = (minutes)/60;
    seconds = fmod(minutes,1);
    minutes = minutes-seconds;
    seconds = seconds*60;
    
    if (minutes==0){
        cout<<endl<<"Completed in "<<seconds<<" seconds."<<endl<<endl;
    }
    else{
        cout<<endl<<"Completed in "<<minutes<<" minutes and "<<seconds<<" seconds."<<endl<<endl;
    }
}

int GetPID(string partName, int kind){

    int ret = 0;
    
    if(kind==0){
        if(partName.compare("Electron")==0){
            ret = PDG_ELECTRON;
        }else if(partName.compare("Positron")==0){
            ret = PDG_POSITRON;
        }else if(partName.compare("Photon")==0){
            ret = PDG_PHOTON;
        }else if(partName.compare("PiPlus")==0){
            ret = PDG_PIPLUS;
        }else if(partName.compare("PiMinus")==0){
            ret = PDG_PIMINUS;
        }else if(partName.compare("KPlus")==0){
            ret = PDG_KPLUS;
        }else if(partName.compare("KMinus")==0){
            ret = PDG_KMINUS;
        }else if(partName.compare("Neutron")==0){
            ret = PDG_NEUTRON;
        }else if(partName.compare("Proton")==0){
            ret = PDG_PROTON;
        }else{
            cerr<<"GetPid(): Unknown PDG particle "<<partName.c_str()<<endl; exit(0);
        }
    }else if(kind==1){
        if(partName.compare("Electron")==0){
            ret = GEANT3_ELECTRON;
        }else if(partName.compare("Positron")==0){
            ret = GEANT3_POSITRON;
        }else if(partName.compare("Photon")==0){
            ret = GEANT3_PHOTON;
        }else if(partName.compare("PiPlus")==0){
            ret = GEANT3_PIPLUS;
        }else if(partName.compare("PiMinus")==0){
            ret = GEANT3_PIMINUS;
        }else if(partName.compare("KPlus")==0){
            ret = GEANT3_KPLUS;
        }else if(partName.compare("KMinus")==0){
            ret = GEANT3_KMINUS;
        }else if(partName.compare("Neutron")==0){
            ret = GEANT3_NEUTRON;
        }else if(partName.compare("Proton")==0){
            ret = GEANT3_PROTON;
        }else{
            cerr<<"GetPid(): Unknown GEANT3 particle "<<partName.c_str()<<endl; exit(0);
        }
    }else{
        cerr<<"GetPID: Unknown analysis channel "<<kind<<endl;
    }
    return ret;
}
