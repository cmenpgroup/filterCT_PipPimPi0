#include "Riostream.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TClasTool.h"
#include "TIdentificator.h"
#include "TMath.h"
#include "massConst.h"
using namespace std;

#define MAX_ELECTRONS 1
#define MAX_PIPLUS 1
#define MAX_PIMINUS 1
#define MAX_PHOTONS 2

#define PDG_PHOTON 22
#define PDG_ELECTRON 11
#define PDG_PIPLUS 211
#define PDG_PIMINUS -211

#define GEANT3_PHOTON 1
#define GEANT3_ELECTRON 3
#define GEANT3_PIPLUS 8
#define GEANT3_PIMINUS 9

//declarations of functions
void PrintAnalysisTime(float tStart, float tStop);
void PrintUsage(char *processName);
int GetPID(string partName, int kind);

typedef struct{
    Float_t EvtNum, ElecVertTarg, Q2, Nu, Xb, W;
    Float_t Xcorr, Ycorr, Zcorr;
    Int_t nElec, nPip, nPim, nGam;
} KINEVAR;

typedef struct{
    int Sector;
    float Charge, Pid, Beta;
    float Px, Py, Pz, Mom, Mass2;
    float X, Y, Z;
    float ECx, ECy, ECz, ECu, ECv, ECw, ECtot, ECin, ECout, ECtime, ECpath;
    float EChit_M2, EChit_M3, EChit_M4, Chi2EC;
    float SCpath, SCtime;
    float CCnphe; 
    float T, Xf, Mx2, Pt, Zh, ThetaPQ, PhiPQ, TimeCorr4;
} PARTVAR;

int main(int argc, char **argv)
{
    extern char *optarg;
    int c;
    extern int optind;
    
    int i, j, k;
    int nRows, tempPid;
    int photonCtr;
    int candCtr = 0;
    int dEvents = 1000; // increment of events for processing print statement
    int MaxEvents = 0; // max. number of events to process
    int nfiles = 0; // number of processed files
    int kind = 0; // initialize bank index for EVNT data
    
    bool bBatchMode = false;    // quiet mode

    char *inFile;
    string outFile = "PipPimPi0.root";
    
    bool topology = false;
    vector<int> partIndex;
    vector<int> partID;
    
    TVector3 *vert;
    TVector3 *ECxyz = new TVector3(0.0,0.0,0.0);
    TVector3 *ECuvw;
    
    float timeStart = clock(); // start time
    
    TClasTool *input = new TClasTool();
    input->InitDSTReader("ROOTDSTR");
  
    TIdentificator *t = new TIdentificator(input);
    
    for (i = 0; i < argc; ++i) cerr << argv[i] << " "; cerr << endl;
    while ((c = getopt(argc,argv, "o:M:D:Sih")) != -1 ) {
        switch (c) {
            case 'o': outFile = optarg; break;
            case 'M': MaxEvents = atoi(optarg); break;
            case 'D': dEvents = atoi(optarg); break;
            case 'S': kind = 1; break;
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

    string kineList = "EvtNum/F:ElecVertTarg/F:Q2/F:Nu/F:Xb/F:W:Xcorr/F:Ycorr/F:Zcorr/F:nElec/I:nPip/I:nPim/I:nGam/I";
 
    string partList = "Sector/I:Charge/F:Pid/F:Beta/F:Px/F:Py/F:Pz/F:Mom/F:Mass2/F:X/F:Y/F:Z/F:ECx/F:ECy/F:ECz/F:ECu/F:ECv/F:ECw/F:ECtot/F:ECin/F:ECout/F:ECtime/F:ECpath/F:EChit_M2/F:EChit_M3/F:EChit_M4/F:Chi2EC/F:SCpath/F:SCtime/F:CCnphe/F:T/F:Xf/F:Mx2/F:Pt/F:Zh/F:ThetaPQ/F:PhiPQ/F:TimeCorr4/F";    
    KINEVAR myKine;
    PARTVAR myPart;
    PARTVAR myElec;
    PARTVAR myPip;
    PARTVAR myPim;
    PARTVAR myPhoton1;
    PARTVAR myPhoton2;

    TTree *dataTree = new TTree("Data","Experimental Data Tree");
    dataTree->Branch("Kinematics",&myKine,kineList.c_str());
    dataTree->Branch("Electron",&myElec,partList.c_str());
    dataTree->Branch("PiPlus",&myPip,partList.c_str());
    dataTree->Branch("PiMinus",&myPim,partList.c_str());
    dataTree->Branch("Photon1",&myPhoton1,partList.c_str());
    dataTree->Branch("Photon2",&myPhoton2,partList.c_str());

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

        myKine.nElec = 0;
        myKine.nPip = 0;
        myKine.nPim = 0;
        myKine.nGam = 0;
        partIndex.clear();
        partID.clear();
        topology = false;
        
        cout<<"Event "<<k+1<<endl;

        if(kind==1){
            nRows = input->GetNRows("GSIM");
        }else{
            nRows = input->GetNRows("EVNT");
        }
        
        if(nRows>0){
            cout<<"Particle 0 "<< t->Id(0,kind) <<" "<<t->GetCategorizationParticle(0,kind)<<endl;
            if(t->GetCategorizationParticle(0, kind) == "electron"){
                myKine.nElec++;
                if(myKine.nElec>0 && myKine.nElec<=MAX_ELECTRONS){
                    partIndex.push_back(0);
                    partID.push_back(GetPID("Electron",kind));
                }
            }
            for (j = 1; j < nRows; j++) {
                tempPid = t -> Id(j,kind);
                cout<<"Particle "<< tempPid <<" "<<t->GetCategorizationParticle(j,kind)<<endl;

                if(tempPid == GetPID("PiPlus",kind)){
                    myKine.nPip++;
                    if(myKine.nPip>0 && myKine.nPip<=MAX_PIPLUS){
                        partIndex.push_back(j);
                        partID.push_back(GetPID("PiPlus",kind));
                    }
                }
                if(tempPid == GetPID("PiMinus",kind)){
                    myKine.nPim++;
                    if(myKine.nPim>0 && myKine.nPim<=MAX_PIMINUS){
                        partIndex.push_back(j);
                        partID.push_back(GetPID("PiMinus",kind));
                    }
                }
                
                if(t->GetCategorizationParticle(j, kind) == "gamma"){
                    myKine.nGam++;
                    if(myKine.nGam>0 && myKine.nGam<=MAX_PHOTONS){
                        partIndex.push_back(j);
                        partID.push_back(GetPID("Photon",kind));
                    }
                }
            }
	    	topology = (myKine.nElec>0 && myKine.nPip>0 && myKine.nPim>0 && myKine.nGam>=MAX_PHOTONS && partIndex.size()==partID.size()); // check event topology

	    	if(topology && t->Q2(kind) > 1. && t->W(kind) > 2. && t->Nu(kind)/5.015 < 0.85) {
                candCtr++;
                myKine.EvtNum = t -> NEvent();
                myKine.ElecVertTarg = t -> ElecVertTarg(kind);
                myKine.Q2 = t -> Q2(kind);
		     	myKine.Nu = t -> Nu(kind);
	       		myKine.Xb = t -> Xb(kind);
        		myKine.W = t -> W(kind);

                if(kind==1){
                    myKine.Xcorr = t->X(0, kind);
                    myKine.Ycorr = t->Y(0, kind);
                    myKine.Zcorr = t->Z(0, kind);
                }else{
                    vert = t->GetCorrectedVert();
                    myKine.Xcorr = vert->X();
                    myKine.Ycorr = vert->Y();
                    myKine.Zcorr = vert->Z();
                }

	    		photonCtr = 0;
        		while (!partIndex.empty()) {
		    		i = partIndex.back(); // retrieve EVNT/GSIM index for each particle
                    partIndex.pop_back(); // erase last entry in the list

                    myPart.Sector = t->Sector(kind);
                    myPart.Charge = t->Charge(kind);
                    myPart.Beta = t->Betta(kind);
                    myPart.Pid = partID.back(); // retrieve PID
                    partID.pop_back(); // erase last entry in the list
                    myPart.Mom = t->Momentum(i,kind);
                    myPart.Px = t->Px(i, kind);
                    myPart.Py = t->Py(i, kind);
                    myPart.Pz = t->Pz(i, kind);
                    myPart.X = t->X(i, kind);
                    myPart.Y = t->Y(i, kind);
                    myPart.Z = t->Z(i, kind);
                    myPart.Mass2 = t->Mass2(i, kind);

                    myPart.ThetaPQ = t -> ThetaPQ(i, kind);
                    myPart.PhiPQ = t -> PhiPQ(i, kind);
                    myPart.Zh = t -> Zh(i, kind);
                    myPart.Pt = TMath::Sqrt(t -> Pt2(i, kind));
                    myPart.Mx2 = t -> Mx2(i, kind);
                    myPart.Xf = t -> Xf(i, kind);
        			myPart.T = t -> T(i, kind);

                    // initialize detector info
                    myPart.ECtot = 0;
                    myPart.ECin = 0;
                    myPart.ECout = 0;
                    myPart.ECx = 0;
                    myPart.ECy = 0;
                    myPart.ECz = 0;
                    myPart.ECu = 0;
                    myPart.ECv = 0;
                    myPart.ECw = 0;
                    myPart.ECtime = 0;
                    myPart.ECpath = 0;
                    myPart.EChit_M2 = 0;
                    myPart.EChit_M3 = 0;
                    myPart.EChit_M4 = 0;
			        myPart.Chi2EC = 0;
                    myPart.SCtime = 0;
                    myPart.SCpath = 0;
                    myPart.CCnphe = 0;
                    myPart.TimeCorr4 = 0;
                    
                    if(kind == 0){
                        myPart.ECtot = TMath::Max(t->Etot(i),t->Ein(i)+t->Eout(i));
                        myPart.ECin = t->Ein(i);
                        myPart.ECout = t->Eout(i);
                        myPart.ECx = t->XEC(i);
                        myPart.ECy = t->YEC(i);
                        myPart.ECz = t->ZEC(i);
                        ECxyz->SetXYZ(t->XEC(i),t->YEC(i),t->ZEC(i));
                        ECuvw = t->XYZToUVW(ECxyz);
                        myPart.ECu = ECuvw->X();
                        myPart.ECv = ECuvw->Y();
                        myPart.ECw = ECuvw->Z();
 				        myPart.ECtime = t->TimeEC(i);
                        myPart.ECpath = t->PathEC(i);
                        myPart.EChit_M2 = t->EChit_Moment2(i);
                        myPart.EChit_M3 = t->EChit_Moment3(i);
                        myPart.EChit_M4 = t->EChit_Moment4(i);
                        myPart.Chi2EC = t->Chi2EC(i);

                        myPart.SCtime = t->TimeSC(i);
                        myPart.SCpath = t->PathSC(i);

                        myPart.CCnphe = t->Nphe(i);

                        if(myPart.Pid == GetPID("Electron",kind)) myPart.TimeCorr4 = t -> TimeCorr4(0.000511,i);
                        if(myPart.Pid == GetPID("PiPlus",kind)) myPart.TimeCorr4 = t -> TimeCorr4(kMassPi_plus,i);
                        if(myPart.Pid == GetPID("PiMinus",kind)) myPart.TimeCorr4 = t -> TimeCorr4(kMassPi_min,i);
                        if(myPart.Pid == GetPID("Photon",kind)) myPart.TimeCorr4 = t -> TimeCorr4(0.0,i);
                    }

                    if(myPart.Pid == GetPID("Electron",kind)) myElec = myPart;
                    if(myPart.Pid == GetPID("PiPlus",kind)) myPip = myPart;
                    if(myPart.Pid == GetPID("PiMinus",kind)) myPim = myPart;
                    if(myPart.Pid == GetPID("Photon",kind)){
                        photonCtr++;
                        if(photonCtr==1) myPhoton1 = myPart;
                        if(photonCtr==2) myPhoton2 = myPart;
                    }
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

    output->Write();
    output->Close();
    
    float timeStop = clock();
    PrintAnalysisTime(timeStart,timeStop);

    return 0;
}

void PrintUsage(char *processName)
{
    cerr << processName << " <options> <filename>\n";
    cerr << "\toptions are:\n";
    cerr << "\t-o<filename>\tROOT output file (def. = PipPimPi0.root).\n";
    cerr << "\t-M#\t\tprocess maximum # of events.\n";
    cerr << "\t-D#\t\tinform user when # of events have been processed (def. = 1000).\n";
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
        }else if(partName.compare("Photon")==0){
            ret = PDG_PHOTON;
        }else if(partName.compare("PiPlus")==0){
            ret = PDG_PIPLUS;
        }else if(partName.compare("PiMinus")==0){
            ret = PDG_PIMINUS;
        }else{
            cerr<<"GetPid(): Unknown PDG particle "<<partName.c_str()<<endl; exit(0);
        }
    }else if(kind==1){
        if(partName.compare("Electron")==0){
            ret = GEANT3_ELECTRON;
        }else if(partName.compare("Photon")==0){
            ret = GEANT3_PHOTON;
        }else if(partName.compare("PiPlus")==0){
            ret = GEANT3_PIPLUS;
        }else if(partName.compare("PiMinus")==0){
            ret = GEANT3_PIMINUS;
        }else{
            cerr<<"GetPid(): Unknown GEANT3 particle "<<partName.c_str()<<endl; exit(0);
        }
    }else{
        cerr<<"GetPID: Unknown analysis channel "<<kind<<endl;
    }
    return ret;
}