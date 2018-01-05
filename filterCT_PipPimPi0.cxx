#include "Riostream.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
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
    int nRows, kind, tempPid;
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
    
    bool topology = false;
    vector<int> partIndex;

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
  
    std::vector<int> pEvtNum;
    std::vector<int> nPart;
    std::vector<int> Sector;
    std::vector<Double_t> Charge;
    std::vector<Double_t> Pid;
    std::vector<Double_t> Beta;
    std::vector<Double_t> Px;
    std::vector<Double_t> Py;
    std::vector<Double_t> Pz;
    std::vector<Double_t> Mom;
    std::vector<Double_t> Mass2;
    std::vector<Double_t> X;
    std::vector<Double_t> Y;
    std::vector<Double_t> Z;
    std::vector<Double_t> ECx;
    std::vector<Double_t> ECy;
    std::vector<Double_t> ECz;
    std::vector<Double_t> ECu;
    std::vector<Double_t> ECv;
    std::vector<Double_t> ECw;
    std::vector<Double_t> ECtot;
    std::vector<Double_t> ECin;
    std::vector<Double_t> ECout;
    std::vector<Double_t> ECtime;
    std::vector<Double_t> ECpath;
    std::vector<Double_t> EChit_M2;
    std::vector<Double_t> EChit_M3;
    std::vector<Double_t> EChit_M4;
    std::vector<Double_t> Chi2EC;
    std::vector<Double_t> SCpath;
    std::vector<Double_t> SCtime;
    std::vector<Double_t> CCnphe;
    std::vector<Double_t> T;
    std::vector<Double_t> Xf;
    std::vector<Double_t> Mx2;
    std::vector<Double_t> Pt;
    std::vector<Double_t> Zh;
    std::vector<Double_t> ThetaPQ;
    std::vector<Double_t> PhiPQ;
    std::vector<Double_t> TimeCorr4;
    
    string kineList = "EvtNum/F:ElecVertTarg/F:Q2/F:Nu/F:Xb/F:W:Xcorr/F:Ycorr/F:Zcorr/F:nElec/I:nPip/I:nPim/I:nGam/I:nProton/I:nNeutron/I:nKp/I:nKm/I:nPositron/I";
    KINEVAR myKine;
    
    TTree *dataTree = new TTree("Data","Experimental Data Tree");
    dataTree->Branch("Kinematics",&myKine,kineList.c_str());

    dataTree->Branch("Sector",&Sector);
    dataTree->Branch("Charge",&Charge);
    dataTree->Branch("Pid",&Pid);
    dataTree->Branch("Beta",&Beta);
    dataTree->Branch("Px",&Px);
    dataTree->Branch("Py",&Py);
    dataTree->Branch("Pz",&Pz);
    dataTree->Branch("Mom",&Mom);
    dataTree->Branch("Mass2",&Mass2);
    dataTree->Branch("X",&X);
    dataTree->Branch("Y",&Y);
    dataTree->Branch("Z",&Z);
    dataTree->Branch("ECx",&ECx);
    dataTree->Branch("ECy",&ECy);
    dataTree->Branch("ECz",&ECz);
    dataTree->Branch("ECu",&ECu);
    dataTree->Branch("ECv",&ECv);
    dataTree->Branch("ECw",&ECw);
    dataTree->Branch("ECtot",&ECtot);
    dataTree->Branch("ECin",&ECin);
    dataTree->Branch("ECout",&ECout);
    dataTree->Branch("ECtime",&ECtime);
    dataTree->Branch("ECpath",&ECpath);
    dataTree->Branch("EChit_M2",&EChit_M2);
    dataTree->Branch("EChit_M3",&EChit_M3),
    dataTree->Branch("EChit_M4",&EChit_M4);
    dataTree->Branch("Chi2EC",&Chi2EC);
    dataTree->Branch("SCpath",&SCpath);
    dataTree->Branch("SCtime",&SCtime);
    dataTree->Branch("CCnphe",&CCnphe);
    dataTree->Branch("T",&T);
    dataTree->Branch("Xf",&Xf);
    dataTree->Branch("Mx2",&Mx2);
    dataTree->Branch("Pt",&Pt);
    dataTree->Branch("Zh",&Zh);
    dataTree->Branch("ThetaPQ",&ThetaPQ);
    dataTree->Branch("PhiPQ",&PhiPQ);
    dataTree->Branch("TimeCorr4",&TimeCorr4);

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
        pEvtNum.clear();
        nPart.clear();
        Sector.clear();
        Charge.clear();
        Pid.clear();
        Beta.clear();
        Px.clear();
        Py.clear();
        Pz.clear();
        Mom.clear();
        Mass2.clear();
        X.clear();
        Y.clear();
        Z.clear();
        ECx.clear();
        ECy.clear();
        ECz.clear();
        ECu.clear();
        ECv.clear();
        ECw.clear();
        ECtot.clear();
        ECin.clear();
        ECout.clear();
        ECtime.clear();
        ECpath.clear();
        EChit_M2.clear();
        EChit_M3.clear();
        EChit_M4.clear();
        Chi2EC.clear();
        SCpath.clear();
        SCtime.clear();
        CCnphe.clear();
        T.clear();
        Xf.clear();
        Mx2.clear();
        Pt.clear();
        Zh.clear();
        ThetaPQ.clear();
        PhiPQ.clear();
        TimeCorr4.clear();
        
        if(nRows>0){
            partIndex.clear(); // clear out the particle list
	    	topology = false; // init. the event topology cut
	    	for (j = 0; j < nRows; j++) {

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
               
//                cout<<"Particle "<< tempPid <<"\t"<<catPid<<endl;
//                if(tempPid == GetPID("Electron",kind)){
                if(catPid.EqualTo("electron")){
                    myKine.nElec++;
                    partIndex.push_back(j);
                    ctr_nElec++;
                }
                if(catPid.EqualTo("high energy pion +") || catPid.EqualTo("low energy pion +") || catPid.EqualTo("pi+")){
//                if(tempPid == GetPID("PiPlus",kind)){
                    myKine.nPip++;
                    partIndex.push_back(j);
                }
                if(catPid.EqualTo("pi-")){
//                if(tempPid == GetPID("PiMinus",kind)){
                    myKine.nPim++;
                    partIndex.push_back(j);
                }
                if(catPid.EqualTo("gamma")){
//                if(tempPid == GetPID("Photon",kind)){
                    myKine.nGam++;
                    partIndex.push_back(j);
                }
                
                if(tempPid == GetPID("Proton",kind)) myKine.nProton++;
                if(tempPid == GetPID("Neutron",kind)) myKine.nNeutron++;
                if(tempPid == GetPID("KPlus",kind)) myKine.nKp++;
                if(tempPid == GetPID("KMinus",kind)) myKine.nKm++;
                if(tempPid == GetPID("Positron",kind)) myKine.nPositron++;
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

                pEvtNum.push_back(t -> NEvent());
                nPart.push_back(partIndex.size()); // save the number of particles
        		while (!partIndex.empty()) {
		    		i = partIndex.back(); // retrieve EVNT index for each particle
                    partIndex.pop_back(); // erase last entry in the list

                    Sector.push_back(t->Sector(i,kind));
                    Charge.push_back(t->Charge(i,kind));
                    Beta.push_back(t->Betta(i,kind));
                    Pid.push_back(t->Id(i,kind));
                    Mom.push_back(t->Momentum(i,kind));
                    Px.push_back(t->Px(i, kind));
                    Py.push_back(t->Py(i, kind));
                    Pz.push_back(t->Pz(i, kind));
                    X.push_back(t->X(i, kind));
                    Y.push_back(t->Y(i, kind));
                    Z.push_back(t->Z(i, kind));
                    Mass2.push_back(t->Mass2(i,kind));
                    ThetaPQ.push_back(t -> ThetaPQ(i, kind));
                    PhiPQ.push_back(t -> PhiPQ(i, kind));
                    Zh.push_back(t -> Zh(i, kind));
                    Pt.push_back(TMath::Sqrt(t -> Pt2(i, kind)));
                    Mx2.push_back(t -> Mx2(i, kind));
                    Xf.push_back(t -> Xf(i, kind));
                    T.push_back(t -> T(i, kind));

                    if(simul_key == false){
                        ECtot.push_back(TMath::Max(t->Etot(i),t->Ein(i)+t->Eout(i)));
                        ECin.push_back(t->Ein(i));
                        ECout.push_back(t->Eout(i));
                        ECx.push_back(t->XEC(i));
                        ECy.push_back(t->YEC(i));
                        ECz.push_back(t->ZEC(i));
                        ECxyz->SetXYZ(t->XEC(i),t->YEC(i),t->ZEC(i));
                        ECuvw = t->XYZToUVW(ECxyz);
                        ECu.push_back(ECuvw->X());
                        ECv.push_back(ECuvw->Y());
                        ECw.push_back(ECuvw->Z());
                        ECtime.push_back(t->TimeEC(i));
                        ECpath.push_back(t->PathEC(i));
                        EChit_M2.push_back(t->EChit_Moment2(i));
                        EChit_M3.push_back(t->EChit_Moment3(i));
                        EChit_M4.push_back(t->EChit_Moment4(i));
                        Chi2EC.push_back(t->Chi2EC(i));
                        SCtime.push_back(t->TimeSC(i));
                        SCpath.push_back(t->PathSC(i));
                        CCnphe.push_back(t->Nphe(i));

                        if(t->Id(i,kind) == GetPID("Electron",kind)) TimeCorr4.push_back(t -> TimeCorr4(0.000511,i));
                        if(t->Id(i,kind) == GetPID("PiPlus",kind)) TimeCorr4.push_back(t -> TimeCorr4(kMassPi_plus,i));
                        if(t->Id(i,kind) == GetPID("PiMinus",kind)) TimeCorr4.push_back(t -> TimeCorr4(kMassPi_min,i));
                        if(t->Id(i,kind) == GetPID("Photon",kind)) TimeCorr4.push_back(t -> TimeCorr4(0.0,i));
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
    cout<<"#(e-) = "<<ctr_nElec<<endl;
                   
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
