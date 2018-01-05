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

#pragma link C++ class vector<float>+;

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

#define MAX_PART 15 // max number of particles in each event

//declarations of functions
void PrintAnalysisTime(float tStart, float tStop);
void PrintUsage(char *processName);
int GetPID(string partName, int kind);

typedef struct{
    Float_t EvtNum, ElecVertTarg, Q2, Nu, Xb, W;
    Float_t Xcorr, Ycorr, Zcorr;
    Int_t nElec, nPip, nPim, nGam, nProton, nNeutron, nKp, nKm, nPositron;
} KINEVAR;

typedef struct{
    Float_t pEvtNum;
    int nPart;
    int Sector[MAX_PART];
    Double_t Charge[MAX_PART];
    Double_t Pid[MAX_PART], Beta[MAX_PART];
    Double_t Px[MAX_PART], Py[MAX_PART], Pz[MAX_PART], Mom[MAX_PART], Mass2[MAX_PART];
    float X[MAX_PART], Y[MAX_PART], Z[MAX_PART];
    float ECx[MAX_PART], ECy[MAX_PART], ECz[MAX_PART], ECu[MAX_PART], ECv[MAX_PART], ECw[MAX_PART];
    float ECtot[MAX_PART], ECin[MAX_PART], ECout[MAX_PART], ECtime[MAX_PART], ECpath[MAX_PART];
    float EChit_M2[MAX_PART], EChit_M3[MAX_PART], EChit_M4[MAX_PART], Chi2EC[MAX_PART];
    float SCpath[MAX_PART], SCtime[MAX_PART];
    float CCnphe[MAX_PART];
    float T[MAX_PART], Xf[MAX_PART], Mx2[MAX_PART], Pt[MAX_PART], Zh[MAX_PART], ThetaPQ[MAX_PART], PhiPQ[MAX_PART], TimeCorr4[MAX_PART];
} PARTVAR;

int main(int argc, char **argv)
{
    extern char *optarg;
    int c;
    extern int optind;
    
    int i, j, k;
    int nRows, kind, tempPid;
    int pCtr;
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
    
    string kineList = "EvtNum/F:ElecVertTarg/F:Q2/F:Nu/F:Xb/F:W:Xcorr/F:Ycorr/F:Zcorr/F:nElec/I:nPip/I:nPim/I:nGam/I:nProton/I:nNeutron/I:nKp/I:nKm/I:nPositron/I";
 
    string partList = "pEvtNum/F:nPart/I:Charge[nPart]/D:Sector[nPart]/I:Pid[nPart]/D:Beta[nPart]/D:Px[nPart]/D:Py[nPart]/D:Pz[nPart]/D:Mom[nPart]/D:Mass2[nPart]/D:X[nPart]/F:Y[nPart]/F:Z[nPart]/F:ECx[nPart]/F:ECy[nPart]/F:ECz[nPart]/F:ECu[nPart]/F:ECv[nPart]/F:ECw[nPart]/F:ECtot[nPart]/F:ECin[nPart]/F:ECout[nPart]/F:ECtime[nPart]/F:ECpath[nPart]/F:EChit_M2[nPart]/F:EChit_M3[nPart]/F:EChit_M4[nPart]/F:Chi2EC[nPart]/F:SCpath[nPart]/F:SCtime[nPart]/F:CCnphe[nPart]/F:T[nPart]/F:Xf[nPart]/F:Mx2[nPart]/F:Pt[nPart]/F:Zh[nPart]/F:ThetaPQ[nPart]/F:PhiPQ[nPart]/F:TimeCorr4[nPart]/F";
    KINEVAR myKine;
    PARTVAR myPart;
    
    TTree *dataTree = new TTree("Data","Experimental Data Tree");
    dataTree->Branch("Kinematics",&myKine,kineList.c_str());
    dataTree->Branch("Particle",&myPart,partList.c_str());
    dataTree->Branch("nPart",&myPart.nPart,"nPart/I");
    dataTree->Branch("Sector",&myPart.Sector,"Sector[nPart]/I");
    dataTree->Branch("Charge",&myPart.Charge,"Charge[nPart]/F");

    vector<int> vSector;
    dataTree->Branch("vSector",&vSector);
    std::vector<Double_t> vCharge;
    std::vector<Double_t> vPid;
    std::vector<Double_t> vBeta;
    std::vector<Double_t> vPx;
    std::vector<Double_t> vPy;
    std::vector<Double_t> vPz;
    std::vector<Double_t> vMom;
    std::vector<Double_t> vMass2;
    std::vector<Double_t> vX;
    std::vector<Double_t> vY;
    std::vector<Double_t> vZ;
    std::vector<Double_t> vECx;
    std::vector<Double_t> vECy;
    std::vector<Double_t> vECz;
    std::vector<Double_t> vECu;
    std::vector<Double_t> vECv;
    std::vector<Double_t> vECw;
    std::vector<Double_t> vECtot;
    std::vector<Double_t> vECin;
    std::vector<Double_t> vECout;
    std::vector<Double_t> vECtime;
    std::vector<Double_t> vECpath;
    std::vector<Double_t> vEChit_M2;
    std::vector<Double_t> vEChit_M3,
    std::vector<Double_t> vEChit_M4;
    std::vector<Double_t> vChi2EC;
    std::vector<Double_t> vSCpath;
    std::vector<Double_t> vSCtime;
    std::vector<Double_t> vCCnphe;
    std::vector<Double_t> vT;
    std::vector<Double_t> vXf;
    std::vector<Double_t> vMx2;
    std::vector<Double_t> vPt;
    std::vector<Double_t> vZh;
    std::vector<Double_t> vThetaPQ;
    std::vector<Double_t> vPhiPQ;
    std::vector<Double_t> vTimeCorr4;
    dataTree->Branch("vCharge",&vCharge);
    dataTree->Branch("vPid",&vPid);
    dataTree->Branch("vBeta",&vBeta);
    dataTree->Branch("vPx",&vPx);
    dataTree->Branch("vPy",&vPy);
    dataTree->Branch("vPz",&vPz);
    dataTree->Branch("vMom",&vMom);
    dataTree->Branch("vMass2",&Mass2);
    dataTree->Branch("vX",&vX);
    dataTree->Branch("vY",&vY);
    dataTree->Branch("vZ",&vZ);
    dataTree->Branch("vECx",&vECx);
    dataTree->Branch("vECy",&vECy);
    dataTree->Branch("vECz",&vECz);
    dataTree->Branch("vECu",&vECu);
    dataTree->Branch("vECv",&vECv);
    dataTree->Branch("vECw",&vECw);
    dataTree->Branch("vECtot",&vECtot);
    dataTree->Branch("vECin",&vECin);
    dataTree->Branch("vECout",&vECout);
    dataTree->Branch("vECtime",&vECtime);
    dataTree->Branch("vECpath",&vECpath);
    dataTree->Branch("vEChit_M2",&vEChit_M2);
    dataTree->Branch("vEChit_M3",&vEChiy_M3),
    dataTree->Branch("vEChit_M4",&vEChiy_M4);
    dataTree->Branch("vChi2EC",&vChi2EC);
    dataTree->Branch("vSCpath",&vSCpath);
    dataTree->Branch("vSCtime",&vSCtime);
    dataTree->Branch("vCCnphe",&vCCnphe);
    dataTree->Branch("vT",&vT);
    dataTree->Branch("vXf",&vXf);
    dataTree->Branch("vMx2",&vMx2);
    dataTree->Branch("vPt",&vPt);
    dataTree->Branch("vZh",&vZh);
    dataTree->Branch("vThetaPQ",&vThetaPQ);
    dataTree->Branch("vPhiPQ",&vPhiPQ);
    dataTree->Branch("vTimeCorr4",&vTimeCorr4);

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
        
//       cout<<"Event "<<k+1<<endl;
        
        memset(&myKine,0,sizeof(myKine)); // init kinematics struct to zeros
        memset(&myPart,0,sizeof(myPart)); // init particle struct to zeros
        vSector.clear();
        vCharge.clear();
        vPid.clear();
        vBeta.clear();
        vPx.clear();
        vPy.clear();
        vPz.clear();
        vMom.clear();
        vMass2.clear();
        vX.clear();
        vY.clear();
        vZ.clear();
        vECx.clear();
        vECy.clear();
        vECz.clear();
        vECu.clear();
        vECv.clear();
        vECw.clear();
        vECtot.clear();
        vECin.clear();
        vECout.clear();
        vECtime.clear();
        vECpath.clear();
        vEChit_M2.clear();
        vEChit_M3.clear();
        vEChit_M4.clear();
        vChi2EC.clear();
        vSCpath.clear();
        vSCtime.clear();
        vCCnphe.clear();
        vT.clear();
        vXf.clear();
        vMx2.clear();
        vPt.clear();
        vZh.clear();
        vThetaPQ.clear();
        vPhiPQ.clear();
        vTimeCorr4.clear();
        
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

	    		pCtr = 0;  // initialize the partcle counter
                myPart.pEvtNum = t -> NEvent();
                myPart.nPart = partIndex.size(); // save the number of particles
        		while (!partIndex.empty()) {
		    		i = partIndex.back(); // retrieve EVNT index for each particle
                    partIndex.pop_back(); // erase last entry in the list

                    vSector.push_back(t->Sector(i,kind));
                    vCharge.push_back(t->Charge(i,kind));
                    vBeta.push_back(t->Betta(i,kind));
                    vPid.push_back(t->Id(i,kind));
                    vMom.push_back(t->Momentum(i,kind));
                    vPx.push_back(t->Px(i, kind));
                    vPy.push_back(t->Py(i, kind));
                    vPz.push_back(t->Pz(i, kind));
                    vX.push_back(t->X(i, kind));
                    vY.push_back(t->Y(i, kind));
                    vZ.push_back(t->Z(i, kind));
                    vMass2.push_back(t->Mass2(i,kind));
                    vThetaPQ.push_back(t -> ThetaPQ(i, kind));
                    vPhiPQ.push_back(t -> PhiPQ(i, kind));
                    vZh.push_back(t -> Zh(i, kind));
                    vPt.push_back(TMath::Sqrt(t -> Pt2(i, kind)));
                    vMx2.push_back(t -> Mx2(i, kind));
                    vXf.push_back(t -> Xf(i, kind));
                    vT.push_back(t -> T(i, kind));
                    
                    myPart.Sector[pCtr] = t->Sector(i,kind);
                    myPart.Charge[pCtr] = t->Charge(i,kind);
                    myPart.Beta[pCtr] = t->Betta(i,kind);
                    myPart.Pid[pCtr] = t->Id(i,kind);
                    myPart.Mom[pCtr] = t->Momentum(i,kind);
                    myPart.Px[pCtr] = t->Px(i, kind);
                    myPart.Py[pCtr] = t->Py(i, kind);
                    myPart.Pz[pCtr] = t->Pz(i, kind);
                    myPart.X[pCtr] = t->X(i, kind);
                    myPart.Y[pCtr] = t->Y(i, kind);
                    myPart.Z[pCtr] = t->Z(i, kind);
                    myPart.Mass2[pCtr] = t->Mass2(i,kind);

                    myPart.ThetaPQ[pCtr] = t -> ThetaPQ(i, kind);
                    myPart.PhiPQ[pCtr] = t -> PhiPQ(i, kind);
                    myPart.Zh[pCtr] = t -> Zh(i, kind);
                    myPart.Pt[pCtr] = TMath::Sqrt(t -> Pt2(i, kind));
                    myPart.Mx2[pCtr] = t -> Mx2(i, kind);
                    myPart.Xf[pCtr] = t -> Xf(i, kind);
        			myPart.T[pCtr] = t -> T(i, kind);
                    
                    if(simul_key == false){
                        vECtot.push_back(TMath::Max(t->Etot(i),t->Ein(i)+t->Eout(i)));
                        vECin.push_back(t->Ein(i));
                        vECout.push_back(t->Eout(i));
                        vECx.push_back(t->XEC(i));
                        vECy.push_back(t->YEC(i));
                        vECz.push_back(t->ZEC(i));
                        ECxyz->SetXYZ(t->XEC(i),t->YEC(i),t->ZEC(i));
                        ECuvw = t->XYZToUVW(ECxyz);
                        vECu.push_back(ECuvw->X());
                        vECv.push_back(ECuvw->Y());
                        vECw.push_back(ECuvw->Z());
                        vECtime.push_back(t->TimeEC(i));
                        vECpath.push_back(t->PathEC(i));
                        vEChit_M2.push_back(t->EChit_Moment2(i));
                        vEChit_M3.push_back(t->EChit_Moment3(i));
                        vEChit_M4.push_back(t->EChit_Moment4(i));
                        vChi2EC.push_back(t->Chi2EC(i));
                        vSCtime.push_back(t->TimeSC(i));
                        vSCpath.push_back(t->PathSC(i));
                        vCCnphe.push_back(t->Nphe(i));
                        
                        myPart.ECtot[pCtr] = TMath::Max(t->Etot(i),t->Ein(i)+t->Eout(i));
                        myPart.ECin[pCtr] = t->Ein(i);
                        myPart.ECout[pCtr] = t->Eout(i);
                        myPart.ECx[pCtr] = t->XEC(i);
                        myPart.ECy[pCtr] = t->YEC(i);
                        myPart.ECz[pCtr] = t->ZEC(i);
                        ECxyz->SetXYZ(t->XEC(i),t->YEC(i),t->ZEC(i));
                        ECuvw = t->XYZToUVW(ECxyz);
                        myPart.ECu[pCtr] = ECuvw->X();
                        myPart.ECv[pCtr] = ECuvw->Y();
                        myPart.ECw[pCtr] = ECuvw->Z();
 				        myPart.ECtime[pCtr] = t->TimeEC(i);
                        myPart.ECpath[pCtr] = t->PathEC(i);
                        myPart.EChit_M2[pCtr] = t->EChit_Moment2(i);
                        myPart.EChit_M3[pCtr] = t->EChit_Moment3(i);
                        myPart.EChit_M4[pCtr] = t->EChit_Moment4(i);
                        myPart.Chi2EC[pCtr] = t->Chi2EC(i);

                        myPart.SCtime[pCtr] = t->TimeSC(i);
                        myPart.SCpath[pCtr] = t->PathSC(i);

                        myPart.CCnphe[pCtr] = t->Nphe(i);

                        if(myPart.Pid[pCtr] == GetPID("Electron",kind)) myPart.TimeCorr4[pCtr] = t -> TimeCorr4(0.000511,i);
                        if(myPart.Pid[pCtr] == GetPID("PiPlus",kind)) myPart.TimeCorr4[pCtr] = t -> TimeCorr4(kMassPi_plus,i);
                        if(myPart.Pid[pCtr] == GetPID("PiMinus",kind)) myPart.TimeCorr4[pCtr] = t -> TimeCorr4(kMassPi_min,i);
                        if(myPart.Pid[pCtr] == GetPID("Photon",kind)) myPart.TimeCorr4[pCtr] = t -> TimeCorr4(0.0,i);
                        vTimeCorr4.push_back(myPart.TimeCorr4[pCtr]);
                    }
                    pCtr++;
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
