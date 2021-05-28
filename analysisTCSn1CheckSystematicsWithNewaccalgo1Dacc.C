#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "bib/TCSclass.h"
#include "bib/TCSfunc.h"
#include "reader.h"
#include "/vol0/pierre/Bureau/Hipo4Ana/Systematics/QAdb/srcC/include/QADB.h"


using namespace std;


Particle RadiativeCorr(Particle vParticle, Particle Photons[], double thetaWin,double thetaWin1, int np){

	for(int i=0;i<np;i++){
		if(Photons[i].Vector.Angle(vParticle.Vector.Vect())*TMath::RadToDeg()<thetaWin && abs(Photons[i].Vector.Theta()-vParticle.Vector.Theta())*TMath::RadToDeg()<thetaWin1 && Photons[i].Vector.P()>0.){
		vParticle.Vector=(vParticle.Vector+Photons[i].Vector);
		//cout<<"in rad corr "<<endl;
		}
	}
 return vParticle;
}


double smear(TH1F histoIn, int bin){

	double maxrange= histoIn.GetXaxis()->GetXmax();
	double minrange= histoIn.GetXaxis()->GetXmin();
	int binNum = (histoIn.GetSize());
	int nstepSmearing=10000;
	TH1F *patate= new TH1F(Form("amp%d",bin),";R^{R};Counts",300,-1.,1.);
	cout<<"start smear "<<endl;
	for(int it=0;it<nstepSmearing;it++){
		TRandom random = TRandom(it);

		double Rc=0.0;
		double R0=0.0;
		double R=0.0;


		for(int i=1; i<binNum-1; i++){
			double meanY=histoIn.GetBinContent(i);
			double phi = histoIn.GetBinCenter(i);
			//cout<<"phi "<<phi<<" mean "<<meanY<<endl;
			double sigmaY=histoIn.GetBinError(i);
			double valueY= random.Gaus(meanY,sigmaY);
			
			Rc+=valueY*cos(phi*TMath::DegToRad());
			R0+=valueY;

		}

		//cout<<bin<<" "<<Rc<<" "<<R0<<endl;

		R=Rc/R0;

		patate->Fill(R);
		//cout<<R<<endl;
		
	}

	
	gStyle->SetOptStat(1111);
	TCanvas *patatec  = new TCanvas("aaddr ","A",3000,2000);
	//patate->Draw();
	patate->Fit("gaus");//,"q0"
	patatec->SaveAs(Form("amp%d.pdf",bin));
	gStyle->SetOptStat(111);
	cout<<"end smear "<<endl;
	cout<<patate->GetFunction("gaus")->GetParameter(1)<<endl;
	cout<<patate->GetFunction("gaus")->GetParameter(2)<<endl;
	return patate->GetFunction("gaus")->GetParameter(2);

}

double smear(TH1F histoIn, TH1F EhistoIn, int bin){

	double maxrange= histoIn.GetXaxis()->GetXmax();
	double minrange= histoIn.GetXaxis()->GetXmin();
	int binNum = (histoIn.GetSize());
	int nstepSmearing=10;
	TH1F *patate= new TH1F(Form("amp%d",bin),"amplitude;R^{R};Counts",300,-2.,2.);
	cout<<"start smear "<<endl;
	for(int it=0;it<nstepSmearing;it++){
		TRandom random = TRandom(it);

		double Rc=0.0;
		double R0=0.0;
		double R=0.0;


		for(int i=1; i<binNum-1; i++){
			double meanY=histoIn.GetBinContent(i);
			double phi = histoIn.GetBinCenter(i);
			double sigmaY=EhistoIn.GetBinContent(i);
			double valueY= random.Gaus(meanY,sigmaY);
			
			Rc+=valueY*cos(phi*TMath::DegToRad());
			R0+=valueY;

		}

		cout<<bin<<" "<<Rc<<" "<<R0<<endl;

		R=Rc/R0;

		patate->Fill(R);
		//cout<<R<<endl;
		
	}

	patate->Fit("gaus","q0");

	TCanvas *patatec  = new TCanvas("aaddr ","A",4000,2000);
	patate->Draw();
	patatec->SaveAs(Form("amp%d.pdf",bin));

	cout<<"end smear "<<endl;
	cout<<patate->GetFunction("gaus")->GetParameter(1)<<endl;
	cout<<patate->GetFunction("gaus")->GetParameter(2)<<endl;
	return patate->GetFunction("gaus")->GetParameter(2);

}

PairValue smearFit(TH1F histoIn, int bin){

	int binNum = (histoIn.GetSize());
	int nstepSmearing=10000;
	TH1F *patate= new TH1F(Form("amp fit %d",bin),"amplitude;Fit;Counts",1000,-1.0,0.5);
	cout<<"start smear fit"<<endl;
	
	PairValue results;
	for(int it=0;it<nstepSmearing;it++){
		TRandom random = TRandom(it);



		TH1F *fitTemp= new TH1F(Form("temp%d%d",it,bin),Form("temp%d",it),10,-180,180);
		for(int i=1; i<binNum-1; i++){
			double meanY=histoIn.GetBinContent(i);
			double sigmaY=histoIn.GetBinError(i);
			
			double valueY= random.Gaus(meanY,sigmaY);
			fitTemp->SetBinContent(i,valueY);
			

		}

		TF1  *f2 = new TF1("f2","[0]*sin(x*[1]/[2])",-180,180);
		f2->FixParameter(1,3.14159264);
		f2->FixParameter(2,180);
		fitTemp->Fit("f2","q0");//->Draw();//
		
		
		double fitresult = fitTemp->GetFunction("f2")->GetParameter(0);
		double fiterror = fitTemp->GetFunction("f2")->GetParError(0);
		
		patate->Fill(fitresult);

		
	}

	
	gStyle->SetOptStat(1111);
	TCanvas *patatec  = new TCanvas("aaddr ","A",4000,2000);
	//patate->Draw();
	patate->Draw();//Fit("gaus","q");
	patatec->SaveAs(Form("ampfit%d.pdf",bin));
	gStyle->SetOptStat(111);
	cout<<"end smear "<<endl;
	cout<<patate->GetMean()<<endl;
	cout<<patate->GetStdDev()<<endl;
	results.mean=patate->GetMean();//->GetFunction("gaus")->GetParameter(1);
	results.sigma=patate->GetStdDev();//->GetFunction("gaus")->GetParameter(2);
	return results;

}


int analysisTCSn1CheckSystematicsWithNewaccalgo1Dacc() {


	//QADB * qa = new QADB();

	gStyle->SetOptStat(111);
	gStyle->SetPalette(55);	
	gStyle->SetLabelSize(.05, "xyz");
	gStyle->SetTitleSize(.05, "xyz");
	gStyle->SetTitleSize(.07,"t");
	/*gStyle->SetFrameLineWidth(2); 
	gStyle->SetLineWidth(2);
	gStyle->SetHistLineWidth(2);*/
	gStyle->SetMarkerStyle(13);

	Int_t argc=gApplication->Argc();
	char** argv=gApplication->Argv();
	bool MC=false;
	double nbrecEvent=0;
	int nbf=0;
	double nEventTCS=0;
	double denom=0;

	double nCD=0;
	double nFD=0;

	int AfterCuts=0;

	TString nameFiles="";

	double sum1=0.0;
	double sumw21=0.0;
	double sum2=0.0;
	double sumw22=0.0;


	TString type="REC";
	int tbin=4;	
	int binPhi=20;
	int binTheta=20;
	double mintheta=0;
	double maxtheta=140;
	double minphi=-180;
	double maxphi=180;

	double sumPhi[23][6]={0.};//one sum for each theta+1slot for number of event per bin//one line for each t
	double ErrorSsumPhi[23][6]={0.};// Statistical error //one error SQUARED for each phi + 1slot for number of event per bin//one line for each t
	double ErrorAcceptance[23][6]={0.};//Acceptance error
	//double tbins[3]={0.340701, 0.485526, 0.636489};

	//TH2D *AcceptancetM2forbinning[3];
	Double_t t1 [5] = {0.15,0.25, 0.34, 0.48,  0.8};//{0.15,0.3, 0.45, 0.6,  0.8};//{0.15,0.247065, 0.327602, 0.452239,  0.8}; ; // ; //
	Double_t t1b [5] = {0.15,0.35, 0.45, 0.55,  0.8};
	Double_t M1 [5] = {1.5, 1.7 ,2.,2.5,3.};//{2.25, 2.5 ,3.,4.,9.};//{2.25, 4. ,5. ,6.,9.}; // 
	/*Double_t t2 [5] = {0.15,0.247065, 0.327602, 0.452239,  0.8}; ; 
	Double_t M2 [5] = {2.25, 2.60732, 3.05203 ,3.82303,9.};
	Double_t t3 [5] = {0.15,0.247065, 0.327602, 0.452239,  0.8}; ; 
	Double_t M3 [5] = {2.25, 2.87779, 3.54478 ,4.50211  ,9.};*/
	TH2D *PhysicstM2forbinning= new TH2D("PhysicstM2forbinning","",4,t1,4,M1);
	/*AcceptancetM2forbinning[1]= new TH2D("AcceptancetM2forbinning1","",4,t2,4,M2);
	AcceptancetM2forbinning[2]= new TH2D("AcceptancetM2forbinning2","",4,t3,4,M3);*/
	//TH1D *AcceptanceEg = new TH1D("AcceptanceEg",";Eg",3,4,11);
	
	Double_t Xi3 [4] = {0.0, 0.12, 0.15  ,0.4};
	TH1D *PhysicsXi = new TH1D("AcceptanceXi",";Xi",3,Xi3);
	
	Double_t Eg1 [4] = {4.,6.4, 8., 10.6};
	TH1D *PhysicsEg1 = new TH1D("Eg1",";Eg",3,Eg1);
	
	TH1D *CheckTbinning = new TH1D("CheckTbinning",";t",4,t1);
	//New Way to do acceptance
	
	
	
	TFile *acceptanceFile = new TFile(TString(argv[4]));
	
	TH2D *AcceptancetM2forbinningHist = (TH2D*)acceptanceFile->Get("Egbin0");  
	TH1D *AcceptanceEg = (TH1D*)acceptanceFile->Get("AcceptanceEg");  
		cout<<"here"<<endl;
	const int nbBinsInT = AcceptancetM2forbinningHist->GetNbinsX();
	cout<<nbBinsInT<<endl;
	const int nbBinsInM = AcceptancetM2forbinningHist->GetNbinsY();
	const int nbBinsInEg = AcceptanceEg->GetNbinsX();
	
	
	TH2D *AcceptancetM2forbinning[nbBinsInEg];
	for(int Egbin=0;Egbin<nbBinsInEg;Egbin++){
		AcceptancetM2forbinning[Egbin]= AcceptancetM2forbinningHist;//new TH2D(Form("AcceptancetM2forbinning%d",Egbin),"",4,t1a,4,M1a);
		AcceptancetM2forbinning[Egbin]->SetName(Form("Egbin%d",Egbin));
		AcceptancetM2forbinning[Egbin]->SetTitle(Form("Eg bin %d;-t;M2",Egbin));
		
	}
	

	
	gROOT->SetBatch(kTRUE);
	TH2D **Acc;
	Acc=new TH2D*[nbBinsInT*nbBinsInM*nbBinsInEg];//[4][3]
	//TH2D *Acc;
	
	cout<<nbBinsInEg<<" "<<nbBinsInM<<" "<<nbBinsInT<<endl;
	for(int Egbin=0;Egbin<nbBinsInEg;Egbin++){
		for(int mass=0;mass<nbBinsInM;mass++){
			for(int tbinning=0;tbinning<nbBinsInT;tbinning++){
			TString numt=TString::Itoa(tbinning,10);
			TString numm=TString::Itoa(mass,10);
			TString numE=TString::Itoa(Egbin,10);
			TString name=numm+numt+numE;
			cout<<name<<endl;
			Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]=(TH2D*)acceptanceFile->Get(name);    
					
					
					 
			}
		}
	}
		
	const int binInPhi=Acc[0]->GetNbinsX();
	const int binInTheta=Acc[0]->GetNbinsY();
	
	cout<<" "<<nbBinsInEg<<" "<<nbBinsInM<<" "<<nbBinsInT<<" "<<binInPhi<<" "<<binInTheta<<endl;
	
	
	TCanvas *AcccanNewAcc0  = new TCanvas("AcccanNewAcc","",7000,5000);
	AcccanNewAcc0->Divide(nbBinsInT,nbBinsInM);
	for(int mass=0;mass<nbBinsInM;mass++){
			for(int tbinning=0;tbinning<nbBinsInT;tbinning++){
			int Egbin =2;
				AcccanNewAcc0->cd(1+tbinning+nbBinsInT*mass);
	Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]->GetYaxis()->SetTitle("#theta (#circ)");
	Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]->GetXaxis()->SetTitle("#phi (#circ)");
	Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]->GetYaxis()->SetTitleOffset(0.9);
	Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]->GetZaxis()->SetLabelSize(0.025);
	Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]->Draw("colz");
	}
	}
	AcccanNewAcc0->SaveAs("acccheck0.pdf");
	Acc[2+nbBinsInT*1+nbBinsInT*nbBinsInM*2]->SaveAs("AccPlotNice.root");
	
	
	TH1F *AcceptanceError = new TH1F("errorAcc","; Error;nb bins",100,0,1);
	TCanvas *AcceptanceErrorAcc  = new TCanvas("AcceptanceErrorAcc","",6500,5000);
	for(int Egbin=0;Egbin<nbBinsInEg;Egbin++){
		for(int mass=0;mass<nbBinsInM;mass++){
			for(int tbinning=0;tbinning<nbBinsInT;tbinning++){
				for(int phiBin=1;phiBin<37;phiBin++){
					for(int thetabin=1;thetabin<14;thetabin++){
					double AccValue = Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]->GetBinContent(phiBin,thetabin);
					double AccError = Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]->GetBinError(phiBin,thetabin);
					if(AccValue>0.05)AcceptanceError->Fill(AccError/AccValue);
					}	
				} 
			}
		}
	}
	AcceptanceError->Draw();
	AcceptanceErrorAcc->SaveAs("AcceptanceErrorAcc.pdf");
		AcceptanceErrorAcc->SaveAs("AcceptanceErrorAcc.root");
	
		//calculate bin volume
	double limitminPhi1 = -40.;
	double limitmaxPhi1 = 40.;
	double limitminTheta1 = 50.;
	double limitmaxTheta1 = 80.;
	
	double limitminPhi2 = 140.;
	double limitmaxPhi2 = 220.;
	double limitminTheta2 = 100.;
	double limitmaxTheta2 = 130.;
	
	/*double limitminPhi1 = -40.;
	double limitmaxPhi1 = 40.;
	double limitminTheta1 = 50.;
	double limitmaxTheta1 = 80.;
	
	double limitminPhi2 = 140.;
	double limitmaxPhi2 = 220.;
	double limitminTheta2 = 100.;
	double limitmaxTheta2 = 130.;*/
	
	/*double limitminPhi1 = -20.;
	double limitmaxPhi1 = 20.;
	double limitminTheta1 = 60.;
	double limitmaxTheta1 = 70.;
	
	double limitminPhi2 = 160.;
	double limitmaxPhi2 = 200.;
	double limitminTheta2 = 110.;
	double limitmaxTheta2 = 120.;*/
	
	double volume1[nbBinsInT*nbBinsInM*nbBinsInEg];
	double volume2[nbBinsInT*nbBinsInM*nbBinsInEg];
	for(int Egbin=0;Egbin<nbBinsInEg;Egbin++){
		for(int mass=0;mass<nbBinsInM;mass++){
			for(int tbinning=0;tbinning<nbBinsInT;tbinning++){
			
				double vol1=0.;
				double vol2=0.;
				for(int phi=0;phi<100;phi++){
					for(int theta=0;theta<100;theta++){	
				 	double phiValue1= limitminPhi1+((limitmaxPhi1-limitminPhi1)/99.)*phi;
				 	double thetaValue1= limitminTheta1+((limitmaxTheta1-limitminTheta1)/99.)*theta;
				 	double phiValue2= limitminPhi2+((limitmaxPhi2-limitminPhi2)/99.)*phi;
				 	if(phiValue2>180.)phiValue2=phiValue2-360.;
				 	double thetaValue2=limitminTheta2+((limitmaxTheta2-limitminTheta2)/99.)*theta;
				 	
				 	//cout<<Egbin<<" "<<mass<<" "<<tbinning<<" "<<phiValue1<<" "<<phiValue2<<" "<<thetaValue1<<" "<<thetaValue2<<" "<<endl;
				 	
				 	double AccValue1=Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]->GetBinContent((Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]->GetXaxis()->FindBin(phiValue1)),(Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]->GetYaxis()->FindBin(thetaValue1)));
				 	double AccValue2=Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin] ->GetBinContent((Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]->GetXaxis()->FindBin(phiValue2)),(Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]->GetYaxis()->FindBin(thetaValue2)));
				 	
				 	double ErrAccValue1=Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]->GetBinError((Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]->GetXaxis()->FindBin(phiValue1)),(Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]->GetYaxis()->FindBin(thetaValue1)));
				 	double ErrAccValue2=Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin] ->GetBinError((Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]->GetXaxis()->FindBin(phiValue2)),(Acc[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]->GetYaxis()->FindBin(thetaValue2)));
				 	if(AccValue1>0.05 && (ErrAccValue1/AccValue1)<0.5 )vol1=vol1+1.;
				 	if(AccValue2>0.05 && (ErrAccValue2/AccValue2)<0.5)vol2=vol2+1.;
					} 
				}
				
				volume1[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]=(vol1/10000.);	 
				volume2[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]=(vol2/10000.);	 
			}
		}
	}
		
		
	for(int Egbin=0;Egbin<nbBinsInEg;Egbin++){
		for(int mass=0;mass<nbBinsInM;mass++){
			for(int tbinning=0;tbinning<nbBinsInT;tbinning++){
				cout<<"volume 1 "<<Egbin<<" "<<mass<<" "<<tbinning<<" "<<volume1[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]<<" "<<endl;
				cout<<"volume 2 "<<Egbin<<" "<<mass<<" "<<tbinning<<" "<<volume2[tbinning+nbBinsInT*mass+nbBinsInT*nbBinsInM*Egbin]<<" "<<endl;
			}
		}
	}
		
	//cout<<"Load acceptance ..."<<argv[4]<<" Finished "<<endl;
		Acc[4]->SaveAs("checkAcc.root");
cout<<"good so far"<<endl;

	cout<<argc<<endl;

	//TChain *tc = new TChain("hipo2root");
	for(Int_t i=0;i<argc;i++){

		if(TString(argv[i]).Contains("MC")){

			type="MC";
			MC=true;
		}
		//if(TString(argv[4]).Contains(".root")){cout<<TString(argv[4]).Contains(".root")<<endl;}

		if(TString(argv[i]).Contains(".hipo")){

			nbf++;
			//cout<<"there"<<endl;
			nameFiles = TString(argv[i]);

		}

		


	}



	//hipo reader
	hipo::reader  reader;
	reader.open(nameFiles);

	hipo::dictionary  factory;
	reader.readDictionary(factory);
	factory.show();
	hipo::event      event;

	//hipo writer

	//hipo::writer  writer;
	//writer.open("essaiwriter.hipo");
/*
	TFile *limitFile = new TFile(TString(argv[5]));
	TGraphErrors *MaxAcc0 = (TGraphErrors*)limitFile->Get("MaxAcc0");                            
	TGraphErrors *MinAcc0 = (TGraphErrors*)limitFile->Get("MinAcc0");   

	TGraphErrors *MaxAcc1 = (TGraphErrors*)limitFile->Get("MaxAcc1");                            
	TGraphErrors *MinAcc1 = (TGraphErrors*)limitFile->Get("MinAcc1");   

	TGraphErrors *MaxAcc2 = (TGraphErrors*)limitFile->Get("MaxAcc2");                            
	TGraphErrors *MinAcc2 = (TGraphErrors*)limitFile->Get("MinAcc2");   

	TGraphErrors *MaxAcc3 = (TGraphErrors*)limitFile->Get("MaxAcc3");                            
	TGraphErrors *MinAcc3 = (TGraphErrors*)limitFile->Get("MinAcc3");
	    
	MaxAcc0->SetLineWidth(2);        
	MinAcc0->SetLineWidth(2);    
	MaxAcc1->SetLineWidth(2);                          
	MinAcc1->SetLineWidth(2);     
	MaxAcc2->SetLineWidth(2);                      
	MinAcc2->SetLineWidth(2);    
	MaxAcc3->SetLineWidth(2);                               
	MinAcc3->SetLineWidth(2);                                                         

	cout<<TString(argv[5])<<endl;

	new TCanvas; MaxAcc0->Draw();                      
*/
	TFile* outFile=new TFile("output1TCS.root","recreate");		

	int nbtc=0;//(tc->GetEntries());
	int nbJPSI=0;

	//tc->GetEntry(0);
	//TTree *tree = (TTree*)tc->GetTree()->CloneTree(0);


	//TMVA PID for Positron
	
	TMVA::Reader *readerTMVA = new TMVA::Reader( "!Color:!Silent" );

	// Create a set of variables and declare them to the reader
	// - the variable names MUST corresponds in name and type to those given in the weight file(s)
	Float_t SFPCAL, SFECIN, SFECOUT;
   	Float_t m2PCAL, m2ECIN, m2ECOUT;
   	readerTMVA->AddVariable( "SFPCAL", &SFPCAL );
   	readerTMVA->AddVariable( "SFECIN", &SFECIN );
   	readerTMVA->AddVariable( "SFECOUT", &SFECOUT );
   	readerTMVA->AddVariable( "m2PCAL", &m2PCAL );
   	readerTMVA->AddVariable( "m2ECIN", &m2ECIN );
   	readerTMVA->AddVariable( "m2ECOUT", &m2ECOUT );
	//Book TMVA method
	TString methodName2 = "MLP method" ;
	TString weightfile2 = "TMVAClassification_MLP6D.weights.xml";
	readerTMVA->BookMVA( methodName2, weightfile2 );




	//Angular Xsection assymetry
	TGraphErrors* Assym = new TGraphErrors();
	Assym->SetTitle("Assymetry in t bins");
	Assym->SetName("Assym");
	Assym->SetMarkerColor(4);
	Assym->SetMarkerStyle(21);

	//Final events
	
	TH1F *FinalEventst = new TH1F("FinalEventst","",50,0,1);
	FinalEventst->GetXaxis()->SetTitle("t (GeV^{2})");
	TH1F *FinalEventsEg = new TH1F("FinalEventsEg","",50,0,11);
	FinalEventsEg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
	TH1F *FinalEventsM = new TH1F("FinalEventsM","",50,1.,3.);
	FinalEventsM->GetXaxis()->SetTitle("M (GeV)");
	
	TH1F *FinalEventsPhi = new TH1F("FinalEventsPhi","",50,-180,180);
	FinalEventsPhi->GetXaxis()->SetTitle("Phi ");
	
	TH1F *FinalEventsP = new TH1F("FinalEventsP","",50,0,1.2);
	FinalEventsP->GetXaxis()->SetTitle("P ");
	
	TH1F *FinalEventsTheta = new TH1F("FinalEventsTheta","",50,5,60);
	FinalEventsTheta->GetXaxis()->SetTitle("Theta");
	
	//Beam spin assymetry
	TGraphErrors *BAssym = new TGraphErrors();
	BAssym->SetTitle("Beam Spin Assymetry");
	BAssym->SetMarkerColor(4);
	BAssym->SetMarkerStyle(21);

	cout<<"there"<<endl;


	TH2F *xihist = new TH2F("xihist","Eg(4-10 GeV),t[0.15-0.8];#xi;Q2",100,0.05,0.35,100,2.,9.);
	TH1F *xihist1 = new TH1F("xihist1","Eg(4-10 GeV),t[0.15-0.8];#xi;Q2",100,0.05,0.35);
	TH2F *thist = new TH2F("thist","Mass[1.5-2],Eg(4-10 GeV),t[0.15-0.8];#xi;-t",100,0.05,0.35,100,0.15,0.8);
	TH2F *thist1 = new TH2F("thist1","Mass[2-3],Eg(4-10 GeV),t[0.15-0.8];#xi;-t",100,0.05,0.35,100,0.15,0.8);

	TH2F* BeforeCuts= new TH2F("Before Cuts", " Before Cuts",200,-400,400,200,-400,400);
	TH2F* BeforeCutsPositron= new TH2F("Before Cuts Positron", " Before Cuts",200,-400,400,200,-400,400);

	TH2F* CheCooElec= new TH2F("CheCooElec", " CheCooElec",200,-110,110,200,-110,110);
	TH2F* CheCooPosi= new TH2F("CheCooElec Positron", " CheCoo Posi",200,-110,110,200,-110,110);

	TH2F* CheCooDiff= new TH2F("CheCooElec Diff Pi0", " CheCoo Posi",200,-110,110,200,-110,110);
	TH2F* CheCooDiffSup= new TH2F("CheCooElec Diff Sup", " CheCoo Posi",200,-110,110,200,-110,110);
	TH2F* CheCooDiffInf= new TH2F("CheCooElec Diff Inf", " CheCoo Posi",200,-110,110,200,-110,110);

	TH1F *histRECPid = new TH1F("RECpid","Histogram of rec pid",100,-15,2500);

	TH1F *histMMass = new TH1F("MMass","Missing mass beam",100,-20,20);
	histMMass->GetXaxis()->SetTitle("MM^{2} (GeV^{2})");
	//histMMass->GetYaxis()->SetTitle("number of events");		
	TH1F *histQ2 = new TH1F("Q2","Virtuality incoming photon (Q2)",100,-20,20);
	histQ2->GetXaxis()->SetTitle("Q^{2}  (GeV^{2})");
	TH1F *histQP2 = new TH1F("M2","Invariant Mass of the lepton pair squared",120,-1,5);	
	histQP2->GetXaxis()->SetTitle("M (GeV)");
	TH1F *histM = new TH1F("M","Invariant Mass of the lepton pair",200,0,4);
	histM->GetXaxis()->SetTitle("M (GeV)");
	TH2D *histME = new TH2D("MvsE","M vs Egamma",120,2,11,120,0,3);
	histME->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
	histME->GetYaxis()->SetTitle("Mass (GeV)");
	TH2D *histMPElec = new TH2D("MvsPElec","M vs PElec",120,0,11,120,0,3);
	histMPElec->GetXaxis()->SetTitle("Momentum Electron (GeV)");
	histMPElec->GetYaxis()->SetTitle("Mass (GeV)");
	TH2D *histMPPosi = new TH2D("MvsPPosi","M vs PPosi",120,0,11,120,0,3);
	histMPPosi->GetXaxis()->SetTitle("Momentum Positron (GeV)");
	histMPPosi->GetYaxis()->SetTitle("Mass (GeV)");
	TH2D *histMME = new TH2D("MMvsE","MM vs Egamma",100,2,11,100,-5,10);
	TH2D *histPMiss = new TH2D("PMiss","Direction of the missing particle",50,-0.5,0.5,50,-0.5,0.5);
	histPMiss->GetXaxis()->SetTitle("Px (%)");
	histPMiss->GetYaxis()->SetTitle("Py (%)");
	TH1F *histT =new TH1F("T rec","t REC",50,-0.2,2);
	histT->GetXaxis()->SetTitle("-t (GeV");
	TH2D *histMMPt = new TH2D("MMPt", "T. Momentum Fraction Vs MM beam", 100,-1,1,50,0,0.1);//-1,1,50,0,0.2);//
	histMMPt->GetXaxis()->SetTitle("Missing mass beam (GeV)");
	histMMPt->GetYaxis()->SetTitle("Transverse momentum fraction beam");
	TH1F *EnergyConservation = new TH1F("Cons. Energy","Cons. Energy",100,-1,1);	
	EnergyConservation->GetXaxis()->SetTitle("sum pz - energy incoming photon (GeV)");
	TH2D *MMQ2 = new TH2D("MMQ2"," Q2 VS Missing mass beam",50,-3,3,50,-3,3);
	MMQ2->GetXaxis()->SetTitle("Missing mass beam (GeV)");
	MMQ2->GetYaxis()->SetTitle("Virtuality incoming photon (GeV)");		

	TH2D *histQ2t = new TH2D("RECq2vst","Q2 vs -t (REC Particles)",50,4,9,30,0,1);
	TH2D *PhiVSPElectron = new TH2D("PhiVSPelectron", "P vs Phi (Electron)",100,-180,180,100,0,10);
	PhiVSPElectron->GetXaxis()->SetTitle("Phi (deg)");
	PhiVSPElectron->GetYaxis()->SetTitle("Momentum (GeV)");
	TH2D *PhiVSPPositron = new TH2D("PhiVSPPositron","P vs Phi (Positron)",100,-180,180,100,0,10);
	PhiVSPPositron->GetXaxis()->SetTitle("Phi (deg)");
	PhiVSPPositron->GetYaxis()->SetTitle("Momentum (GeV)");
	TH2D *PhiVSPProton = new TH2D("PhiVSPProton","P vs Phi (Proton)",100,-180,180,100,0,5);
	PhiVSPProton->GetXaxis()->SetTitle("Phi (deg)");
	PhiVSPProton->GetYaxis()->SetTitle("Momentum (GeV)");
	TH2D *ThetaVSPElectron = new TH2D("ThetaVSPelectron","Theta vs P (Electron)",100,0,10,100,0,50);
	ThetaVSPElectron->GetYaxis()->SetTitle("#theta (deg)");
	ThetaVSPElectron->GetXaxis()->SetTitle("Momentum (GeV)");
	TH2D *ThetaVSPPositron = new TH2D("ThetaVSPPositron","Theta vs P (Positron)",100,0,10,100,0,50);
	ThetaVSPPositron->GetYaxis()->SetTitle("#theta (deg)");
	ThetaVSPPositron->GetXaxis()->SetTitle("Momentum (GeV)");

	TH2D *histMPhiProton = new TH2D("MvsPhiProton","M vs PhiProton",120,-190,190,120,0,3);
	TH2D *histMPhiElectron = new TH2D("MvsPhiElectron","M vs PhiElec",120,-190,190,120,0,3);
	TH2D *histMPhiPositron = new TH2D("MvsPhiPositron","M vs PhiPosi",120,-190,190,120,0,3);
	TH2D *histMThetaElectron = new TH2D("MvsThetaElectrogdgdgn","M vs ThetaElec",100,0,50,120,0,3);
	TH2D *histAngleMass= new TH2D("MvsThetaElectron","M vs ThetaElec",100,0,0,120,0,3);
	TH2D *histMThetaPositron = new TH2D("MvsThetaPositron","M vs ThetaPosi",100,0,50,120,0,3);
	TH2D *histMThetaProton = new TH2D("MvsThetaProton","M vs ThetaProt",100,10,80,120,0,3);

	TH2D *histMVElectron = new TH2D("MvsVElectron","M vs Vx diff pair",100,-10,10,120,0,3);
	TH2D *histMVPositron = new TH2D("MvsVPositron","M vs Vy diff pair",100,-10,10,120,0,3);
	TH2D *histMV = new TH2D("MvsVProton","M vs Vz diff pair",100,-10,10,120,0,3);

	TH2D *histMCheElectron = new TH2D("MvsCheElectron","M vs CheElec",100,0,40,120,0,3);
	TH2D *histMChePositron = new TH2D("MvsChePositron","M vs ChePosi",100,0,40,120,0,3);
	TH2D *histMCheProton = new TH2D("MvsCheProton","M vs CheProt",100,0,40,120,0,3);

	TH2D *histCheCheElectron = new TH2D("MvsCheElectron","ChePosi vs CheElec",100,0,40,100,0,40);
	TH2D *histCheCheElectronInf = new TH2D("MvsCheElectrondgdgdd","ChePosi vs CheElec",100,0,40,100,0,40);
	TH2D *histCheCheElectronSup = new TH2D("MvsCheElectronaaaa","ChePosi vs CheElec",100,0,40,100,0,40);

	TH2D *PhiVSThetaElectron = new TH2D("PhiVSThetaelectron","Theta vs Phi (Electron)",100,-180,180,100,-0.1,50);
	PhiVSThetaElectron->GetXaxis()->SetTitle("Phi (deg)");
	PhiVSThetaElectron->GetYaxis()->SetTitle("Theta (deg)");
	TH2D *PhiVSThetaPositron = new TH2D("PhiVSThetaPositron","Theta vs Phi (Positron)",100,-180,180,100,-0.1,50);
	PhiVSThetaPositron->GetXaxis()->SetTitle("Phi (deg)");
	PhiVSThetaPositron->GetYaxis()->SetTitle("Theta (deg)");
	TH2D *PhiVSThetaProton = new TH2D("ThetaVSThetaProton","Theta vs Phi (Proton)",100,-180,180,100,-0.1,70);
	PhiVSThetaProton->GetXaxis()->SetTitle("Phi (deg)");
	PhiVSThetaProton->GetYaxis()->SetTitle("Theta (deg)");
	TH2D *ThetaVSPProton = new TH2D("ThetaVSPProton","Theta vs P (Proton)",100,0,5,100,0,70);
	ThetaVSPProton->GetYaxis()->SetTitle("#theta (deg)");
	ThetaVSPProton->GetXaxis()->SetTitle("Momentum (GeV)");

	TH1F *EhistMMass = new TH1F("EMMass","Missing mass beam",100,-1,1);
	EhistMMass->GetXaxis()->SetTitle("MM^{2} (GeV^{2})");
	TH1F *EhistQ2 = new TH1F("EQ2","Virtuality incoming photon (Q2)",80,-0.1,0.3);
	EhistQ2->GetXaxis()->SetTitle("Q^{2}  (GeV^{2})");
	TH1F *EhistQ21 = new TH1F("EQ21","Virtuality incoming photon (Q2)",80,-0.1,0.3);
	EhistQ21->GetXaxis()->SetTitle("Q^{2}  (GeV^{2})");		
	TH1F *EhistQP2 = new TH1F("EM","Invariant Mass of the lepton pair",120,-1,5);
	EhistQP2->GetXaxis()->SetTitle("M (GeV)");
	//TH1F *EhistM = new TH1F("EM1","",150,0.0,3.3);
	TH1F *EhistM = new TH1F("EM1","",200,0.0,4.0);
	EhistM->GetXaxis()->SetTitle("M (GeV)");

	TH1F *EhistMAfterFid = new TH1F("EM1AfterFid","Invariant Mass of the lepton pair AfterFid",200,0,4);
	EhistMAfterFid->GetXaxis()->SetTitle("M (GeV)");

	TH1F *EhistMprim = new TH1F("EM1prim","Invariant Mass of the lepton pair",200,0,3);
	EhistMprim->GetXaxis()->SetTitle("M (GeV)");

	TH1F *EhistEssai= new TH1F("EM2","Invariant Mass of the lepton pair with scattered electron",120,0,3);
	EhistEssai->GetXaxis()->SetTitle("M (GeV)");

	EhistMprim->SetLineColor(1);
	TH2D *EhistMSectorDiff = new TH2D("EM1SD","Sector Diff vs Invariant Mass of the lepton pair",120,0,3,6,-1,5);
	EhistMSectorDiff->GetXaxis()->SetTitle("M (GeV)");
	EhistMSectorDiff->GetYaxis()->SetTitle("Sector Diff)");
	TH2D *EhistME = new TH2D("EMvsE","M vs Egamma",100,2,11,100,0,2.5);
	EhistME->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
	EhistME->GetYaxis()->SetTitle("Mass (GeV)");
	TH2D *EhistMPElec = new TH2D("EMvsPElec","M vs P Electron",100,0,11,100,0,2.5);
	EhistMPElec->GetXaxis()->SetTitle("Momentum (GeV)");
	EhistMPElec->GetYaxis()->SetTitle("Mass (GeV)");
	TH2D *EhistMPPosi = new TH2D("EMvsPPosi","M vs P Positron",100,0,11,100,0,2.5);
	EhistMPPosi->GetXaxis()->SetTitle("Momentum (GeV)");
	EhistMPPosi->GetYaxis()->SetTitle("Mass (GeV)");
	TH2D *EhistMME = new TH2D("EMMvsE","MM vs Egamma",100,2,11,100,-5,10);
	TH2D *EhistMPProton = new TH2D("EMvsPProton","M vs P Proton",100,0,2,100,0,2.5);
	EhistMPProton->GetXaxis()->SetTitle("Momentum (GeV)");
	EhistMPProton->GetYaxis()->SetTitle("Mass (GeV)");
	TH2D *EhistMXProton = new TH2D("EMvsXProton","M vs XProton",100,-10,10,100,0,2.5);
	TH2D *EhistMRun = new TH2D("EMvsRun","M vs RUN",65,4013,4078,30,0,2.5);		
	EhistMRun->GetXaxis()->SetTitle("Run Num.");
	EhistMRun->GetYaxis()->SetTitle("Mass (GeV)");
	TH2D *EhistMPhiProton = new TH2D("EMvsPhiProton","M vs PhiProton",120,-190,190,120,0,3);
	TH2D *EhistPMiss = new TH2D("EPMiss","Direction of the missing particle",50,-0.5,0.5,50,-0.5,0.5);
	EhistPMiss->GetXaxis()->SetTitle("Px (%)");
	EhistPMiss->GetYaxis()->SetTitle("Py (%)");
	TH1F *EhistT =new TH1F("ET rec","t REC",50,-0.2,2);
	EhistT->GetXaxis()->SetTitle("-t (GeV");	
	TH2D *EhistMMPt = new TH2D("EMMPt", "Pt Vs MM", 100,-5,5,50,0,1);
	EhistMMPt->GetXaxis()->SetTitle("Missing mass beam (GeV)");
	EhistMMPt->GetYaxis()->SetTitle("Transverse momentum fraction beam");
	TH1F *EEnergyConservation = new TH1F("ECons. Energy","Cons. Energy",100,-1,1);
	EEnergyConservation->GetXaxis()->SetTitle("#sum P_{z} - E_{#gamma} (GeV)");	
	TH1F *EEnergyConservationEvents = new TH1F("ECons. EnergyEvents","Cons. Energy",100,-1,1);
	EEnergyConservationEvents->GetXaxis()->SetTitle("#sum P_{z} - E_{#gamma} (GeV)");		
	TH2D *EMMQ2 = new TH2D("EMMQ2","MM Q2",50,-15,15,50,-6,6);	
	EMMQ2->GetXaxis()->SetTitle("Missing mass beam (GeV)");
	EMMQ2->GetYaxis()->SetTitle("Virtuality incoming photon (GeV)");

	TH2D *EPhiVSPElectron = new TH2D("EPhiVSPelectron", "P vs Phi (Electron)",100,-180,180,100,0,10);
	EPhiVSPElectron->GetXaxis()->SetTitle("Phi (deg)");
	EPhiVSPElectron->GetYaxis()->SetTitle("Momentum (GeV)");
	TH2D *EPhiVSPPositron = new TH2D("EPhiVSPPositron","P vs Phi (Positron)",100,-180,180,100,0,10);
	EPhiVSPPositron->GetXaxis()->SetTitle("Phi (deg)");
	EPhiVSPPositron->GetYaxis()->SetTitle("Momentum (GeV)");
	TH2D *EPhiVSPProton = new TH2D("EPhiVSPProton", "P vs Phi (Proton)",100,-180,180,100,0,5);
	EPhiVSPProton->GetXaxis()->SetTitle("Phi (deg)");
	EPhiVSPProton->GetYaxis()->SetTitle("Momentum (GeV)");
	TH2D *EThetaVSPElectron = new TH2D("EThetaVSPelectron","Theta vs P (Electron)",100,0,10,100,0,50);
	EThetaVSPElectron->GetYaxis()->SetTitle("#theta (deg)");
	EThetaVSPElectron->GetXaxis()->SetTitle("Momentum (GeV)");	
	TH2D *EThetaVSPPositron = new TH2D("EThetaVSPPositron","Theta vs P (Positron)",100,0,10,100,0,50);
	EThetaVSPPositron->GetYaxis()->SetTitle("#theta (deg)");
	EThetaVSPPositron->GetXaxis()->SetTitle("Momentum (GeV)");

	TH2D *EPhiVSThetaElectron = new TH2D("EPhiVSThetaelectron","Theta vs Phi (Electron)",100,-180,180,100,-0.1,50);
	EPhiVSThetaElectron->GetXaxis()->SetTitle("Phi (deg)");
	EPhiVSThetaElectron->GetYaxis()->SetTitle("Theta (deg)");
	TH2D *EPhiVSThetaPositron = new TH2D("EPhiVSThetaPositron","Theta vs Phi (Positron)",100,-180,180,100,-0.1,50);
	EPhiVSThetaPositron->GetXaxis()->SetTitle("Phi (deg)");
	EPhiVSThetaPositron->GetYaxis()->SetTitle("Theta (deg)");		
	TH2D *EPhiVSThetaProton = new TH2D("EThetaVSThetaProton","Theta vs Phi (Proton)",100,-180,180,100,-0.1,70);
	EPhiVSThetaProton->GetXaxis()->SetTitle("Phi (deg)");
	EPhiVSThetaProton->GetYaxis()->SetTitle("Theta (deg)");
	TH2D *EThetaVSPProton = new TH2D("EThetaVSPProton","Theta vs P (Proton)",100,0,5,50,0,70);
	EThetaVSPProton->GetYaxis()->SetTitle("#theta (deg)");
	EThetaVSPProton->GetXaxis()->SetTitle("Momentum (GeV)");

	TH2D *PhiVSThetaCM = new TH2D("PhiVSThetaCMm","Theta vs Phi (CM)",binPhi,minphi,maxphi,binTheta,mintheta,maxtheta);
	TH2D *PhiVSThetaCM1 = new TH2D("PhiVSThetaCM1m","Theta vs Phi (CM)",binPhi,minphi,maxphi,binTheta,mintheta,maxtheta);
	TH2D *PhiVSThetaCM2 = new TH2D("PhiVSThetaCM2m","Theta vs Phi (CM)",binPhi,minphi,maxphi,binTheta,mintheta,maxtheta);
	TH2D *PhiVSThetaCM3 = new TH2D("PhiVSThetaCM3m","Theta vs Phi (CM)",binPhi,minphi,maxphi,binTheta,mintheta,maxtheta);

	TH1F *vertexPair = new TH1F("vertexPair","Vertex diff pair",50,-10,10);
	TH1F *vertexElecP = new TH1F("vertexElecP","Vertex Elec p",50,-10,10);
	TH1F *vertexPosiP = new TH1F("vertexPosP","Vertex Posi p",50,-10,10);

	TH1F *vertexElec = new TH1F("vertexElec","Vertex elec",50,-10,10);
	TH1F *vertexPosi = new TH1F("vertexPosi","Vertex posi ",50,-10,10);
	TH1F *vertexProt = new TH1F("vertexProt","Vertex proton",50,-10,10);

	TH2F *vertexProtTheta = new TH2F("vertexProtTheta","Vertex proton vs Theta",100,0,70,50,-10,10);

	TH2F *vertexProtThetaCD = new TH2F("vertexProtThetaCD","Vertex proton vs Theta CD",100,0,70,50,-10,10);
	TH2F *vertexProtThetaFD = new TH2F("vertexProtThetaFD","Vertex proton vs Theta FD",100,0,70,50,-10,10);

	TH1F *EvertexPair = new TH1F("EvertexPair","Vertex diff pair",50,-10,10);
	TH1F *EvertexElecP = new TH1F("EvertexElecP","Vertex Elec p",50,-10,10);
	TH1F *EvertexPosiP = new TH1F("EvertexPosP","Vertex Posi p",50,-10,10);

	TH1F *vertexTimePair =new TH1F("vertexTimePair","vertex time pair",150,-2,2);
	vertexTimePair->GetXaxis()->SetTitle("Vertex time difference (ns)");
	TH1F *vertexTimeP =new TH1F("vertexTimeP","vertex time proton",150,-2,2);
	vertexTimeP->GetXaxis()->SetTitle("Vertex time difference (ns)");
	TH2D *vertexTimePP =new TH2D("vertexTimePP","vertex time proton vs p proton",100,0,3,100,-2,2);
	vertexTimePP->GetYaxis()->SetTitle("Vertex time difference (ns)");
	vertexTimePP->GetXaxis()->SetTitle("Momentum (GeV)");

	TH1F *Pproton =new TH1F("Pproton","Pproton",100,0.,100);
	Pproton->GetXaxis()->SetTitle("P");


	TH1F *EvertexTimePair =new TH1F("EvertexTimePair","Vertex time pair",150,-2,2);
	EvertexTimePair->GetXaxis()->SetTitle("Vertex time difference (ns)");
	TH2F *vertexTimePairMass =new TH2F("vertexTimePairMass","",100,0,3,150,-2,2);
	
	TH1F *EvertexTimeP =new TH1F("EvertexTimeP","vertex time proton",100,-0.5,0.5);
	EvertexTimeP->GetXaxis()->SetTitle("Vertex time difference (ns)");
	TH2D *EvertexTimePPCD =new TH2D("EvertexTimePP","vertex time proton vs p proton CD",60,0.25,1,100,-0.5,0.5);
	EvertexTimePPCD->GetYaxis()->SetTitle("Vertex time difference (ns)");
	EvertexTimePPCD->GetXaxis()->SetTitle("Momentum (GeV)");
	TH2D *EvertexTimePPFD =new TH2D("EvertexTimePPFD","vertex time proton vs p proton FD",60,0.25,1,100,-0.5,0.5);
	EvertexTimePPFD->GetYaxis()->SetTitle("Vertex time difference (ns)");
	EvertexTimePPFD->GetXaxis()->SetTitle("Momentum (GeV)");

	//TCS distribution
	TH2D *TCSEThetaVSPElectron = new TH2D("EThetaVSPelectron","Theta vs P (Electron)",100,0,10,100,0,50);
	TCSEThetaVSPElectron->GetYaxis()->SetTitle("#theta (deg)");
	TCSEThetaVSPElectron->GetXaxis()->SetTitle("Momentum (GeV)");
	TH2D *TCSEThetaVSPPositron = new TH2D("EThetaVSPPositron","Theta vs P (Positron)",100,0,10,100,0,50);
	TCSEThetaVSPPositron->GetYaxis()->SetTitle("#theta (deg)");
	TCSEThetaVSPPositron->GetXaxis()->SetTitle("Momentum (GeV)");

	TH2D *TCSEThetaVSPProton = new TH2D("EThetaVSPProton","Theta vs P (Proton)",50,0,1.5,50,0,70);
	TCSEThetaVSPProton->GetYaxis()->SetTitle("#theta (deg)");
	TCSEThetaVSPProton->GetXaxis()->SetTitle("Momentum (GeV)");

	TH2D *PPosiPElec = new TH2D("PPosiPElec","All events;P (e+);P (e-)",50,0,10,50,0,10);
	TH2D *PPosiPElec1 = new TH2D("PPosiPElec1","M>1.5 GeV;P (e+);P (e-)",50,0,10,50,0,10);
	//Beam spin assymetry
	TH1F *HPosi =new TH1F("HPosi","HPosi",10,-180,180);
	TH1F *HNega =new TH1F("HNega","HNega",10,-180,180);

	//PID histogram
	TH2D *SFelectron = new TH2D("SFelectron","SF vs P (Electron)",80,0,8,100,0.,0.35);
	SFelectron->GetXaxis()->SetTitle("Momentum (Gev)");
	SFelectron->GetYaxis()->SetTitle("Sampling fraction");
	TH2D *SFpositron = new TH2D("SFpositron","SF vs P (Positron)",80,0,8,100,0.,0.35);
	SFpositron->GetXaxis()->SetTitle("Momentum (Gev)");
	SFpositron->GetYaxis()->SetTitle("Sampling fraction");

	TH2D *corrSFelectron = new TH2D("corrSFelectron","corrSF vs P (Electron)",80,0,8,100,0.1,0.35);
	TH2D *corrSFpositron = new TH2D("corrSFpositron","corrSF vs P (Positron)",80,0,8,100,0.1,0.35);

	TH2D *ECelectron = new TH2D("ECelectron","ECout vs ECin (Electron)",50,0,0.7,50,0,0.7);
	TH2D *ECpositron = new TH2D("ECpositron","ECout vs ECin (Positron)",50,0,0.7,50,0,0.7);
	TH2D *EECelectron = new TH2D("EECelectron","ECout vs ECin (Electron)",50,0,0.7,50,0,0.7);
	TH2D *EECpositron = new TH2D("EECpositron","ECout vs ECin (Positron)",50,0,0.7,50,0,0.7);

	TH1F *CheElectron = new TH1F("Cherenkov Electron", "Cherenkov Electron",40,0,40);
	CheElectron->GetXaxis()->SetTitle("number of photon");
	TH1F *ChePositron = new TH1F("Cherenkov Positron", "Cherenkov Positron",40,0,40);
	ChePositron->GetXaxis()->SetTitle("number of photon");

	TH1F *ECheElectron = new TH1F("ECherenkov Electron", "Cherenkov Electron",40,0,40);
	ECheElectron->GetXaxis()->SetTitle("number of photon");
	TH1F *EChePositron = new TH1F("ECherenkov Positron", "Cherenkov Positron",40,0,40);
	EChePositron->GetXaxis()->SetTitle("number of photon");		

	TH2D *BetaProton = new TH2D("Beta P Proton"," Beta P Proton ",100,0,2.5,50,0,1.2);
	BetaProton->GetXaxis()->SetTitle("Momentum (GeV)");
	BetaProton->GetYaxis()->SetTitle("#beta");
	TH2D *EBetaProton = new TH2D("EBeta P Proton"," Beta P Proton ",100,0,2.5,50,0,1.2);
	EBetaProton->GetXaxis()->SetTitle("Momentum (GeV)");
	EBetaProton->GetYaxis()->SetTitle("#beta");

	TH2D *SFUelectron = new TH2D(" SF Electron vs U", " SF Electron vs U",150,0,450,50,0.,0.35);
	SFUelectron->GetXaxis()->SetTitle("Position (cm)");
	SFUelectron->GetYaxis()->SetTitle("Sampling Fraction Electron");
	TH2D *SFVelectron = new TH2D(" SF Electron vs V", " SF Electron vs V",150,0,450,50,0.,0.35);
	SFVelectron->GetXaxis()->SetTitle("Position (cm)");
	SFVelectron->GetYaxis()->SetTitle("Sampling Fraction Electron");
	TH2D *SFWelectron = new TH2D(" SF Electron vs W", " SF Electron vs W",150,0,450,50,0.,0.35);
	SFWelectron->GetXaxis()->SetTitle("Position (cm)");
	SFWelectron->GetYaxis()->SetTitle("Sampling Fraction Electron");

	TH2D *SFUpositron = new TH2D(" SF Positron vs U", " SF Positron vs U",150,0,450,50,0.,0.35);
	SFUpositron->GetXaxis()->SetTitle("Position (cm)");
	SFUpositron->GetYaxis()->SetTitle("Sampling Fraction Positron");
	TH2D *SFVpositron = new TH2D(" SF Positron vs V", " SF Positron vs V",150,0,450,50,0.,0.35);
	SFVpositron->GetXaxis()->SetTitle("Position (cm)");
	SFVpositron->GetYaxis()->SetTitle("Sampling Fraction Positron");
	TH2D *SFWpositron = new TH2D(" SF Positron vs W", " SF Positron vs W",150,0,450,50,0.,0.35);
	SFWpositron->GetXaxis()->SetTitle("Position (cm)");
	SFWpositron->GetYaxis()->SetTitle("Sampling Fraction Positron");

	TH2D *FiducialCutPCAL = new TH2D(" FiducialPCAL", "FiducialPCAL",500,-500,500,500,-500,500);
	TH2D *EFiducialCutPCAL = new TH2D(" EFiducialPCAL", "FiducialPCAL",500,-500,500,500,-500,500);

	TH1D *VertexOtherPart = new TH1D("VertexOtherPart", "VertexOtherPart",500,-30,30);
	TH1D *Chi2OtherPart = new TH1D("Chi2OtherPart", "Chi2OtherPart",500,-40,40);

	TH2F *EgVSPeletron = new TH2F("EgVSPeletron","EgVSPeletron",100,0,10,100,3,10);
	TH2F *tVSPeletron = new TH2F("tVSPeletron","tVSPeletron",100,0,10,100,0.0,1.);
	//Pio Study
	TH1F *FuckingPi0 = new TH1F("Pi0","Pi0",50,0.0,0.3);	
	TH2F *VertexPio = new TH2F("Pi0vertex","Pi0vertex",50,-7,3,50,-7,3);
	TH2F *VertexAutre = new TH2F("vertex autre","vertex autre",50,-7,3,50,-7,3);	
	TH2F *VertexAutreInf = new TH2F("veddddrtex autre","vertex autre",50,-7,3,50,-7,3);	

	TH2D *MomentumPairPio = new TH2D(" MomentumPairPio", "MomentumPairPio",100,0,10,100,0,10);
	TH2D *MomentumPairAutre = new TH2D(" MomentumPairAutre", "MomentumPairAutre",100,0,10,100,0,10);
	TH2D *MomentumPairAutreInf = new TH2D(" MomentumPaddddirAutre", "MomentumPairAutre",100,0,10,100,0,10);

	TH2D *ThetaPairPio = new TH2D("ThetaPairPio", "ThetaPairPio",100,0,40,100,0,40);
	TH2D *ThetaPairAutre = new TH2D("ThetaPairAutre", "ThetaPairAutre",100,0,40,100,0,40);
	TH2D *ThetaPairAutreInf = new TH2D("ThetaPairdddAutre", "ThetaPairAutre",100,0,40,100,0,40);

	TH2D *TimeChePio = new TH2D("TimeChePio", "TimeChePio",100,50,180,100,50,180);
	TH2D *TimeChePioSup = new TH2D("TimeChePioSup", "TimeChePioSup",100,50,180,100,50,180);
	TH2D *TimeChePioInf = new TH2D("TimeChePioInf", "TimeChePioInf",100,50,180,100,50,180);
	TH1D *TimeCheDiff = new TH1D("TimeCheDiff", "TimeCheDiff",100,-10,10);
	TH2D *TimeCheDiffMass = new TH2D("TimeCheDiffMass", "TimeCheDiffMass",100,0,3,100,-20,20);

	TH2D *SFMassElec = new TH2D("SFMass", "SFMassElec",100,0,3,100,0.,0.3);
	TH2D *SFMassPosi = new TH2D("SFMass", "SFMassPosi",100,0,3,100,0.,0.3);

	TH1F *tCheck = new TH1F("tCheck","",100,-1,1);
	tCheck->GetXaxis()->SetTitle("#Delta t");
	
	TH1F *tCheckEvents = new TH1F("tCheckEvents","",100,-1,1);
	tCheckEvents->GetXaxis()->SetTitle("#Delta t");
	
	TH1F *tCheckBefore = new TH1F("tCheckBefore","",100,-1,1);
	tCheckBefore->GetXaxis()->SetTitle("#Delta t");

	//Check Posi PCAl/Ecal sf
	TH2D *SFPCALECALPosiavant = new TH2D("SFPCALECALPosi", "SFPCALECALPosi",100,0,0.3,100,0.,0.3);
	TH2D *SFPCALECALPosiapres = new TH2D("SFPCALECALPosia", "SFPCALECALPosia",100,0,0.3,100,0.,0.3);

	TH2D *SFPCALECALElecavant = new TH2D("SFPCALECALElec", "SFPCALECALElec",100,0,0.3,100,0.,0.3);
	TH2D *SFPCALECALElecapres = new TH2D("SFPCALECALEleca", "SFPCALECALEleca",100,0,0.3,100,0.,0.3);

	TH2D *SFPpcalzeroPosi= new TH2D("SFPpcalzeroPosi", "SFPpcalzeroPosi",100,0,11,100,0.,0.3);
	TH2D *SFPpcalzeroElec= new TH2D("SFPpcalzeroElec", "SFPpcalzeroElec",100,0,11,100,0.,0.3);
	
	TH2F *Q2t= new TH2F("tmass",";t;mass",100,0.15,0.8,100,2,9.);
	
	//Missing particle
	
	TH2F *PhiMissingParticle= new TH2F("PhiMissingParticle",";phi;Pt/P",100,-180.,180,100,0,0.7);
	
	//chi2 Part
	
	TH1F *Chi2ElectronAvant= new TH1F("Chi2ElectronAvant","Electron (Before exclu.); #chi^{2}",100,-6,6);
	TH1F *Chi2PositronAvant= new TH1F("Chi2PositronAvant","Positron (Before exclu.); #chi^{2}",100,-6,6);
	TH1F *Chi2ProtonAvant= new TH1F("Chi2ProtonAvant","Proton CD (Before exclu.); #chi^{2}",100,-10,10);
	TH1F *Chi2ProtonFDAvant= new TH1F("Chi2ProtonFDAvant","Proton FD (Before exclu.); #chi^{2}",100,-10,10);
	
	
	TH1F *Chi2Electron= new TH1F("Chi2Electron","Electron; #chi^{2}",100,-6,6);
	TH1F *Chi2Positron= new TH1F("Chi2Positron","Positron; #chi^{2}",100,-6,6);
	TH1F *Chi2Proton= new TH1F("Chi2Proton","Proton CD; #chi^{2}",100,-10,10);
	TH1F *Chi2ProtonFD= new TH1F("Chi2ProtonFD","Proton FD; #chi^{2}",100,-10,10);
	
	// GDA check
	
	TH1F *Udistribution= new TH1F("Udistribution",";u",100,0.,13.);
	TH2F *UdistributionVsTheta= new TH2F("UdistributionVsTheta",";u;theta",100,0.,1.5,100,0.,60);
	
	//Radiative corrections
	
	TH1F *diffPhidiffThetaRad= new TH1F("diffPhidiffThetaRad",";Cone Angle",100,0,180);//,100,0,10);
	TH1F *diffPTRad= new TH1F("diffPTRad","diffPTRad",100,-10,00);//,100,0,10);
	TH1F *diffPTNORad= new TH1F("diffPTRad","diffPTRad",100,-10,00);//,100,0,10);
	TH1F *diffPhidiffThetaRad1= new TH1F("diffPhidiffThetaRad1","diffPhidiffThetaRad1",100,-180,180);//,100,0,10);
	TH2F *diffPhidiffThetaRad1vsTheta= new TH2F("diffPhidiffThetaRad1vsTheta","diffPhidiffThetaRad1vsTheta",100,-10,10,100,5,40);
	TH2F *diffPhidiffThetaRad1vsP= new TH2F("diffPhidiffThetaRad1vsP",";#Delta #theta;P",100,-10,10,100,0.5,10);
	TH2F *diffPhidiffThetaRad1vsPhi= new TH2F("diffPhidiffThetaRad1vsPhi","diffPhidiffThetaRad1vsPhi",100,-10,10,100,-120,120);
	TH1F *DiffRadCor= new TH1F("DiffRadCor",";#Delta P/P",100,0,1);//,100,0,10);

	//Plot binning
	TH2F *EgVST = new TH2F("EgVST","EgVST",80,3,11,80,0.05,1.1);
	TH2F *MVST = new TH2F("MVST","MVST",80,1.4,3.1,80,0.05,1.1);
	TH2F *MVSEg = new TH2F("MVST","MVST",80,1.4,3.1,80,2,11);

	TH1F *Positronpt= new TH1F("PositronPt","PositronPt",100,0.,4.);

	TH1F *pos= new TH1F("pos","",4,t1);
	TH1F *neg= new TH1F("neg","neg",4,t1);
	pos->Sumw2();
	neg->Sumw2();
	TH1F *tmean = new TH1F("tmean","tmean",4,0,4);
	TH1F *tmeanN = new TH1F("tmeanN","tmeanN",4,0,4);
	
	//for the exclu cut syst
	TH1F *posExclusyst= new TH1F("posExclusyst","",1,0.15,0.8);
	TH1F *negExclusyst= new TH1F("negExclusyst","negExclusyst",1,0.15,0.8);
	posExclusyst->Sumw2();
	negExclusyst->Sumw2();
	
	//AFB in mass bins
	
	TH1F *pos1= new TH1F("pos1","",4,t1);
	TH1F *neg1= new TH1F("neg1","neg1",4,t1);
	pos1->Sumw2();
	neg1->Sumw2();
	TH1F *tmean1 = new TH1F("tmean","tmean",4,0,4);
	TH1F *tmeanN1 = new TH1F("tmeanN","tmeanN",4,0,4);
	
	TH1F *pos2= new TH1F("pos2","",4,t1b);
	TH1F *neg2= new TH1F("neg2","neg2",4,t1b);
	pos2->Sumw2();
	neg2->Sumw2();
	TH1F *tmean2 = new TH1F("tmean2","tmean2",4,0,4);
	TH1F *tmeanN2 = new TH1F("tmeanN2","tmeanN2",4,0,4);
	//for the exclu cut syst
	TH1F *pos2Exclusyst= new TH1F("pos2Exclusyst","",1,0.15,0.8);
	TH1F *neg2Exclusyst= new TH1F("neg2Exclusyst","neg2Exclusyst",1,0.15,0.8);
	pos2Exclusyst->Sumw2();
	neg2Exclusyst->Sumw2();
	
	//AFB in high mass
	
	TH1F *postHighMass= new TH1F("postHighMass","",4,t1);
	TH1F *negtHighMass= new TH1F("negtHighMass","negt",4,t1);
	postHighMass->Sumw2();
	negtHighMass->Sumw2();
	TH1F *tmeanHighMass = new TH1F("tmeanHighMass","tmean",4,0,4);
	TH1F *tmeanNHighMass = new TH1F("tmeanNHighMass","tmeanN",4,0,4);
	
	TH1F *posXiHighMass= new TH1F("posXiHighMass","",3,Xi3);
	TH1F *negXiHighMass= new TH1F("negXiHighMass","",3,Xi3);
	posXiHighMass->Sumw2();
	negXiHighMass->Sumw2();
	TH1F *XimeanHighMass = new TH1F("XimeanHighMass","Ximean",3,0,3);
	TH1F *XimeanNHighMass = new TH1F("XimeanHighMass","XimeanN",3,0,3);
	
	TH1F *posEgHighMass= new TH1F("posEgHighMass","",3,4,10);
	TH1F *negEgHighMass= new TH1F("negEgHighMass","",3,4,10);
	posEgHighMass->Sumw2();
	negEgHighMass->Sumw2();
	TH1F *EgmeanHighMass = new TH1F("EgmeanHighMass","Egmean",3,0,3);
	TH1F *EgmeanNHighMass = new TH1F("EgmeanNHighMass","EgmeanN",3,0,3);
	
	
	////////////
	TH1F *posXi= new TH1F("posXi","",3,Xi3);
	TH1F *negXi= new TH1F("negXi","",3,Xi3);
	posXi->Sumw2();
	negXi->Sumw2();
	TH1F *Ximean = new TH1F("Ximean","Ximean",3,0,3);
	TH1F *XimeanN = new TH1F("Ximean","XimeanN",3,0,3);
	
	TH1F *posEg= new TH1F("posEg","",3,4,10);
	TH1F *negEg= new TH1F("negEg","",3,4,10);
	posEg->Sumw2();
	negEg->Sumw2();
	TH1F *Egmean = new TH1F("Egmean","Egmean",3,0,3);
	TH1F *EgmeanN = new TH1F("EgmeanN","EgmeanN",3,0,3);
	
	TH1F *posQ2= new TH1F("posQ2","",4,M1);
	TH1F *negQ2= new TH1F("negQ2","",4,M1);
	posQ2->Sumw2();
	negQ2->Sumw2();
	TH1F *Q2mean = new TH1F("Q2mean","Q2mean",4,0,4);
	TH1F *Q2meanN = new TH1F("Q2mean","Q2meanN",4,0,4);
	
	TH1F *posQ2Phi= new TH1F("posQ2Phi","",2,1.00,1.04);
	TH1F *negQ2Phi= new TH1F("negQ2Phi","",2,1.00,1.04);
	posQ2Phi->Sumw2();
	negQ2Phi->Sumw2();
	
	
	TH1F *pospola= new TH1F("pospola","",10,-180,180);
	TH1F *negpola= new TH1F("negpola","",10,-180,180);
	pospola->Sumw2();
	negpola->Sumw2();
	
	
	TH1F *RratioNum[4];
	TH1F *RratioDenom[4];
	for(int tbin=0;tbin<4;tbin++){
		RratioNum[tbin] = new TH1F(Form("t bin %d num",tbin), Form("t bin %d num",tbin) ,10,-180,180);
		RratioDenom[tbin] = new TH1F(Form("t bin %d denom",tbin), Form("t bin %d denom",tbin) ,10,-180,180);
		RratioNum[tbin]->Sumw2();
		RratioDenom[tbin]->Sumw2();
	}
	TH1F *tmeanRratio = new TH1F("tmeanRratio","tmeanRratio",4,0,4);
	TH1F *tmeanRratioN = new TH1F("tmeanRratioN","tmeanRratioN",4,0,4);
	
	TH1F *RratioNumHighMass[4];
	TH1F *RratioDenomHighMass[4];
	for(int tbin=0;tbin<4;tbin++){
		RratioNumHighMass[tbin] = new TH1F(Form("t bin %d numHighMass",tbin), Form("t bin %d num",tbin) ,10,-180,180);
		RratioDenomHighMass[tbin] = new TH1F(Form("t bin %d denomHighMass",tbin), Form("t bin %d denom",tbin) ,10,-180,180);
		RratioNumHighMass[tbin]->Sumw2();
		RratioDenomHighMass[tbin]->Sumw2();
	}
	TH1F *tmeanRratioHighMass = new TH1F("tmeanRratioHighMass","tmeanRratio",4,0,4);
	TH1F *tmeanRratioNHighMass = new TH1F("tmeanRratioNHighMass","tmeanRratioN",4,0,4);
	
	TH1F *RratioNumRafo[8];
	TH1F *RratioDenomRafo[8];
	for(int tbin=0;tbin<8;tbin++){
		RratioNumRafo[tbin] = new TH1F(Form("t bin %d num Rafo",tbin), Form("t bin %d num",tbin) ,10,-180,180);
		RratioDenomRafo[tbin] = new TH1F(Form("t bin %d denom Rafo",tbin), Form("t bin %d denom",tbin) ,10,-180,180);
		RratioNumRafo[tbin]->Sumw2();
		RratioDenomRafo[tbin]->Sumw2();
	}
	TH1F *tmeanRratioRafo = new TH1F("tmeanRratioRafo","tmeanRratio",8,0,8);
	TH1F *tmeanRratioNRafo = new TH1F("tmeanRratioNRafo","tmeanRratioN",8,0,8);
	Double_t t1Rafo [9] = {0.15,0.22, 0.29, 0.36, 0.43,0.50, 0.57,0.64,0.8};//{
	Double_t M1Rafo [5] = {2.25, 3.5 ,4.,5.,9.};
	TH2D *PhysicstM2forbinningRafo= new TH2D("PhysicstM2forbinningRafo","",8,t1Rafo,4,M1Rafo);
	
	TH1F *RratioNumXi[3];
	TH1F *RratioDenomXi[3];
	for(int Xibin=0;Xibin<3;Xibin++){
		RratioNumXi[Xibin] = new TH1F(Form("Xi bin %d num",Xibin), Form("Xi bin %d num",Xibin) ,10,-180,180);
		RratioDenomXi[Xibin] = new TH1F(Form("Xi bin %d denom",Xibin), Form("Xi bin %d denom",Xibin) ,10,-180,180);
		RratioNumXi[Xibin]->Sumw2();
		RratioDenomXi[Xibin]->Sumw2();
	}
	TH1F *XimeanRratio = new TH1F("XimeanRratio","XimeanRratio",3,0,3);
	TH1F *XimeanRratioN = new TH1F("XimeanRratioN","XimeanRratioN",3,0,3);
	
	TH1F *BSApos[4];
	TH1F *BSAneg[4];
	for(int tbin=0;tbin<4;tbin++){
		BSApos[tbin] = new TH1F(Form("t bin %d pos",tbin), Form("t bin %d pos",tbin) ,10,0,360);
		BSAneg[tbin] = new TH1F(Form("t bin %d neg",tbin), Form("t bin %d neg",tbin) ,10,0,360);
		BSApos[tbin]->Sumw2();
		BSAneg[tbin]->Sumw2();
	}
	
	//exclu syst
	TH1F *BSAposExclusyst[1];
	TH1F *BSAnegExclusyst[1];
	for(int tbin=0;tbin<1;tbin++){
		BSAposExclusyst[tbin] = new TH1F(Form("t bin %d pos Exclusyst",tbin), Form("t bin %d pos Exclusyst",tbin) ,10,0,360);
		BSAnegExclusyst[tbin] = new TH1F(Form("t bin %d neg Exclusyst",tbin), Form("t bin %d neg Exclusyst",tbin) ,10,0,360);
		BSAposExclusyst[tbin]->Sumw2();
		BSAnegExclusyst[tbin]->Sumw2();
	}
	
	
	TH1F *tmeanBSA = new TH1F("tmeanBSA","tmeanBSA",4,0,4);
	TH1F *tmeanBSAN = new TH1F("tmeanBSAN","tmeanBSAN",4,0,4);
	
	TH1F *Q2meanBSA = new TH1F("Q2meanBSA","Q2meanBSA",1,0,11);
	
	TH1F *BSAposXi[3];
	TH1F *BSAnegXi[3];
	for(int Xibin=0;Xibin<3;Xibin++){
		BSAposXi[Xibin] = new TH1F(Form("Xi bin %d pos",Xibin), Form("Xi bin %d pos",Xibin) ,10,-180,180);
		BSAnegXi[Xibin] = new TH1F(Form("Xi bin %d neg",Xibin), Form("Xi bin %d neg",Xibin) ,10,-180,180);
		BSAposXi[Xibin]->Sumw2();
		BSAnegXi[Xibin]->Sumw2();
	}
	TH1F *XimeanBSA = new TH1F("XimeanBSA","XimeanBSA",3,0,3);
	TH1F *XimeanBSAN = new TH1F("XimeanBSAN","XimeanBSAN",3,0,3);
	
	TH1F *BSAposM[4];
	TH1F *BSAnegM[4];
	for(int Mbin=0;Mbin<4;Mbin++){
		BSAposM[Mbin] = new TH1F(Form("M bin %d pos",Mbin), Form("M bin %d pos",Mbin) ,10,-180,180);
		BSAnegM[Mbin] = new TH1F(Form("M bin %d neg",Mbin), Form("M bin %d neg",Mbin) ,10,-180,180);
		BSAposM[Mbin]->Sumw2();
		BSAnegM[Mbin]->Sumw2();
	}
	TH1F *MmeanBSA = new TH1F("MmeanBSA","MmeanBSA",4,0,4);
	TH1F *MmeanBSAN = new TH1F("MmeanBSAN","MmeanBSAN",4,0,4);
	
	TH2D *PolarizationTransfer= new TH2D("PolarizationTransfer", "PolarizationTransfer",100,0,11,100,0.,1.);

	//Fiducial cuts
	TH2D *ThetaPhiCut[24];
	for(int sect=0;sect<6;sect++){

		for(int i=0;i<4;i++){

			char* title = new char[10];
			sprintf(title, "sect_%d_pbin%d",sect,i); 
			ThetaPhiCut[(sect*4)+i] = new TH2D(title, title ,150,-100,50,50,0,50);
		} 

	}



	TF1 *pBetaProton= new TF1("pbprot","x/sqrt(0.938*0.938+x*x)",0,2.5);

	TH2D *ChePhi = new TH2D("ChePhi", "ChePhi",200,-180,180,200,0,40);

	TF1 *FitFuncM = new TF1("FitFunc","0.263*(0.985+(-0.036/x)+(0.002/(x*x)))-(5*0.0166)",0,3);

	
	TH1F *MLPplot= new TH1F("MLPplot","",100,-0.05,1.05);
	
	
	//Calculate mean of other variables
	
	TH1F *AFBt0_Eg= new TH1F("AFBt0_Eg","AFBt0_Eg",10,0,11);
	AFBt0_Eg->Sumw2();
	TH1F *AFBt0_Q2= new TH1F("AFBt0_Q2","AFBt0_Q2",10,0,11);
	AFBt0_Q2->Sumw2();
											
	TH1F *AFBEg0_t= new TH1F("AFBEg0_t","AFBEg0_t",10,0,11);
	AFBEg0_t->Sumw2();
	TH1F *AFBEg0_Q2= new TH1F("AFBEg0_Q2","AFBEg0_Q2",10,0,11);
	AFBEg0_Q2->Sumw2();
							
	TH1F *AFBQ20_Eg= new TH1F("AFBQ20_Eg","AFBQ20_Eg",10,0,11);
	AFBQ20_Eg->Sumw2();
	TH1F *AFBQ20_t= new TH1F("AFBQ20_t","AFBQ20_t",10,0,11);
	AFBQ20_t->Sumw2();
							
	TH1F *AFBXi0_Q2= new TH1F("AFBXi0_Q2","AFBXi0_Q2",10,0,11);
	AFBXi0_Q2->Sumw2();
	TH1F *AFBXi0_Eg= new TH1F("AFBXi0_Eg","AFBXi0_Eg",10,0,11);
	AFBXi0_Eg->Sumw2();
	TH1F *AFBXi0_t= new TH1F("AFBXi0_t","AFBXi0_t",10,0,11);
	AFBXi0_t->Sumw2();
								
	TH1F *AFBt1_Eg= new TH1F("AFBt1_Eg","AFBt1_Eg",10,0,11);
	AFBt1_Eg->Sumw2();
	TH1F *AFBt1_Q2= new TH1F("AFBt1_Q2","AFBt1_Q2",10,0,11);
	AFBt1_Q2->Sumw2();
	
	TH1F *AFBt2_Eg= new TH1F("AFBt2_Eg","AFBt2_Eg",10,0,11);
	AFBt2_Eg->Sumw2();
	TH1F *AFBt2_Q2= new TH1F("AFBt2_Q2","AFBt2_Q2",10,0,11);
	AFBt2_Q2->Sumw2();
							
	TH1F *RratioHighMass_Eg= new TH1F("RratioHighMass_Eg","RratioHighMass_Eg",10,0,11);
	RratioHighMass_Eg->Sumw2();
	TH1F *RratioHighMass_Q2= new TH1F("RratioHighMass_Q2","RratioHighMass_Q2",10,0,11);
	RratioHighMass_Q2->Sumw2();
						
	TH1F *Rratio_Eg= new TH1F("Rratio_Eg","Rratio_Eg",10,0,11);
	Rratio_Eg->Sumw2();
	TH1F *Rratio_Q2= new TH1F("Rratio_Q2","Rratio_Q2",10,0,11);
	Rratio_Q2->Sumw2();
	TH1F *Rratio_t= new TH1F("Rratio_t","Rratio_t",10,0,11);
	Rratio_t->Sumw2();
	
	TH1F *BSA_Eg= new TH1F("BSA_Eg","BSA_Eg",10,0,11);
	BSA_Eg->Sumw2();
	TH1F *BSA_Q2= new TH1F("BSA_Q2","BSA_Q2",10,0,11);
	BSA_Q2->Sumw2();
	TH1F *BSA_t= new TH1F("BSA_t","BSA_t",10,0,11);
	BSA_t->Sumw2();
	
	TH1F* BSA_theta = new TH1F("BSA_theta","BSA_theta",10,0,180);
	
	TH1F *RunCheck= new TH1F("RunCheck","RunCheck",10000,5000,5500);
	
	
	//1D exclusivity cuts
	TH1F *Pt1D =new TH1F("Pt1D","Pt/P",100,-0.05,0.2);
	TH1F *MassScat =new TH1F("MassScat","M^2 (GeV)",100,-1.5,1.5);
	
	TH1F *Pt1D1 =new TH1F("Pt1D1","Pt/P",100,-0.05,0.2);
	TH1F *MassScat1 =new TH1F("MassScat1","M^2 (GeV)",100,-1.5,1.5);
	
	
	// check moyenne kinemaics
	TH2F *EgVSQ2 =new TH2F("EgVSQ2","EgVSQ2",20,4,10.6,20,1.5,3);
	
	///////
	//PID and
	//Detector ID
	//////

	double me = 0.000511;
	double mp = 0.938;
	double ebeam = 10.604;

	double Pi =3.14159269;

	bool PIDflag=true;
	bool SFdist=true;

	int LTCC=16;
	int HTCC=15;
	int ECAL=7;

	//ECAL layer
	int PCAL=1;
	int ECIN=4;
	int ECOUT=7;
	
	//chi2cut proton
    	double meanFD = 0.26;
   	double meanCD = 0.81;
   	double sigmaFD = 1.207;
   	double sigmaCD = 1.972;
	
	TF1 *f1P = new TF1("f1P","[0]+[1]*x+[2]*x*x",0,1);
	f1P->SetParameters(0.153319,-0.298968,0.1607);
	TF1 *f2P = new TF1("f1P","[0]+[1]*x+[2]*x*x",0,1);
	f2P->SetParameters(0.0398946,-0.0748125,0.0395764);
	TF1 *f3P = new TF1("f1P","[0]+[1]*x",0,1);
	f3P->SetParameters(0.0292947,-0.0577956);

	/*TF1 *f1Pproton = new TF1("f1Pproton","[0]+[1]*x",-180.,180.);
	f1Pproton->SetParameters(-0.0067175,-1.27761e-05);
	TF1 *f2Pproton = new TF1("f2Pproton","[0]+[1]*x",-180.,180.);
	f2Pproton->SetParameters(-0.0712359,0.000864772);
	TF1 *f3Pproton = new TF1("f3Pproton","[0]+[1]*x",-180.,180.);
	f3Pproton->SetParameters(-0.0809343,0.00049828);
*/

	//Inbending
	TF1 *f1Pproton = new TF1("f1Pproton","[0]+[1]*x",-180.,180.);
	f1Pproton->SetParameters(-0.00146068,-2.13735e-05);
	TF1 *f2Pproton = new TF1("f2Pproton","[0]+[1]*x",-180.,180.);
	f2Pproton->SetParameters(-0.0608671,0.000849025);
	TF1 *f3Pproton = new TF1("f3Pproton","[0]+[1]*x",-180.,180.);
	f3Pproton->SetParameters(-0.0670748,0.000419003);
	/*
	//outbending
	TF1 *f1Pproton = new TF1("f1Pproton","[0]+[1]*x",-180.,180.);
	f1Pproton->SetParameters(-0.0111,0.00025);
	TF1 *f2Pproton = new TF1("f2Pproton","[0]+[1]*x",-180.,180.);
	f2Pproton->SetParameters(-0.11497,0.00065);
	TF1 *f3Pproton = new TF1("f3Pproton","[0]+[1]*x",-180.,180.);
	f3Pproton->SetParameters(-0.1802,0.00144);
*/
/*
TF1 *f1Pproton = new TF1("f1Pproton","[0]+[1]*x",-130.,10.);
	f1Pproton->SetParameters(-0.0314093,-5.85409e-05);
	TF1 *f2Pproton = new TF1("f2Pproton","[0]+[1]*x",-100.,40.);
	f2Pproton->SetParameters(-0.0924035,0.000822015);
	TF1 *f3Pproton = new TF1("f3Pproton","[0]+[1]*x",25.,160.);
	f3Pproton->SetParameters(-0.0973604,0.00046165);
	*/	
	cout<<" here "<<endl;


	hipo::bank EVENT(factory.getSchema("REC::Event"));
	hipo::bank PART(factory.getSchema("REC::Particle"));
	hipo::bank SCIN(factory.getSchema("REC::Scintillator"));
	hipo::bank CHE(factory.getSchema("REC::Cherenkov"));
	hipo::bank CALO(factory.getSchema("REC::Calorimeter"));

	hipo::bank RUN(factory.getSchema("RUN::config"));

	hipo::bank MCPART(factory.getSchema("MC::Particle"));
	hipo::bank MCEVENT(factory.getSchema("MC::Event"));
	
	hipo::bank TRACK(factory.getSchema("REC::Track"));
	hipo::bank TRAJ(factory.getSchema("REC::Traj"));


	int seed = atoi(argv[5]);
	TRandom *polarizationGene = new TRandom(seed);
	TRandom3 *ChoiceEvent = new TRandom3(seed);
	
	int corrrad=0;
	int nentriesChe=0;
	int test =0;
	int nbEvent=0;
	int IDrunC=0;
	
	bool Rafo=true;
	while (reader.next() ){

		nbEvent++;
		//if(nbEvent<3500000)continue;

		reader.read(event);

		//writer.writeEvent(event);



		event.getStructure(MCPART);
		event.getStructure(MCEVENT);

		event.getStructure(RUN);

		event.getStructure(PART);
		event.getStructure(SCIN);
		event.getStructure(CHE);
		event.getStructure(CALO);
		event.getStructure(EVENT);

		event.getStructure(TRAJ);
		event.getStructure(TRACK);

		int IDevent=RUN.getInt("event",0);
		int IDrun=RUN.getInt("run",0);
		
		
		
		/*if(IDrun!=IDrunC){cout<<IDrun<<endl;
		IDrunC=IDrun;}*/
		
		
		//if(IDrun<5800)continue;
		//if(IDrun==5581)cout<<"5581"<<endl;
		
		
//cout<<"here"<<endl;
		if(PART.getSize()<1)continue; 

		//PART.show();

		//RunCheck->Fill(IDrun);
		
		//Choose events
		//double choice = ChoiceEvent->Uniform(100.);
		double factionEvent = 0.015*5000000.;
		//if(nbEvent>factionEvent*(seed+1) || nbEvent<factionEvent*(seed))continue;
		
		nentriesChe+=CHE.getRows();


		//il faudra changer la place de ces variables
		TLorentzVector vBeam;
		TLorentzVector vMissing;
		double MMassBeam;

		Particle vProton;
		Particle vElectron;
		Particle vPositron;

		TLorentzVector vRestProton;
		TLorentzVector vPhoton;
		
		TLorentzVector vPhotonAlongZ;
		
		double phiCM;
		double thetaCM;

		double Q2;

		double qp2;
		double t;
		double Epho;

		int pid;
		double p;

		int recem;
		int recep;
		int recp;
		int tcs;

		int recplus;
		int recmoins;

		recem=0;
		recep=0;
		recp=0;

		recplus=0;
		recmoins=0;
		vRestProton.SetPxPyPzE(0.,0.,0.,mp);
		vBeam.SetPxPyPzE(0.,0.,ebeam,ebeam);

		
		int polarization=EVENT.getInt("helicity",0);
		polarization = -1.*polarization;
		
		//int polarization=polarizationGene->Integer(2);
		//if(polarization==0)polarization=-1;
		//cout<<polarization<<endl;

		double MCfluxBH = MCEVENT.getFloat("ebeam",0);
		double MCpsfBH = MCEVENT.getFloat("ptarget",0);
		double MCcsBH = MCEVENT.getFloat("pbeam",0);
		
		double MCcsTOT = MCEVENT.getFloat("weight",0);

		//double w=MCCrossSectionBH[0];//*(3.*17./250.);
		double w=1;//MCfluxBH[0]*MCpsfBH[0]*MCcsBH[0];
		//double w=MCfluxBH*MCpsfBH*MCcsBH;//-2.*(MCcsTOT-(MCfluxBH*MCpsfBH*MCcsBH));//MCcsTOT;//(MCfluxBH*MCpsfBH*MCcsBH);//////

/*
				Particle MCvProton;
				Particle MCvElectron1;
				Particle MCvElectron2;
				Particle MCvPositron;

		int npMC=MCPART.getRows();
//cout<<npMC<<" "<<PART.getRows()<<endl;

				//if(npMC!=3)continue;
				for( int i=0; i < npMC; i++ ){
					int MCpid=MCPART.getInt("pid",i);
					float  px = MCPART.getFloat("px",i);
					float  py = MCPART.getFloat("py",i);
					float  pz = MCPART.getFloat("pz",i);
					float  vx = MCPART.getFloat("vx",i);
					float  vy = MCPART.getFloat("vy",i);
					float  vz = MCPART.getFloat("vz",i);
					if(MCpid==-11){
						MCvPositron.Vector.SetXYZM(px,py,pz,me);
						MCvPositron.index=i;
						MCvPositron.pid=-11;
						MCvPositron.vertex.x=vx;
						MCvPositron.vertex.y=vy;
						MCvPositron.vertex.z=vz;
						

					}
					if(MCpid==11  ){//pair electron
						MCvElectron2.Vector.SetXYZM(px,py,pz,me);
						MCvElectron2.index=i;
						MCvElectron2.pid=11;
						MCvElectron2.vertex.x=vx;
						MCvElectron2.vertex.y=vy;
						MCvElectron2.vertex.z=vz;
						

					}
					
					
					
					if(MCpid==2212){
						MCvProton.Vector.SetXYZM(px,py,pz,mp);
						MCvProton.index=i;
						MCvProton.pid=2212;
						MCvProton.vertex.x=vx;
						MCvProton.vertex.y=vy;
						MCvProton.vertex.z=vz;
						
					}


				}

*/
		/*
		   if(w>200)continue;
		   if(w<0.)continue;
		 */
		nbrecEvent+=w;
		int nbELEC=0;

		//double w=psf[0]*crs_BH[0]*flux_factor[0];

		int np=PART.getRows();
		Particle Photons[np];

		for( int i=0; i < np; i++ ){


			pid=PART.getInt("pid",i);
			float  px = PART.getFloat("px",i);
			float  py = PART.getFloat("py",i);
			float  pz = PART.getFloat("pz",i);
			float  beta = PART.getFloat("beta",i);
			int status = abs(PART.getInt("status",i));
			int charge = PART.getInt("charge",i);
			float  chi2 = PART.getFloat("chi2pid",i);
			float  vx = PART.getFloat("vx",i);
			float  vy = PART.getFloat("vy",i);
			float  vz = PART.getFloat("vz",i);
			float  vt = PART.getFloat("vt",i);

			if(charge>0)recplus++;
			if(charge<0)recmoins++;

			if(pid==-11 ){
				vPositron.Vector.SetXYZM(px,py,pz,me);
				vPositron.index=i;
				vPositron.pid=-11;
				vPositron.beta=beta;
				vPositron.status=status;
				vPositron.chi2=chi2;
				vPositron.vertex.x=vx;
				vPositron.vertex.y=vy;
				vPositron.vertex.z=vz;
				vPositron.vt=vt;
				recep++;


			}

			if(pid==11 ){
				if(status>2000){
					//cout<<status<<endl;
					vElectron.Vector.SetXYZM(px,py,pz,me);
					vElectron.index=i;
					vElectron.pid=11;
					vElectron.beta=beta;
					vElectron.status=status;
					vElectron.chi2=chi2;
					vElectron.vertex.x=vx;
					vElectron.vertex.y=vy;
					vElectron.vertex.z=vz;
					vElectron.vt=vt;
					recem++;


				}
			}

			if(pid==2212){
				vProton.Vector.SetXYZM(px,py,pz,mp);
				vProton.index=i;
				vProton.pid=2212;
				vProton.beta=beta;
				vProton.status=status;
				vProton.chi2=chi2;
				vProton.vertex.x=vx;
				vProton.vertex.y=vy;
				vProton.vertex.z=vz;
				vProton.vt=vt;
				recp++;
			}
			
			if(pid==22 ){
				Photons[i].Vector.SetXYZM(px,py,pz,0.);
				Photons[i].index=i;
				Photons[i].pid=22;
				Photons[i].status=status;
				Photons[i].beta=beta;
			}


		}

		if(!( recem==1  && recep==1  && recp==1 ))continue;
		
		vElectron=ApplyECcuts(vElectron, CALO);
		vPositron=ApplyECcuts(vPositron, CALO);
		vProton=ApplyECcuts(vProton, CALO);

		//cout<<PassElectron<<" "<<PassPositron<<" "<<PassProton<<endl;
		//if(recem==1)cout<<vElectron.Energy(ECAL,PCAL)<<endl;
		vector<Particle> Particles={vPositron,vElectron,vProton};				

				

		if(PIDflag){
			//cout<<" "<<endl;
			CalorimeterResp Calo;
			CheResp Che;
			ScinResp Scin;
			for (int i=0;i<3 ; i++){

				/*for(int c=0;c<CALO.getRows();c++){

					int Calopindex = CALO.getInt("pindex",c);
					int Calosector = CALO.getInt("sector",c);
					int Calolayer = CALO.getInt("layer",c);
					int Calodetector = CALO.getInt("detector",c);
					float Caloenergy = CALO.getFloat("energy",c);
					float Calox = CALO.getFloat("x",c);
					float Caloy = CALO.getFloat("y",c);
					float Caloz = CALO.getFloat("z",c);
					float Calou = CALO.getFloat("lu",c);
					float Calov = CALO.getFloat("lv",c);
					float Calow = CALO.getFloat("lw",c);
					float Calodu = CALO.getFloat("du",c);
					float Calodv = CALO.getFloat("dv",c);
					float Calodw = CALO.getFloat("dw",c);
					float Calom2u = CALO.getFloat("m2u",c);
					float Calom2v = CALO.getFloat("m2v",c);
					float Calom2w = CALO.getFloat("m2w",c);
					float Calom3u = CALO.getFloat("m3u",c);
					float Calom3v = CALO.getFloat("m3v",c);
					float Calom3w = CALO.getFloat("m3w",c);

					if(Calopindex==(Particles[i].index)){
						Calo.detector=Calodetector;
							//cout<<"index "<<Particles[i].index<<" ener "<<Caloenergy<<endl;
						Calo.pindex=Calopindex;
						Calo.sector=Calosector;
						Calo.layer=Calolayer;
						Calo.energy=Caloenergy;
						Calo.x=Calox;
						Calo.y=Caloy;
						Calo.z=Caloz;
						Calo.u=Calou;
						Calo.v=Calov;
						Calo.w=Calow;
						Calo.du=Calodu;
						Calo.dv=Calodv;
						Calo.dw=Calodw;
						Calo.m2u=Calom2u;
						Calo.m2v=Calom2v;
						Calo.m2w=Calom2w;
						Calo.m3u=Calom3u;
						Calo.m3v=Calom3v;
						Calo.m3w=Calom3w;
						Particles[i].Calorimeter.push_back(Calo);
						//cout<<"pid "<<Particles[i].pid<<" ener "<<Particles[i].Energy(ECAL,PCAL)<<endl;	
						//			if(Calo.layer==1 && Calo.energy<0.06)cout<<Particles[i].pid<<" "<<Calo.layer<<" "<<Calo.energy<<endl;
					}
				}*/
				for(int c=0;c<CHE.getRows();c++){
					int Chepindex = CHE.getInt("pindex",c);
					int Chedetector = CHE.getInt("detector",c);
					int Chesector = CHE.getInt("sector",c);
					float Chenphe = CHE.getFloat("nphe",c);
					float Chetime = CHE.getFloat("time",c);
					float Chechi2 = CHE.getFloat("chi2",c);
					float Chex = CHE.getFloat("x",c);
					float Chey = CHE.getFloat("y",c);
					float Chez = CHE.getFloat("z",c);
					//cout<<Chechi2<<endl;
					if(Chepindex==(Particles[i].index)){
						Che.detector=Chedetector;
						Che.pindex=Chepindex;
						Che.sector=Chesector;
						Che.nphe=Chenphe;
						Che.time=Chetime;
						Che.chi2=Chechi2;
						Che.x=Chex;
						Che.y=Chey;
						Che.z=Chez;
						Particles[i].Cherenkov.push_back(Che);

					}
				}
				for(int c=0;c<SCIN.getRows();c++){
					int Scindetector = SCIN.getInt("detector",c);
					int Scinpindex = SCIN.getInt("pindex",c);
					float Scintime = SCIN.getFloat("time",c);
					float Scinpath = SCIN.getFloat("path",c);
					float Scinenergy = SCIN.getFloat("energy",c);
					int Scinsector = SCIN.getInt("sector",c);

					if(Scinpindex==(Particles[i].index)){
						Scin.detector=Scindetector;
						Scin.pindex=Scinpindex;
						Scin.t=Scintime;
						Scin.path=Scinpath;
						Scin.energy=Scinenergy;
						Scin.sector=Scinsector;
						//cout<<Scinsector[c]<<endl;
						if(Particles[i].Scintillator.energy<Scinenergy){Particles[i].Scintillator=Scin;};

					}
				}

			}

		}

		vElectron=Particles[1];
		vPositron=Particles[0];
		vProton=Particles[2];


		
		
		//Radcoor
		bool plotradcorr=false;
		if(plotradcorr){	
		int inloop=0;
		for(int i=0;i<np;i++){
			if(/*Photons[i].Vector.Angle(vElectron.Vector.Vect())*TMath::RadToDeg()<10. && */Photons[i].Vector.P()>0. /*&& inloop==0*/){
			inloop=1;
			//cout<<(Photons[i].Vector.Theta()*TMath::RadToDeg())<<" "<<(vElectron.Vector.Theta()*TMath::RadToDeg())<<endl;
				if(/*abs(Photons[i].Vector.Theta()-vElectron.Vector.Theta())*TMath::RadToDeg()<1.5*//*vElectron.Vector.P()<2.*/true){
					diffPhidiffThetaRad->Fill(Photons[i].Vector.Angle(vElectron.Vector.Vect())*TMath::RadToDeg());
					diffPhidiffThetaRad1->Fill((Photons[i].Vector.Theta()-vElectron.Vector.Theta())*TMath::RadToDeg());
					diffPhidiffThetaRad1vsTheta->Fill((Photons[i].Vector.Theta()-vElectron.Vector.Theta())*TMath::RadToDeg(),vElectron.Vector.Theta()*TMath::RadToDeg());
					diffPhidiffThetaRad1vsP->Fill((Photons[i].Vector.Theta()-vElectron.Vector.Theta())*TMath::RadToDeg(),vElectron.Vector.P());
					diffPhidiffThetaRad1vsPhi->Fill((Photons[i].Vector.Theta()-vElectron.Vector.Theta())*TMath::RadToDeg(),(Photons[i].Vector.Phi()-vElectron.Vector.Phi())*TMath::RadToDeg());
					/*
					diffPhidiffThetaRad->Fill(Photons[i].Vector.Angle(vPositron.Vector.Vect())*TMath::RadToDeg());
					diffPhidiffThetaRad1->Fill((Photons[i].Vector.Theta()-vPositron.Vector.Theta())*TMath::RadToDeg());
					diffPhidiffThetaRad1vsTheta->Fill((Photons[i].Vector.Theta()-vPositron.Vector.Theta())*TMath::RadToDeg(),vPositron.Vector.Theta()*TMath::RadToDeg());
					diffPhidiffThetaRad1vsP->Fill((Photons[i].Vector.Theta()-vPositron.Vector.Theta())*TMath::RadToDeg(),vPositron.Vector.P());*/
					
				}
				
			}
		}
		}
		
		
		//TMVA pid
					
		double SFPCALa =(vPositron.Energy(ECAL,PCAL))/vPositron.Vector.P();       
		double SFECALin =(vPositron.Energy(ECAL,ECIN))/vPositron.Vector.P();     
		double SFECALout =(vPositron.Energy(ECAL,ECOUT))/vPositron.Vector.P(); 

		SFPCAL=SFPCALa;
                SFECIN=SFECALin;
                SFECOUT=SFECALout;
                        
		double M2PCAL=-1.;
		double M2ECIN=-1.;
		double M2ECOUT=-1.;
			
                for(int i=0;i<vPositron.Calorimeter.size();i++){
                        if(PCAL==vPositron.Calorimeter[i].layer){
                         M2PCAL=(vPositron.Calorimeter[i].m2u+vPositron.Calorimeter[i].m2v+vPositron.Calorimeter[i].m2w)/3.;
                  	}
                        if(ECIN==vPositron.Calorimeter[i].layer){
                               M2ECIN=(vPositron.Calorimeter[i].m2u+vPositron.Calorimeter[i].m2v+vPositron.Calorimeter[i].m2w)/3.;
                         }
                        if(ECOUT==vPositron.Calorimeter[i].layer){
                               M2ECOUT=(vPositron.Calorimeter[i].m2u+vPositron.Calorimeter[i].m2v+vPositron.Calorimeter[i].m2w)/3.;
                        }
                 }
                        
                  m2PCAL=M2PCAL;
                  m2ECIN=M2ECIN;
                  m2ECOUT=M2ECOUT;
					
		//score on which to classify the positrons
		double MLPscore = readerTMVA->EvaluateMVA ( "MLP method" );
		//cout<<"first mlp score "<<MLPscore<<endl;
		if(vPositron.Vector.P()>4. && MLPscore<0.5)continue;
		
		/*
		
		///////////////////// Electron MVA
					
		SFPCALa =(vElectron.Energy(ECAL,PCAL))/vElectron.Vector.P();       
		SFECALin =(vElectron.Energy(ECAL,ECIN))/vElectron.Vector.P();     
		SFECALout =(vElectron.Energy(ECAL,ECOUT))/vElectron.Vector.P(); 

		SFPCAL=SFPCALa;
                SFECIN=SFECALin;
                SFECOUT=SFECALout;
                        
		M2PCAL=-1.;
		M2ECIN=-1.;
		M2ECOUT=-1.;
			
                for(int i=0;i<vElectron.Calorimeter.size();i++){
                        if(PCAL==vElectron.Calorimeter[i].layer){
                         M2PCAL=(vElectron.Calorimeter[i].m2u+vElectron.Calorimeter[i].m2v+vElectron.Calorimeter[i].m2w)/3.;
                  	}
                        if(ECIN==vElectron.Calorimeter[i].layer){
                               M2ECIN=(vElectron.Calorimeter[i].m2u+vElectron.Calorimeter[i].m2v+vElectron.Calorimeter[i].m2w)/3.;
                         }
                        if(ECOUT==vElectron.Calorimeter[i].layer){
                               M2ECOUT=(vElectron.Calorimeter[i].m2u+vElectron.Calorimeter[i].m2v+vElectron.Calorimeter[i].m2w)/3.;
                        }
                 }
                        
                  m2PCAL=M2PCAL;
                  m2ECIN=M2ECIN;
                  m2ECOUT=M2ECOUT;
					
		
		MLPscore = readerTMVA->EvaluateMVA ( "MLP method" );
		*/
		//if(/*vPositron.Vector.P()>4. &&*/ MLPscore<0.5)continue;
		////////////////
		
		/*
		//TMVA electron
		
		SFPCALa =(vElectron.Energy(ECAL,PCAL))/vElectron.Vector.P();       
		SFECALin =(vElectron.Energy(ECAL,ECIN))/vElectron.Vector.P();     
		SFECALout =(vElectron.Energy(ECAL,ECOUT))/vElectron.Vector.P(); 

		SFPCAL=SFPCALa;
                SFECIN=SFECALin;
                SFECOUT=SFECALout;
                        
		M2PCAL=-1.;
		M2ECIN=-1.;
		M2ECOUT=-1.;
			
                for(int i=0;i<vElectron.Calorimeter.size();i++){
                        if(PCAL==vElectron.Calorimeter[i].layer){
                         M2PCAL=(vElectron.Calorimeter[i].m2u+vElectron.Calorimeter[i].m2v+vElectron.Calorimeter[i].m2w)/3.;
                  	}
                        if(ECIN==vElectron.Calorimeter[i].layer){
                               M2ECIN=(vElectron.Calorimeter[i].m2u+vElectron.Calorimeter[i].m2v+vElectron.Calorimeter[i].m2w)/3.;
                         }
                        if(ECOUT==vElectron.Calorimeter[i].layer){
                               M2ECOUT=(vElectron.Calorimeter[i].m2u+vElectron.Calorimeter[i].m2v+vElectron.Calorimeter[i].m2w)/3.;
                        }
                 }
                        
                  m2PCAL=M2PCAL;
                  m2ECIN=M2ECIN;
                  m2ECOUT=M2ECOUT;
					
		//score on which to classify the positrons
		
		MLPscore = readerTMVA->EvaluateMVA ( "MLP method" );
		//cout<<" "<<MLPscore<<endl;
		if(vElectron.Vector.P()>4. && MLPscore<0.5)continue;
		*/
		/////////////////////////////////////
		
		
		
		double STT = EVENT.getFloat("startTime",0);
		double protonVT=VT(vProton.Vector.P(),mp,vProton.Scintillator.t,vProton.Scintillator.path,STT);
		
		Particle oldElectron = vElectron;
		double avant=vElectron.Vector.P();
		//cout<<"avant "<<vElectron.Vector.P()<<" status "<<vElectron.status<<endl;
		bool RadCorr = true;//true;//
		if(RadCorr){
		//cout<<"here "<<endl;
		vElectron = RadiativeCorr(vElectron,Photons,10.,1.5,np);
		vPositron = RadiativeCorr(vPositron,Photons,10.,1.5,np);
		
		}
		
		if(avant!=vElectron.Vector.P())DiffRadCor->Fill((vElectron.Vector.P()-avant)/avant);
		
		if(avant!=vElectron.Vector.P()){corrrad++;
		//cout<<"avant "<<avant<<"après "<<vElectron.Vector.P()<<" "<<(avant-vElectron.Vector.P())/avant<<" status "<<vElectron.status<<endl;
		}

		//MC corr first
		double PP=vProton.Vector.P();
		double newPP = 0.0;
		if(vProton.status<4000 && vProton.Vector.Theta()*TMath::RadToDeg()>27.)newPP=PP+(f1P->Eval(PP));
		if(vProton.status<4000 && vProton.Vector.Theta()*TMath::RadToDeg()<27.)newPP=PP+(f2P->Eval(PP));
		if(vProton.status>4000)newPP=PP+(f3P->Eval(PP));
		
		vProton.Vector.SetRho(newPP);
		//cout<<vProton.Vector.P()<<endl;

		//cout<<"new event"<<endl;
		//cout<<vProton.Vector.P()<<endl;
		
		bool correctMom=true;//true;//true;
		if(recem==1  && recep==1 && recp==1  && vElectron.passEC && vPositron.passEC && vProton.passEC && correctMom){
		Track Region1Track;
		for(int i=0;i<TRACK.getRows();i++){
		
			


			int sector = 0;
			float chi2 = -1000.0;
			int NDF = 1;
			int indexTrack=-5000;
		
			if(TRACK.getInt("status",i)>0 && TRACK.getInt("detector",i)==5 ){
				sector=TRACK.getInt("sector",i);
				chi2=TRACK.getFloat("chi2",i);
				NDF=TRACK.getInt("NDF",i);
				indexTrack=TRACK.getInt("index",i);
			}




			float Nchi2=chi2/NDF;


			for(int c=0;c<TRAJ.getRows();c++){

				int detector=0;
				int pindex = TRAJ.getInt("pindex",c);
				int index = TRAJ.getInt("index",c);
				int layer = TRAJ.getInt("layer",c);
				detector = TRAJ.getInt("detector",c);
				float x = TRAJ.getFloat("x",c);
				float y = TRAJ.getFloat("y",c);
				float z = TRAJ.getFloat("z",c);

				//Region1
				if(index==(indexTrack) && detector==5 && layer==12 && pindex==vProton.index ){
				//cout<<layer<<endl;
					Region1Track.detector=detector;
					Region1Track.pindex=index;
					Region1Track.sector=sector;
					Region1Track.layer=layer;
					Region1Track.x=x;
					Region1Track.y=y;
					Region1Track.z=z;
					Region1Track.chi2=Nchi2;

				}


			}
			}
		
		
		

		double x=Region1Track.x;
		double y=Region1Track.y;
		double z=Region1Track.z;
		double Nchi2=Region1Track.chi2;
		float d = sqrt(x*x+y*y);
		float phi = (TMath::ATan2(y,x)*TMath::RadToDeg());
		float theta = TMath::ATan(d/z)*TMath::RadToDeg();

		
							
		double PP=vProton.Vector.P();
		//if(correctMom && vProton.status>4000)cout<<" PP "<<PP<< " status "<<vProton.status<<" phi "<<phi<<endl;
		double newPP = 0.0;
		if(phi>150. || phi<-90.){
			double phi1=phi;
					if(phi>0.){phi1=phi-270.;}
					else {phi1=phi+90.;}
					newPP=PP*(1.-(f1Pproton->Eval(phi1)));
					}
		if(phi>-90. && phi<30.){newPP=PP*(1.-(f2Pproton->Eval(phi)));}
		if(phi<150. && phi>30.){newPP=PP*(1.-(f3Pproton->Eval(phi)));}
		if(correctMom && vProton.status>4000)vProton.Vector.SetRho(newPP);
		//cout<<vProton.Vector.P()<<endl;
		
		}
		
		
		
		
		
		//if(recem==1)cout<<vElectron.Energy(ECAL,PCAL)<<endl;
		//if( recem==1  && recep==1 && recp==1)PPosiPElec->Fill(vPositron.Vector.P(),vElectron.Vector.P());
		if( /*avant!=vElectron.Vector.P() &&*/ recem==1  && recep==1 && recp==1  && vElectron.passEC && vPositron.passEC && vProton.passEC /*&& recplus==2 && recmoins==1*/ /*&& vElectron.nphe(HTCC)>1
				&& vPositron.nphe(HTCC)>1 */){
				
				//vElectron= oldElectron;

			/*cout<<" NEW EVENT"<<endl;
			  cout<<" PARTICLE "<<endl;
			  PART.show();
			  cout<<endl;
			  cout<<endl;
			  cout<<endl;
			  cout<<endl;*/
			/*cout<<" CHERENKOV "<<endl;
			  CHE.show();
			  cout<<endl;
			  cout<<endl;
			  cout<<endl;
			  cout<<endl;
			 */
			 
			 
			 
			 //cout<<vProton.Vector.P()<<endl;
			/* vProton.Vector=MCvProton.Vector;
			 vElectron.Vector=MCvElectron2.Vector;
			 vPositron.Vector=MCvPositron.Vector;*/
			 //cout<<"apres cut "<<vProton.Vector.P()<<endl;
			 
			 
			
			 
			double enerEPCAL=(vElectron.Energy(ECAL,PCAL));
			double enerPPCAL=(vPositron.Energy(ECAL,PCAL));
			//if(enerEPCAL<0.030)continue;
			//if(enerPPCAL<0.030)continue;

			ChePhi->Fill(vPositron.Vector.Phi()*TMath::RadToDeg(),vPositron.nphe(HTCC));
			
			
			qp2 = (vPositron.Vector+vElectron.Vector).M2();

			t = (vRestProton-vProton.Vector).M2();
			
			
			vPhoton = vProton.Vector+vPositron.Vector+vElectron.Vector-vRestProton;
			Epho=vPhoton.E();
			
			vPhotonAlongZ.SetXYZM(0.0,0.0,vPhoton.Pz(),0.0);
			double tFromPhoton = -1.*(vPositron.Vector+vElectron.Vector-vPhotonAlongZ).M2();
			//cout<<"t check "<<tFromPhoton<<" "<<(-t)<<endl;
			tCheckBefore->Fill(tFromPhoton+t,w);
					
					
			vMissing=vBeam-vPhoton;
			MMassBeam = (vMissing).M2();
			double Q2prim=(vPhoton).M2();
			Q2=2*ebeam*vMissing.E()*(1.-cos(vMissing.Theta()));			

			histPMiss->Fill(vMissing.Px()/vMissing.P(),vMissing.Py()/vMissing.P(),w);
			histMMass->Fill(MMassBeam,w);
			
			Pt1D->Fill(vMissing.Pt()/vMissing.P(),w);
			MassScat->Fill(MMassBeam);
			
			if(abs(MMassBeam)<0.4)Pt1D1->Fill(vMissing.Pt()/vMissing.P(),w);
			if(abs(vMissing.Pt()/vMissing.P())<0.05)MassScat1->Fill(MMassBeam,w);
			
			histQ2->Fill(Q2,w);			
			histQP2->Fill(qp2,w);
			histM->Fill((vPositron.Vector+vElectron.Vector).M(),w);

			histME->Fill(Epho,(vPositron.Vector+vElectron.Vector).M(),w); 
			histMPElec->Fill(vElectron.Vector.P(),(vPositron.Vector+vElectron.Vector).M(),w);
			histMPPosi->Fill(vPositron.Vector.P(),(vPositron.Vector+vElectron.Vector).M(),w);
			histMPhiProton->Fill(vProton.Vector.Phi()*TMath::RadToDeg(),(vPositron.Vector+vElectron.Vector).M(),w);					
			histMME->Fill(Epho,MMassBeam,w);
			histMMPt->Fill(MMassBeam, vMissing.Pt()/vMissing.P(),w);
			histT->Fill(-t,w);

			EnergyConservation->Fill(vProton.Vector.Pz()+vPositron.Vector.Pz()+vElectron.Vector.Pz()-vPhoton.E(),w);					

			MMQ2->Fill(MMassBeam,Q2,w);

			PhiVSPElectron-> Fill(vElectron.Vector.Phi()*TMath::RadToDeg(),vElectron.Vector.P(),w);
			PhiVSPPositron-> Fill(vPositron.Vector.Phi()*TMath::RadToDeg(),vPositron.Vector.P(),w);
			if(vProton.status>4000)PhiVSPProton-> Fill(vProton.Vector.Phi()*TMath::RadToDeg(),vProton.Vector.P(),w);
			ThetaVSPElectron-> Fill(vElectron.Vector.P(),vElectron.Vector.Theta()*TMath::RadToDeg(),w);
			ThetaVSPPositron-> Fill(vPositron.Vector.P(),vPositron.Vector.Theta()*TMath::RadToDeg(),w);

			PhiVSThetaElectron-> Fill(vElectron.Vector.Phi()*TMath::RadToDeg(),vElectron.Vector.Theta()*TMath::RadToDeg(),w);
			PhiVSThetaPositron-> Fill(vPositron.Vector.Phi()*TMath::RadToDeg(),vPositron.Vector.Theta()*TMath::RadToDeg(),w);
			PhiVSThetaProton-> Fill(vProton.Vector.Phi()*TMath::RadToDeg(),vProton.Vector.Theta()*TMath::RadToDeg(),w);

			ThetaVSPProton->Fill(vProton.Vector.P(),vProton.Vector.Theta()*TMath::RadToDeg(),w);


			double SFelectronPCAL =(vElectron.Energy(ECAL,PCAL))/vElectron.Vector.P();
			double SFpositronPCAL =(vPositron.Energy(ECAL,PCAL))/vPositron.Vector.P();
			double SFelectronECALin =(vElectron.Energy(ECAL,ECIN))/vElectron.Vector.P();
			double SFpositronECALin =(vPositron.Energy(ECAL,ECIN))/vPositron.Vector.P();
			double SFelectronECALout =(vElectron.Energy(ECAL,ECOUT))/vElectron.Vector.P();
			double SFpositronECALout =(vPositron.Energy(ECAL,ECOUT))/vPositron.Vector.P();

			ECelectron->Fill(vElectron.Energy(ECAL,PCAL),vElectron.Energy(ECAL,ECIN)+vElectron.Energy(ECAL,ECOUT));
			ECpositron->Fill(vPositron.Energy(ECAL,PCAL),vPositron.Energy(ECAL,ECIN)+vPositron.Energy(ECAL,ECOUT));

			double enerE=(vElectron.Energy(ECAL,PCAL)+vElectron.Energy(ECAL,ECIN)+vElectron.Energy(ECAL,ECOUT));
			double enerP=(vPositron.Energy(ECAL,PCAL)+vPositron.Energy(ECAL,ECIN)+vPositron.Energy(ECAL,ECOUT));
			double correnerE=0.263*(0.985+(-0.036/enerE)+(0.002/(enerE*enerE)));
			double correnerP=0.263*(0.985+(-0.036/enerP)+(0.002/(enerP*enerP)));
			SFelectron->Fill(vElectron.Vector.P(), enerE/vElectron.Vector.P());
			SFpositron->Fill(vPositron.Vector.P(), enerP/vPositron.Vector.P());


			if(vPositron.Energy(ECAL,ECIN)>0.0)SFPCALECALPosiavant->Fill(vPositron.Energy(ECAL,PCAL)/vPositron.Vector.P(), vPositron.Energy(ECAL,ECIN)/vPositron.Vector.P());
			if(vElectron.Energy(ECAL,ECIN)>0.0)SFPCALECALElecavant->Fill(vElectron.Energy(ECAL,PCAL)/vElectron.Vector.P(), vElectron.Energy(ECAL,ECIN)/vElectron.Vector.P());

			corrSFelectron->Fill(vElectron.Vector.P(), correnerE/vElectron.Vector.P());
			corrSFpositron->Fill(vPositron.Vector.P(), correnerP/vPositron.Vector.P());

			CheElectron->Fill(vElectron.nphe(HTCC));
			ChePositron->Fill(vPositron.nphe(HTCC));			

			BetaProton->Fill(vProton.Vector.P(),vProton.beta);




			//vertex analysis
			double elecVT=VT(vElectron.Vector.P(),me,vElectron.Scintillator.t,vElectron.Scintillator.path,STT);
			double posiVT=VT(vPositron.Vector.P(),me,vPositron.Scintillator.t,vPositron.Scintillator.path,STT);
			//double protonVT=VT(vProton.Vector.P(),mp,vProton.Scintillator.t,vProton.Scintillator.path,STT);

			vertexTimePair->Fill(elecVT-posiVT,w);
			vertexTimeP->Fill(protonVT,w);
			vertexTimePP->Fill(vProton.Vector.P(),protonVT,w);

			vertexPair->Fill(vElectron.vertex.z-vPositron.vertex.z,w);
			vertexElecP->Fill(vElectron.vertex.z-vProton.vertex.z,w);
			vertexPosiP->Fill(vPositron.vertex.z-vProton.vertex.z,w);


			vertexProtTheta->Fill(vProton.Vector.Theta()*TMath::RadToDeg(),vProton.vertex.z,w);


			histMPhiElectron->Fill(vElectron.Vector.Phi()*TMath::RadToDeg(),(vPositron.Vector+vElectron.Vector).M(),w);
			histMPhiPositron->Fill(vPositron.Vector.Phi()*TMath::RadToDeg(),(vPositron.Vector+vElectron.Vector).M(),w);


			//if(vElectron.nphe(HTCC)<1)cout<< vElectron.nphe(HTCC) <<endl;
			int SectorDiff=abs(vElectron.Scintillator.sector-vPositron.Scintillator.sector);
			if(SectorDiff==4)SectorDiff=2;
			if(SectorDiff==5)SectorDiff=1;

			if((vProton.status/1000)==2)vertexProtThetaFD->Fill(vProton.Vector.Theta()*TMath::RadToDeg(),vProton.vertex.z,w);
			if((vProton.status/1000)==4)vertexProtThetaCD->Fill(vProton.Vector.Theta()*TMath::RadToDeg(),vProton.vertex.z,w);
			vertexElec->Fill(vElectron.vertex.z,w);
			vertexPosi->Fill(vPositron.vertex.z,w);
			vertexProt->Fill(vProton.vertex.z,w);
			
			
			
			Chi2ElectronAvant->Fill(vElectron.chi2,w);
			Chi2PositronAvant->Fill(vPositron.chi2,w);
			if(vProton.status>4000)Chi2ProtonAvant->Fill(vProton.chi2,w);
			if(vProton.status<4000)Chi2ProtonFDAvant->Fill(vProton.chi2,w);
			
			
			
			if(abs(MMassBeam)<0.4)PhiMissingParticle->Fill(vMissing.Phi()*TMath::RadToDeg(),abs(vMissing.Pt()/vMissing.P()),w);
			
			//TMVA pid
					
			/*double SFPCALa =(vPositron.Energy(ECAL,PCAL))/vPositron.Vector.P();       
			double SFECALin =(vPositron.Energy(ECAL,ECIN))/vPositron.Vector.P();     
			double SFECALout =(vPositron.Energy(ECAL,ECOUT))/vPositron.Vector.P(); 

			SFPCAL=SFPCALa;
                        SFECIN=SFECALin;
                       	SFECOUT=SFECALout;
                        
			double M2PCAL=-1.;
			double M2ECIN=-1.;
			double M2ECOUT=-1.;
			
                        	for(int i=0;i<vPositron.Calorimeter.size();i++){
                                if(PCAL==vPositron.Calorimeter[i].layer){
                                        M2PCAL=(vPositron.Calorimeter[i].m2u+vPositron.Calorimeter[i].m2v+vPositron.Calorimeter[i].m2w)/3.;
                                }
                                if(ECIN==vPositron.Calorimeter[i].layer){
                                	M2ECIN=(vPositron.Calorimeter[i].m2u+vPositron.Calorimeter[i].m2v+vPositron.Calorimeter[i].m2w)/3.;
                                }
                                if(ECOUT==vPositron.Calorimeter[i].layer){
                                	M2ECOUT=(vPositron.Calorimeter[i].m2u+vPositron.Calorimeter[i].m2v+vPositron.Calorimeter[i].m2w)/3.;
                                }
                        }
                        
                        m2PCAL=M2PCAL;
                        m2ECIN=M2ECIN;
                        m2ECOUT=M2ECOUT;
					
			//score on which to classify the positrons
			double MLPscore = readerTMVA->EvaluateMVA ( "MLP method" );*/
			
			
			if(
					vElectron.Vector.P()>1. 
					&& vPositron.Vector.P()>4.
					&& abs(vMissing.Pt()/vMissing.P())<0.05
					&& abs(MMassBeam)<0.4
					&& enerE/vElectron.Vector.P()>0.15
					&& enerP/vPositron.Vector.P()>0.15
					

			){
			
			//if(MLPscore<0.6 && MLPscore>0.5 && vPositron.Vector.P()>4.)cout<<"t "<<(-t)<<" QP2 "<<qp2<<" mlp score "<<MLPscore<<endl;
			MLPplot->Fill(MLPscore,w);
			}
			
			
			double cutchi2 = 3.;
                    bool chi2cutProton = (vProton.status>4000 && abs(vProton.chi2-meanCD)<(cutchi2*sigmaCD))||(vProton.status<4000 && abs(vProton.chi2-meanFD)<(cutchi2*sigmaFD));
                    
                      bool chi2cutPositron = abs(vPositron.chi2)<3.;

			
			if(

					//(vElectron.Energy(ECAL,PCAL)/vElectron.Vector.P())>0.05
					//&& (vPositron.Energy(ECAL,PCAL)/vPositron.Vector.P())>0.05 
					vElectron.Vector.P()>1. 
					//&& ((vPositron.Vector.P()>1. && vPositron.Vector.P()<=4. ) || (vPositron.Vector.P()>4. && MLPscore>0.5))
					&& vPositron.Vector.P()>1.
					
					&& abs(vMissing.Pt()/vMissing.P())<0.05//0.04//0.05//0.025
					//&& abs(vMissing.M())<0.1
					&& abs(MMassBeam)<0.4//0.3//0.4//0.2
					//&& MMassBeam<0.7 && MMassBeam>-0.1
					&& enerE/vElectron.Vector.P()>0.15
					&& enerP/vPositron.Vector.P()>0.15
					
					//&& chi2cutPositron
					//&& chi2cutProton
					//&& abs(vElectron.TimeChe(HTCC)-vPositron.TimeChe(HTCC))<5.
					//&& abs(vElectron.Calorimeter[0].sector-vPositron.Calorimeter[0].sector)>1
					//&& vProton.Vector.P()>0.1
					//&& abs(vProton.chi2)<5.
					//&& SectorDiff==0
					//&& abs(vElectron.TimeChe(HTCC)-vPositron.TimeChe(HTCC))<5.
					//&& (vPositron.Vector+vElectron.Vector).M()>0.1
					//&& (vPositron.Vector+vElectron.Vector).M()<0.2
					//&& abs(Q2)<0.5
					//&& vElectron.nphe(HTCC)>5.
					//&& vPositron.nphe(HTCC)>5.

					//fiducial cut acceptance
					
					 

				){
				
				/*if((qa->Golden(IDrun,IDevent))==false){
				//cout<<"not golden"<<endl;
				continue;
				}*/
		
		
//cout<<"electron "<<vElectron.Calorimeter[0].sector<<" "<<vElectron.TimeChe(HTCC)<<endl;
//cout<<"positron "<<vPositron.Calorimeter[0].sector<<endl;
				Positronpt->Fill(vPositron.Vector.Pt());
		//cout<<"coucou "<<Epho<<" "<<-t<<" "<<qp2<<endl;
Q2t->Fill(-t,qp2);
					if(vPositron.Energy(ECAL,ECIN)>0.0)SFPCALECALPosiapres->Fill(vPositron.Energy(ECAL,PCAL)/vPositron.Vector.P(), vPositron.Energy(ECAL,ECIN)/vPositron.Vector.P());
					if(vElectron.Energy(ECAL,ECIN)>0.0)SFPCALECALElecapres->Fill(vElectron.Energy(ECAL,PCAL)/vElectron.Vector.P(), vElectron.Energy(ECAL,ECIN)/vElectron.Vector.P());

					
					//cout<<"event"<<endl;
					//cout<<"nb positif "<< recplus<<" nb negatif "<< recmoins<<endl;


					if((vPositron.Vector+vElectron.Vector).M()>3.)nbJPSI++;
					
					/*if((vPositron.Vector+vElectron.Vector).M()>1.5 && vPositron.Vector.P()>4.)*/
					
					
					//if(enerP/vPositron.Vector.P()<0.15)cout<<"positron "<<enerP/vPositron.Vector.P()<<" "<<vPositron.Vector.P()<<endl;
					//if(enerE/vElectron.Vector.P()<0.15){cout<<"electron "<<vElectron.status<<" "<<enerE<<" "<<vElectron.Energy(ECAL,PCAL)<<" "<<vElectron.Vector.P()<<endl;CALO.show();PART.show();}
					if(true/*vPositron.Energy(ECAL,ECIN)==0.0*/)SFPpcalzeroPosi->Fill(vPositron.Vector.P(),enerP/vPositron.Vector.P());
					if(true/*vElectron.Energy(ECAL,ECIN)==0.0*/)SFPpcalzeroElec->Fill(vElectron.Vector.P(),enerE/vElectron.Vector.P());
					
					
					
					
					
					
					int np=PART.getRows();
					bool noOtherPart = true;
					for( int i=0; i < np; i++ ){
						int pid=PART.getInt("pid",i);
						int status = abs(PART.getInt("status",i));
						float  chi2 = PART.getFloat("chi2pid",i);
						float  vz = PART.getFloat("vz",i);
						if(pid==2112 || pid==-11 || pid ==2212 || pid==22 || pid==11 || pid==0 )continue;
						//cout<<pid<<" status "<<status<<" vz "<<vz<<" chi2 "<<chi2<<endl;
						noOtherPart=false;
						VertexOtherPart->Fill(vz);
						Chi2OtherPart->Fill(chi2);
					}


					float UlimitsA[6] = {45 , 53 , 55 , 39 , 53 , 47};
					float VlimitsA[6] = { 27 , 11 , 19 , 15 , 19 , 15};
					float WlimitsA[6] = { 17 , 23 , 23 , 11 , 19 , 17};






					float UlimitsB[6] = { 128 , 112 , 132 , 84 , 116 , 108};
					float VlimitsB[6] = { 29 , 17 , 21 , 13 , 25 , 19};
					float WlimitsB[6] = { 19 , 25 , 23 , 11 , 23 , 13};

					int passElec=1;
					int passPosi=1;

					for(int i=0; i<vElectron.Calorimeter.size();i++){
						if(vElectron.Calorimeter[i].layer==PCAL){	
							float u=vElectron.Calorimeter[i].u;
							float v=vElectron.Calorimeter[i].v;
							float w=vElectron.Calorimeter[i].w;
							int sector=vElectron.Calorimeter[i].sector;
							if((vPositron.Vector+vElectron.Vector).M()>0.1 && (vPositron.Vector+vElectron.Vector).M()<0.2 && abs(vElectron.TimeChe(HTCC)-vPositron.TimeChe(HTCC))<4.)BeforeCuts->Fill(vElectron.Calorimeter[i].x,vElectron.Calorimeter[i].y);
							//cout<<u<<" "<<v<<" "<<w<<endl;
							if(u<UlimitsA[sector-1] || v<VlimitsA[sector-1] || w<WlimitsA[sector-1]){
								passElec=0;
							}
						}
					}

					for(int i=0; i<vPositron.Calorimeter.size();i++){
						if(vPositron.Calorimeter[i].layer==PCAL){	
							float u=vPositron.Calorimeter[i].u;
							float v=vPositron.Calorimeter[i].v;
							float w=vPositron.Calorimeter[i].w;
							int sector=vPositron.Calorimeter[i].sector;
							if((vPositron.Vector+vElectron.Vector).M()>0.1 && (vPositron.Vector+vElectron.Vector).M()<0.2 && abs(vElectron.TimeChe(HTCC)-vPositron.TimeChe(HTCC))<4.)BeforeCutsPositron->Fill(vPositron.Calorimeter[i].x,vPositron.Calorimeter[i].y);
							//cout<<u<<" "<<v<<" "<<w<<endl;
							if(u<UlimitsB[sector-1] || v<VlimitsB[sector-1] || w<WlimitsB[sector-1]){
								passPosi=0;
							}
						}
					}

					if((vPositron.Vector+vElectron.Vector).M()>0.1 && (vPositron.Vector+vElectron.Vector).M()<0.2){
						CheCooElec->Fill(vElectron.XChe(HTCC),vElectron.YChe(HTCC));
						CheCooPosi->Fill(vPositron.XChe(HTCC),vPositron.YChe(HTCC));
					}
					//cout<<pass<<endl;	
					if(/*passElec==1 &&  passPosi==1 && */abs(vElectron.TimeChe(HTCC)-vPositron.TimeChe(HTCC))<4.){	

						if(qp2>4. && qp2<9.)AfterCuts++;}	

					



					histAngleMass->Fill(vElectron.Vector.Vect().Angle(vPositron.Vector.Vect())*TMath::RadToDeg(),(vPositron.Vector+vElectron.Vector).M(),w);
					histMThetaElectron->Fill(vElectron.Vector.Theta()*TMath::RadToDeg(),(vPositron.Vector+vElectron.Vector).M(),w);
					histMThetaPositron->Fill(vPositron.Vector.Theta()*TMath::RadToDeg(),(vPositron.Vector+vElectron.Vector).M(),w);
					histMThetaProton->Fill(vProton.Vector.Theta()*TMath::RadToDeg(),(vPositron.Vector+vElectron.Vector).M(),w);


					histMVElectron->Fill(vElectron.vertex.x-vPositron.vertex.x,(vPositron.Vector+vElectron.Vector).M(),w);
					histMVPositron->Fill(vElectron.vertex.y-vPositron.vertex.y,(vPositron.Vector+vElectron.Vector).M(),w);
					histMV->Fill(vElectron.vertex.z-vPositron.vertex.z,(vPositron.Vector+vElectron.Vector).M(),w);

					histMCheElectron->Fill(vElectron.nphe(HTCC),(vPositron.Vector+vElectron.Vector).M(),w);
					histMChePositron->Fill(vPositron.nphe(HTCC),(vPositron.Vector+vElectron.Vector).M(),w);


					if((vPositron.Vector+vElectron.Vector).M()>0.1 && (vPositron.Vector+vElectron.Vector).M()<0.2){
						histCheCheElectron->Fill(vElectron.nphe(HTCC),vPositron.nphe(HTCC),w);
						VertexPio->Fill(vElectron.vertex.z,vPositron.vertex.z);
						MomentumPairPio->Fill(vElectron.Vector.P(),vPositron.Vector.P());
						ThetaPairPio->Fill(vElectron.Vector.Theta()*TMath::RadToDeg(),vPositron.Vector.Theta()*TMath::RadToDeg());
						TimeChePio->Fill(vElectron.TimeChe(HTCC),vPositron.TimeChe(HTCC));
						CheCooDiff->Fill(vElectron.XChe(HTCC)-vPositron.XChe(HTCC),vElectron.YChe(HTCC)-vPositron.YChe(HTCC));

						//cout<<IDrun<<" event "<<IDevent<<" position in the file "<<nbEvent<<endl;
						/*cout<<endl;
						  cout<<"Event "<<(vPositron.Vector+vElectron.Vector).M()<<endl;
						  PART.show();
						  CHE.show();
						  cout<<vElectron.TimeChe(HTCC)<<endl;
						  cout<<vElectron.Chi2Che(HTCC)<<endl;*/
					}
					if( (vPositron.Vector+vElectron.Vector).M()>0.0 && (vPositron.Vector+vElectron.Vector).M()<0.3){
						FuckingPi0->Fill((vPositron.Vector+vElectron.Vector).M(),w);
					}

					if( (vPositron.Vector+vElectron.Vector).M()>0.2 ){
						VertexAutre->Fill(vElectron.vertex.z,vPositron.vertex.z);
						MomentumPairAutre->Fill(vElectron.Vector.P(),vPositron.Vector.P());
						ThetaPairAutre->Fill(vElectron.Vector.Theta()*TMath::RadToDeg(),vPositron.Vector.Theta()*TMath::RadToDeg());
						histCheCheElectronSup->Fill(vElectron.nphe(HTCC),vPositron.nphe(HTCC),w);
						TimeChePioSup->Fill(vElectron.TimeChe(HTCC),vPositron.TimeChe(HTCC));	
						CheCooDiffSup->Fill(vElectron.XChe(HTCC)-vPositron.XChe(HTCC),vElectron.YChe(HTCC)-vPositron.YChe(HTCC));
					}

					if((vPositron.Vector+vElectron.Vector).M()<0.1 ){
						VertexAutreInf->Fill(vElectron.vertex.z,vPositron.vertex.z);
						MomentumPairAutreInf->Fill(vElectron.Vector.P(),vPositron.Vector.P());
						ThetaPairAutreInf->Fill(vElectron.Vector.Theta()*TMath::RadToDeg(),vPositron.Vector.Theta()*TMath::RadToDeg());
						histCheCheElectronInf->Fill(vElectron.nphe(HTCC),vPositron.nphe(HTCC),w);	
						CheCooDiffInf->Fill(vElectron.XChe(HTCC)-vPositron.XChe(HTCC),vElectron.YChe(HTCC)-vPositron.YChe(HTCC));
						TimeChePioInf->Fill(vElectron.TimeChe(HTCC),vPositron.TimeChe(HTCC));	

					}

					 /* if(abs(vElectron.Calorimeter[0].sector-vPositron.Calorimeter[0].sector)==0)*/EhistM->Fill((vPositron.Vector+vElectron.Vector).M(),w);

					TimeCheDiff->Fill(vElectron.TimeChe(HTCC)-vPositron.TimeChe(HTCC));
					TimeCheDiffMass->Fill((vPositron.Vector+vElectron.Vector).M(),vElectron.TimeChe(HTCC)-vPositron.TimeChe(HTCC));


					double enerEpcal=vElectron.Energy(ECAL,PCAL);
					double enerPpcal=vPositron.Energy(ECAL,PCAL);
					SFMassElec->Fill((vPositron.Vector+vElectron.Vector).M(),enerEpcal/vElectron.Vector.P());
					if(vPositron.Energy(ECAL,ECIN)>0.0)SFMassPosi->Fill((vPositron.Vector+vElectron.Vector).M(),enerPpcal/vPositron.Vector.P());
					double Mscat=(vPositron.Vector+vMissing).M();
					EhistEssai->Fill(Mscat,w);

					EECelectron->Fill(vElectron.Energy(ECAL,PCAL),vElectron.Energy(ECAL,ECIN)+vElectron.Energy(ECAL,ECOUT),w);
					EECpositron->Fill(vPositron.Energy(ECAL,PCAL),vPositron.Energy(ECAL,ECIN)+vPositron.Energy(ECAL,ECOUT),w);						
					EhistPMiss->Fill(vMissing.Px()/vMissing.P(),vMissing.Py()/vMissing.P(),w);
					EhistMMass->Fill(MMassBeam,w);
					EhistQ2->Fill(Q2,w);
					EhistQP2->Fill(qp2,w);
					
					EhistMAfterFid->Fill((vPositron.Vector+vElectron.Vector).M(),w);


					EhistMSectorDiff->Fill((vPositron.Vector+vElectron.Vector).M(),SectorDiff,w);

					EhistMMPt->Fill(MMassBeam, vMissing.Pt()/vMissing.P(),w);
					EhistT->Fill(-t,w);
					EhistME->Fill(Epho,(vPositron.Vector+vElectron.Vector).M(),w);
					EhistMPElec->Fill(vElectron.Vector.P(),(vPositron.Vector+vElectron.Vector).M(),w);
					EhistMPPosi->Fill(vPositron.Vector.P(),(vPositron.Vector+vElectron.Vector).M(),w);
					EhistMPProton->Fill(vProton.Vector.P(),(vPositron.Vector+vElectron.Vector).M(),w);
					EhistMME->Fill(Epho,MMassBeam,w);
					EhistMXProton->Fill(vProton.chi2,(vPositron.Vector+vElectron.Vector).M(),w);
					//	EhistMRun->Fill(RUNn[0],(vPositron.Vector+vElectron.Vector).M());
					EhistMPhiProton->Fill(vProton.Vector.Phi()*TMath::RadToDeg(),(vPositron.Vector+vElectron.Vector).M(),w);
					EEnergyConservation->Fill(vProton.Vector.Pz()+vPositron.Vector.Pz()+vElectron.Vector.Pz()-vPhoton.E(),w);
					
					
					//vPhotonAlongZ.SetXYZM(0.0,0.0,vPhoton.Pz(),0.0);
					//double tFromPhoton = -1.*(vPositron.Vector+vElectron.Vector-vPhotonAlongZ).M2();
					//cout<<"t check "<<tFromPhoton<<" "<<(-t)<<endl;
					tCheck->Fill(tFromPhoton+t,w);
					
					
					ECheElectron->Fill(vElectron.nphe(HTCC),w);
					EChePositron->Fill(vPositron.nphe(HTCC),w);

					EMMQ2->Fill(MMassBeam,Q2,w);

															
					EThetaVSPElectron-> Fill(vElectron.Vector.P(),vElectron.Vector.Theta()*TMath::RadToDeg(),w);
					EThetaVSPPositron-> Fill(vPositron.Vector.P(),vPositron.Vector.Theta()*TMath::RadToDeg(),w);

					EPhiVSThetaElectron-> Fill(vElectron.Vector.Phi()*TMath::RadToDeg(),vElectron.Vector.Theta()*TMath::RadToDeg(),w);
					EPhiVSThetaPositron-> Fill(vPositron.Vector.Phi()*TMath::RadToDeg(),vPositron.Vector.Theta()*TMath::RadToDeg(),w);
					EPhiVSThetaProton-> Fill(vProton.Vector.Phi()*TMath::RadToDeg(),vProton.Vector.Theta()*TMath::RadToDeg(),w);

					EThetaVSPProton->Fill(vProton.Vector.P(),vProton.Vector.Theta()*TMath::RadToDeg(),w);

					EBetaProton->Fill(vProton.Vector.P(),vProton.beta,w);



					//vertex analysis
					EvertexTimePair->Fill(vElectron.vt-vPositron.vt,w);
					vertexTimePairMass->Fill((vPositron.Vector+vElectron.Vector).M(),vElectron.vt-vPositron.vt,w);
					if(true /*abs(vProton.chi2)<3. && vProton.Vector.P()<0.5*/){
						EvertexTimeP->Fill(vProton.vt-STT,w);

						if((vProton.status/1000)==4)EvertexTimePPCD->Fill(vProton.Vector.P(),/*elecVT-*/vProton.vt-STT,w);
						if((vProton.status/1000)==2)EvertexTimePPFD->Fill(vProton.Vector.P(),/*elecVT-*/vProton.vt-STT,w);
					}

					Pproton->Fill(vProton.chi2,w);


					

					EvertexPair->Fill(vElectron.vertex.z-vPositron.vertex.z,w);
					EvertexElecP->Fill(vElectron.vertex.z-vProton.vertex.z,w);
					EvertexPosiP->Fill(vPositron.vertex.z-vProton.vertex.z,w);

					ThetaPhi cm;
					cm=CM(vElectron,vPositron,vProton);

					

					double qp=(vPositron.Vector+vElectron.Vector).M();
					double b=2*(vElectron.Vector-vPositron.Vector)*(vRestProton-vProton.Vector);
					
					
					

					//double MaxAcc;
					//double MinAcc;

					double x,x1;

					int bin_phi2;

					//if(bint(t,tbins)==1)cout<<"t "<<t<<" "<<bint(t,tbins)<<" bin "<<phi_bin<<" phi "<<phi_<<" cm phi "<<cm.phi<<" theta "<<theta_<<" theta bin "<< (Acc[bint(t,tbins)]->GetYaxis()->FindBin(cm.theta))-1 <<" acc Value "<<AccValue<<endl;

					/*if(bint(t,tbins)==0){
						MaxAcc0->GetPoint(phi_bin,x,MaxAcc);
						MinAcc0->GetPoint(phi_bin,x1,MinAcc);
						//cout<<bint(t,tbins)<<" bin "<<phi_bin<<" min Acc "<<MinAcc<<" max Acc "<<MaxAcc<<endl;
					}

					if(bint(t,tbins)==1){
						MaxAcc1->GetPoint(phi_bin,x,MaxAcc);
						MinAcc1->GetPoint(phi_bin,x1,MinAcc);
						//cout<<bint(t,tbins)<<" bin "<<phi_bin<<" min Acc "<<MinAcc<<" max Acc "<<MaxAcc<<endl;
					}

					if(bint(t,tbins)==2){
						MaxAcc2->GetPoint(phi_bin,x,MaxAcc);
						MinAcc2->GetPoint(phi_bin,x1,MinAcc);
						//cout<<bint(t,tbins)<<" bin "<<phi_bin<<" min Acc "<<MinAcc<<" max Acc "<<MaxAcc<<endl;
					}

					if(bint(t,tbins)==3){
						MaxAcc3->GetPoint(phi_bin,x,MaxAcc);
						MinAcc3->GetPoint(phi_bin,x1,MinAcc);
						//cout<<bint(t,tbins)<<" bin "<<phi_bin<<" min Acc "<<MinAcc<<" max Acc "<<MaxAcc<<endl;
					}
*/



					if( 
					  qp2>1.5*1.5
					  && qp2<9.

						
					  ){
					 double u = (vPhoton-vProton.Vector).M2();
					Udistribution->Fill(-u);
					UdistributionVsTheta->Fill(-t,vProton.Vector.Theta()*TMath::RadToDeg());
					//cout<<u<<endl;
					}

					
					
					double phiEssai=true;
					if(phiEssai){
						if( 
						  qp2>0.7*0.7
						  && qp2<1.3*1.3
						  && Epho <6.
						  && Epho>5.
						  ){
					
					
					int binningEg = 0;
					
					int binningt;
					if(-t<0.15){binningt=0;}
					else if(-t>0.8){binningt=3;}
					else binningt = (AcceptancetM2forbinning[0]->GetXaxis()->FindBin(-t))-1;
					int binningMass = 0;
					
				
					int binningtPhysics = (PhysicstM2forbinningRafo->GetXaxis()->FindBin(-t))-1;
					
						
					int phi_bin = (Acc[0]->GetXaxis()->FindBin(cm.phi))-1;
					
					double phi_ = cm.phi;
					double theta_ = cm.theta;

					double CorrVolume1 = volume1[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg];	 
					double CorrVolume2 = volume2[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg];
					
					double AccValue = Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetBinContent((Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetXaxis()->FindBin(cm.phi)),(Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetYaxis()->FindBin(cm.theta)));

					double AccError = Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetBinError(Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetXaxis()->FindBin(cm.phi),Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetYaxis()->FindBin(cm.theta));
					
					
							if(phi_>-40. && phi_<40. && theta_<80. && theta_>50.){
							posQ2Phi->Fill(sqrt(qp2),w/(AccValue*CorrVolume1));
							}
					
					
							if( theta_<130. && theta_>100. && (phi_>140. || phi_<-140.) ){
							negQ2Phi->Fill(sqrt(qp2),w/(AccValue*CorrVolume2));
							}
					
						}
					}
					
					
					
					
					
					if(Rafo){
					if( Epho<5.
					  && Epho>2.
					  && -t<0.8
					  && -t>0.15
					  && qp2>1.1*1.1
					  && qp2<1.7*1.7
					  ){
					
					
					int binningEg = 0;
					int binningt = (AcceptancetM2forbinning[0]->GetXaxis()->FindBin(-t))-1;
					int binningMass = 0;
					
				
					int binningtPhysics = (PhysicstM2forbinningRafo->GetXaxis()->FindBin(-t))-1;
					
						
					int phi_bin = (Acc[0]->GetXaxis()->FindBin(cm.phi))-1;
					
					double phi_ = cm.phi;
					double theta_ = cm.theta;


					double L0=qp2*qp2*sin(theta_*TMath::DegToRad())*sin(theta_*TMath::DegToRad())/4.;
					double L=(((qp2-t)*(qp2-t))-(b*b))/4.;
					//cout<<"L0 "<<L0<<" L "<<L<<" ratio "<<L/L0<<endl;
					

					double AccValue = Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetBinContent((Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetXaxis()->FindBin(cm.phi)),(Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetYaxis()->FindBin(cm.theta)));

					double AccError = Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetBinError(Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetXaxis()->FindBin(cm.phi),Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetYaxis()->FindBin(cm.theta));
					
					
					
					
						if(AccValue>0.05 && (AccError/AccValue)<0.5){
							if(cm.theta>45. && cm.theta<135.){
							//cout<<t<<"  "<<binningtPhysics<<endl;
								RratioNumRafo[binningtPhysics]->Fill(cm.phi,(L/L0)*(w/AccValue));
								RratioDenomRafo[binningtPhysics]->Fill(cm.phi,(L/L0)*(w/AccValue));
								tmeanRratioRafo->Fill(binningtPhysics,-t*(w/AccValue));
								tmeanRratioNRafo->Fill(binningtPhysics,(w/AccValue));
							}
						}
					
					}
					}

					if( Epho<10.6
					  && Epho>2.
					  && -t<0.8
					  && -t>0.15
					  && qp2>0.6*0.6
					  && qp2<9.
					  ){
					
					
					
					}
					
					
					if( Epho<10.6
					  && Epho>4.
					  && -t<0.8
					  && -t>0.15
					  && qp2>1.5*1.5
					  && qp2<9.
					  //&& abs(Q2)<0.05
							//&& theta_<(MaxAcc-2.5)
							//&& theta_>(MinAcc+2.5)
					  ){
					 //cout<<"dafuCKKKKKK"<<endl;
					 
					
					/*EgVST->Fill(Epho,-t);
					MVST->Fill(sqrt(qp2),-t);
					MVSEg->Fill(sqrt(qp2),Epho);*/
					EPhiVSPElectron-> Fill(vElectron.Vector.Phi()*TMath::RadToDeg(),vElectron.Vector.P(),w);
					EPhiVSPPositron-> Fill(vPositron.Vector.Phi()*TMath::RadToDeg(),vPositron.Vector.P(),w);
					EPhiVSPProton-> Fill(vProton.Vector.Phi()*TMath::RadToDeg(),vProton.Vector.P(),w);						
	
					Chi2Electron->Fill(vElectron.chi2,w);
					Chi2Positron->Fill(vPositron.chi2,w);
					//double chi2WithMomCor = 
					double protonVT1=VT(vProton.Vector.P(),mp,vProton.Scintillator.t,vProton.Scintillator.path,STT);
					if(vProton.status>4000)Chi2Proton->Fill(vProton.chi2/*(protonVT1-(-1.*vProton.vertex.z+vElectron.vertex.z)/29.92)/0.1*/,w);
					if(vProton.status<4000)Chi2ProtonFD->Fill(vProton.chi2,w);
					EEnergyConservationEvents->Fill(vProton.Vector.Pz()+vPositron.Vector.Pz()+vElectron.Vector.Pz()-vPhoton.E(),w);
					tCheckEvents->Fill(tFromPhoton+t,w);
	
	
					 EhistQ21->Fill(Q2,w);
					 
					 
					PPosiPElec->Fill(vPositron.Vector.P(),vElectron.Vector.P());
					EgVSPeletron->Fill(vElectron.Vector.P(),Epho,w);
					tVSPeletron->Fill(vElectron.Vector.P(),-t,w);
				
					TCSEThetaVSPElectron-> Fill(vElectron.Vector.P(),vElectron.Vector.Theta()*TMath::RadToDeg(),w);
					TCSEThetaVSPPositron-> Fill(vPositron.Vector.P(),vPositron.Vector.Theta()*TMath::RadToDeg(),w);
					/*if(vProton.status>4000)*/TCSEThetaVSPProton-> Fill(vProton.Vector.P(),vProton.Vector.Theta()*TMath::RadToDeg(),w);
						
					
					int binningEg = (AcceptanceEg->GetXaxis()->FindBin(Epho))-1;
					int binningt = (AcceptancetM2forbinning[binningEg]->GetXaxis()->FindBin(-t))-1;
					int binningMass = (AcceptancetM2forbinning[binningEg]->GetYaxis()->FindBin(qp2))-1;
					
					
					int binningEgPhysics = (PhysicsEg1->GetXaxis()->FindBin(Epho))-1;
					int binningtPhysics = (PhysicstM2forbinning->GetXaxis()->FindBin(-t))-1;
					
					int binningtPhysicsHighMass = (pos2->GetXaxis()->FindBin(-t))-1;
					
					int binningMassPhysics = (PhysicstM2forbinning->GetYaxis()->FindBin(sqrt(qp2)))-1;
					
					//int binningMass0 = (AcceptancetM2forbinning[0]->GetYaxis()->FindBin(qp2))-1;
					
						
					int phi_bin = (Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetXaxis()->FindBin(cm.phi))-1;
					


					//double phi_ = Acc[bint(t,tbins)]->GetXaxis()->GetBinCenter(phi_bin+1);
					//double theta_ = Acc[bint(t,tbins)]->GetYaxis()->GetBinCenter(Acc[bint(t,tbins)]->GetYaxis()->FindBin(cm.theta));



					double phi_ = cm.phi;
					double theta_ = cm.theta;

					double s=(vRestProton+vPhoton).M2();

					double xi = qp2/(2*(s-mp*mp)-qp2);
					//cout<<xi<<endl;
					int binningXiPhysics = (PhysicsXi->GetXaxis()->FindBin(xi))-1;

					double L0=qp2*qp2*sin(theta_*TMath::DegToRad())*sin(theta_*TMath::DegToRad())/4.;
					double L=(((qp2-t)*(qp2-t))-(b*b))/4.;
					//cout<<"L0 "<<L0<<" L "<<L<<" ratio "<<L/L0<<endl;
					
					
					

					//w=w*(1+0.1*polarization*(L0/L)*( 1.+cos(theta_*TMath::DegToRad())*cos(theta_*TMath::DegToRad()) )*(1./sin(theta_*TMath::DegToRad()))*sin(cm.phi*TMath::DegToRad()));


					//w=w*(1-0.0*cos(cm.phi*TMath::DegToRad()));


					double AccValue = Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetBinContent((Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetXaxis()->FindBin(cm.phi)),(Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetYaxis()->FindBin(cm.theta)));

					double AccError = Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetBinError(Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetXaxis()->FindBin(cm.phi),Acc[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg]->GetYaxis()->FindBin(cm.theta));
					 
					 
					  if(MLPscore<0.6 &&vPositron.Vector.P()>4. && xi<0.2 && xi>0.1 && qp2>4. )cout<<"t "<<(-t)<<" QP2 "<<qp2<<" mlp score "<<MLPscore<<"  Accvalue "<<AccValue<<endl;
					 //cout<<nbEvent<<" "<<AccValue<<endl;
					// cout<<binningEg<<" "<<binningt<<" "<<binningMass<<" "<<phi_bin<<endl;
					 //cout<<binningEgPhysics<<" "<<binningtPhysics<<" "<<binningMassPhysics<<" "<<phi_bin<<endl;
					  nEventTCS+=w;
						denom+=(w*w);
						/*if(qp2>4.)*/
						
						FinalEventst->Fill(-t,w);
						FinalEventsPhi->Fill(vProton.Vector.Phi()*TMath::RadToDeg(),w);
						FinalEventsTheta->Fill(vProton.Vector.Theta()*TMath::RadToDeg(),w);
						FinalEventsP->Fill(vProton.Vector.P(),w);
						FinalEventsEg->Fill(Epho,w);
						FinalEventsM->Fill(sqrt(qp2),w);
						if(vProton.status>4000)xihist1->Fill(xi,w);
						xihist->Fill(xi,qp2,w);
						if(sqrt(qp2)>1.5 && sqrt(qp2)<2.0)thist->Fill(xi,-t,w);
						if(sqrt(qp2)>2.0 && sqrt(qp2)<3.0)thist1->Fill(xi,-t,w);
	
						/*if(sqrt(qp2)>1.5 && sqrt(qp2)<2.0 && xi>0.11 && x<0.2){
							sum1+=w;
							sumw21+=(w*w);
						}
						
						if(sqrt(qp2)>2. && sqrt(qp2)<3.0 && xi>0.11 && x<0.2){
							sum2+=w;
							sumw22+=(w*w);
						}*/
	
						//cout<<L<<" "<<L0<<" "<<AccValue<<endl;
						if((vProton.status/1000)==4)nCD++;
						if((vProton.status/1000)==2)nFD++;
						
						if((L/L0)<0.)cout<<(L/L0)<<endl;
						
						//sumPhi[phi_bin][binningt]+=(L/L0)*(w/AccValue);
						//sumPhi[binPhi][binningt]+=(w/AccValue);
						//ErrorSsumPhi[phi_bin][binningt]+=(((w*w)/(AccValue*AccValue))*(L/L0)*(L/L0));//*((1./AccValue)+(w*(L/L0))*AccError);
						ErrorAcceptance[phi_bin][binningt]+=0.0;
					
						double polaT =polarizationTransfer(10.6, Epho, vPhoton.Theta());
						//if(Epho>9.)cout<<Epho<<" "<<polarizationTransfer(10.6, Epho)<<endl;
						
						
						//cout<<"dafuCKKKKKK 3 "<<binningXiPhysics<<endl;
					//cout<<"dafuCKKKKKK 2"<<endl;
						if(AccValue>0.05 && (AccError/AccValue)<0.5){
						
						if(cm.theta>45. && cm.theta<135.){
						
						if(sqrt(qp2)>2.){
						
						RratioNumHighMass[binningtPhysicsHighMass]->Fill(cm.phi,(L/L0)*(w/AccValue));
						RratioDenomHighMass[binningtPhysicsHighMass]->Fill(cm.phi,(L/L0)*(w/AccValue));
						tmeanRratioHighMass->Fill(binningtPhysicsHighMass,-t*(w/AccValue));
						tmeanRratioNHighMass->Fill(binningtPhysicsHighMass,(w/AccValue));
						
						RratioHighMass_Eg->Fill(Epho,(w/AccValue));
						RratioHighMass_Q2->Fill(qp2,(w/AccValue));
						
						}
						
						
						RratioNum[binningtPhysics]->Fill(cm.phi,(L/L0)*(w/AccValue));
						RratioDenom[binningtPhysics]->Fill(cm.phi,(L/L0)*(w/AccValue));
						tmeanRratio->Fill(binningtPhysics,-t*(w/AccValue));
						tmeanRratioN->Fill(binningtPhysics,(w/AccValue));
						
						RratioNumXi[binningXiPhysics]->Fill(cm.phi,(L/L0)*(w/AccValue));
						RratioDenomXi[binningXiPhysics]->Fill(cm.phi,(L/L0)*(w/AccValue));
						XimeanRratio->Fill(binningXiPhysics,xi*(w/AccValue));
						XimeanRratioN->Fill(binningXiPhysics,(w/AccValue));
						
						Rratio_Eg->Fill(Epho,(w/AccValue));
						Rratio_Q2->Fill(qp2,(w/AccValue));
						Rratio_t->Fill(-t,(w/AccValue));
						
						}
						//cout<<"dafuCKKKKKK 2"<<endl;
						PolarizationTransfer->Fill(Epho,polaT);//Epho,polarizationTransfer(10.6, Epho,vPhoton.Theta()));
						if(polarization==1)pospola->Fill(cm.phi,w/(AccValue*polaT));
						if(polarization==-1)negpola->Fill(cm.phi,w/(AccValue*polaT));
						
						double PhiEssai = cm.phi;
						if(PhiEssai<0.0)PhiEssai=PhiEssai+360.;
						
						if(polarization==1)BSApos[binningtPhysics]->Fill(PhiEssai,w/(AccValue*polaT));
						if(polarization==-1)BSAneg[binningtPhysics]->Fill(PhiEssai,w/(AccValue*polaT));
						
						if(polarization==1)BSAposExclusyst[0]->Fill(PhiEssai,w/(AccValue*polaT));
						if(polarization==-1)BSAnegExclusyst[0]->Fill(PhiEssai,w/(AccValue*polaT));
						
						/*if(polarization==1)BSApos[binningtPhysics]->Fill(cm.phi,w/(AccValue*polaT));
						if(polarization==-1)BSAneg[binningtPhysics]->Fill(cm.phi,w/(AccValue*polaT));*/
						
						if(polarization==1)BSAposXi[binningXiPhysics]->Fill(cm.phi,w/(AccValue*polaT));
						if(polarization==-1)BSAnegXi[binningXiPhysics]->Fill(cm.phi,w/(AccValue*polaT));
						
						
						if(polarization==1)BSAposM[binningMassPhysics]->Fill(cm.phi,w/(AccValue*polaT));
						if(polarization==-1)BSAnegM[binningMassPhysics]->Fill(cm.phi,w/(AccValue*polaT));
						
						tmeanBSA->Fill(binningtPhysics,-t*(w/AccValue));
						tmeanBSAN->Fill(binningtPhysics,(w/AccValue));
						
						XimeanBSA->Fill(binningXiPhysics,xi*(w/AccValue));
						XimeanBSAN->Fill(binningXiPhysics,(w/AccValue));
						
						MmeanBSA->Fill(binningMassPhysics,sqrt(qp2)*(w/AccValue));
						MmeanBSAN->Fill(binningMassPhysics,(w/AccValue));
						
						Q2meanBSA->Fill(sqrt(qp2),(w/AccValue));
						
						BSA_theta->Fill(theta_,(w/AccValue));
						
						BSA_Eg->Fill(Epho,(w/AccValue));
						BSA_Q2->Fill(sqrt(qp2),(w/AccValue));
						BSA_t->Fill(-t,(w/AccValue));
						
						}
					
						if(AccValue>0.05 && (AccError/AccValue)<0.5 /*&& AccValue>0.05*/){
						
						
						double CorrVolume1 = volume1[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg];	 
						double CorrVolume2 = volume2[binningt+nbBinsInT*binningMass+nbBinsInT*nbBinsInM*binningEg];
						//cout<<binningEg<<" "<<binningMass<<" "<<binningt<<" "<<CorrVolume1<<" "<<CorrVolume2<<endl;
						//if(CorrVolume1<0.05)CorrVolume1=1.;
						//if(CorrVolume2<0.05)CorrVolume2=1.;
						//if(theta_<50.)PPosiPElec1->Fill(vPositron.Vector.P(),vElectron.Vector.P(),1);
						//if(theta_>50. && theta_<60.)PPosiPElec1->Fill(vPositron.Vector.P(),vElectron.Vector.P(),-1);
						if(phi_>-40. && phi_<40. && theta_<80. && theta_>50.){
						
							if(binningXiPhysics==0)PPosiPElec1->Fill(vPositron.Vector.P(),vElectron.Vector.P(),1);
							
							EgVSQ2->Fill(Epho,sqrt(qp2),w/(AccValue*CorrVolume1));
							
							pos->Fill(-t,w/(AccValue*CorrVolume1));
							posExclusyst->Fill(-t,w/(AccValue*CorrVolume1));
							tmean->Fill(binningtPhysics,-t*(w/(AccValue*CorrVolume1)));
							tmeanN->Fill(binningtPhysics,(w/(AccValue*CorrVolume1)));
							AFBt0_Eg->Fill(Epho,w/(AccValue*CorrVolume1));
							AFBt0_Q2->Fill(sqrt(qp2),w/(AccValue*CorrVolume1));
							
							posEg->Fill(Epho,w/(AccValue*CorrVolume1));
							Egmean->Fill(binningEgPhysics,Epho*(w/(AccValue*CorrVolume1)));
							EgmeanN->Fill(binningEgPhysics,(w/(AccValue*CorrVolume1)));
							AFBEg0_t->Fill(-t,w/(AccValue*CorrVolume1));
							AFBEg0_Q2->Fill(sqrt(qp2),w/(AccValue*CorrVolume1));
							
							posQ2->Fill(sqrt(qp2),w/(AccValue*CorrVolume1));
							Q2mean->Fill(binningMassPhysics,sqrt(qp2)*(w/(AccValue*CorrVolume1)));
							Q2meanN->Fill(binningMassPhysics,(w/(AccValue*CorrVolume1)));
							AFBQ20_Eg->Fill(Epho,w/(AccValue*CorrVolume1));
							AFBQ20_t->Fill(-t,w/(AccValue*CorrVolume1));
							
							posXi->Fill(xi,w/(AccValue*CorrVolume1));
							Ximean->Fill(binningXiPhysics,xi*(w/(AccValue*CorrVolume1)));
							XimeanN->Fill(binningXiPhysics,(w/(AccValue*CorrVolume1)));
							AFBXi0_Q2->Fill(sqrt(qp2),w/(AccValue*CorrVolume1));
							AFBXi0_Eg->Fill(Epho,w/(AccValue*CorrVolume1));
							AFBXi0_t->Fill(-t,w/(AccValue*CorrVolume1));
							
							if(  qp2<4.){
							
							
							
							pos1->Fill(-t,w/(AccValue*CorrVolume1));
							tmean1->Fill(binningtPhysics,-t*(w/(AccValue*CorrVolume1)));
							tmeanN1->Fill(binningtPhysics,(w/(AccValue*CorrVolume1)));
							AFBt1_Eg->Fill(Epho,w/(AccValue*CorrVolume1));
							AFBt1_Q2->Fill(sqrt(qp2),w/(AccValue*CorrVolume1));
							
							sum1+=w;
							sumw21+=(w*w);
							}
							
							if(/*xi<0.2 && xi>0.1 &&*/ qp2>4.){
							//if(-t<0.25)cout<<"EVENT"<<endl;
							
							if(binningtPhysicsHighMass==0)cout<<"EVENT in bin t0 high mass "<<(-t)<<" P proton "<<vProton.Vector.P()<<" Theta proton "<<vProton.Vector.Theta()*TMath::RadToDeg()<<" beta proton "<<vProton.Vector.Beta()<<endl;
							//if(binningtPhysicsHighMass==2)cout<<"EVENT in bin "<<(-t)<<" Acc "<<AccValue*CorrVolume1<<endl;
							
							CheckTbinning->Fill(-t);
							pos2->Fill(-t,w/(AccValue*CorrVolume1));
							pos2Exclusyst->Fill(-t,w/(AccValue*CorrVolume1));
							tmean2->Fill(binningtPhysicsHighMass,-t*(w/(AccValue*CorrVolume1)));
							tmeanN2->Fill(binningtPhysicsHighMass,(w/(AccValue*CorrVolume1)));
							AFBt2_Eg->Fill(Epho,w/(AccValue*CorrVolume1));
							AFBt2_Q2->Fill(sqrt(qp2),w/(AccValue*CorrVolume1));
							
							sum2+=w;
							sumw22+=(w*w);
							}
							
							if(qp2>4.){
							postHighMass->Fill(-t,w/(AccValue*CorrVolume1));
							
							
							
							posXiHighMass->Fill(xi,w/(AccValue*CorrVolume1));
							posEgHighMass->Fill(Epho,w/(AccValue*CorrVolume1));
							tmeanHighMass->Fill(binningtPhysics,-t*(w/(AccValue*CorrVolume1)));
							tmeanNHighMass->Fill(binningtPhysics,(w/(AccValue*CorrVolume1)));
							EgmeanHighMass->Fill(binningEgPhysics,Epho*(w/(AccValue*CorrVolume1)));
							EgmeanNHighMass->Fill(binningEgPhysics,(w/(AccValue*CorrVolume1)));
							XimeanHighMass->Fill(binningXiPhysics,xi*(w/(AccValue*CorrVolume1)));
							XimeanNHighMass->Fill(binningXiPhysics,(w/(AccValue*CorrVolume1)));
							}
							
						}
						
						if( theta_<130. && theta_>100. && (phi_>140. || phi_<-140.) ){
						if(binningXiPhysics==0)PPosiPElec1->Fill(vPositron.Vector.P(),vElectron.Vector.P(),-1);
						
						
							EgVSQ2->Fill(Epho,sqrt(qp2),w/(AccValue*CorrVolume2));
							neg->Fill(-t,w/(AccValue*CorrVolume2));
							negExclusyst->Fill(-t,w/(AccValue*CorrVolume2));
							tmean->Fill(binningtPhysics,-t*(w/(AccValue*CorrVolume2)));
							tmeanN->Fill(binningtPhysics,(w/(AccValue*CorrVolume2)));
							AFBt0_Eg->Fill(Epho,w/(AccValue*CorrVolume2));
							AFBt0_Q2->Fill(sqrt(qp2),w/(AccValue*CorrVolume2));
							
							negEg->Fill(Epho,w/(AccValue*CorrVolume2));
							Egmean->Fill(binningEgPhysics,Epho*(w/(AccValue*CorrVolume2)));
							EgmeanN->Fill(binningEgPhysics,(w/(AccValue*CorrVolume2)));
							AFBEg0_t->Fill(-t,w/(AccValue*CorrVolume2));
							AFBEg0_Q2->Fill(sqrt(qp2),w/(AccValue*CorrVolume2));
							
							
							negQ2->Fill(sqrt(qp2),w/(AccValue*CorrVolume2));
							Q2mean->Fill(binningMassPhysics,sqrt(qp2)*(w/(AccValue*CorrVolume2)));
							Q2meanN->Fill(binningMassPhysics,(w/(AccValue*CorrVolume2)));
							AFBQ20_Eg->Fill(Epho,w/(AccValue*CorrVolume2));
							AFBQ20_t->Fill(-t,w/(AccValue*CorrVolume2));
							
							
							negXi->Fill(xi,w/(AccValue*CorrVolume2));
							Ximean->Fill(binningXiPhysics,xi*(w/(AccValue*CorrVolume2)));
							XimeanN->Fill(binningXiPhysics,(w/(AccValue*CorrVolume2)));
							AFBXi0_Q2->Fill(sqrt(qp2),w/(AccValue*CorrVolume2));
							AFBXi0_Eg->Fill(Epho,w/(AccValue*CorrVolume2));
							AFBXi0_t->Fill(-t,w/(AccValue*CorrVolume2));
							
							
							if( qp2<4.){
							neg1->Fill(-t,w/(AccValue*CorrVolume2));
							tmean1->Fill(binningtPhysics,-t*(w/(AccValue*CorrVolume2)));
							tmeanN1->Fill(binningtPhysics,(w/(AccValue*CorrVolume2)));
							
							AFBt1_Eg->Fill(Epho,w/(AccValue*CorrVolume2));
							AFBt1_Q2->Fill(sqrt(qp2),w/(AccValue*CorrVolume2));
							
							sum1+=w;
							sumw21+=(w*w);
							}
							
							if(/*xi<0.2 && xi>0.1 &&*/ qp2>4.){
							//if(-t<0.25)cout<<"EVENT"<<endl;
							
							if(binningtPhysicsHighMass==0)cout<<"EVENT in bin t0 high mass "<<(-t)<<" P proton "<<vProton.Vector.P()<<" Theta proton "<<vProton.Vector.Theta()*TMath::RadToDeg()<<" beta proton "<<vProton.Vector.Beta()<<endl;
							
							CheckTbinning->Fill(-t);
							neg2->Fill(-t,w/(AccValue*CorrVolume2));
							neg2Exclusyst->Fill(-t,w/(AccValue*CorrVolume2));
							tmean2->Fill(binningtPhysicsHighMass,-t*(w/(AccValue*CorrVolume2)));
							tmeanN2->Fill(binningtPhysicsHighMass,(w/(AccValue*CorrVolume2)));
							
							AFBt2_Eg->Fill(Epho,w/(AccValue*CorrVolume2));
							AFBt2_Q2->Fill(sqrt(qp2),w/(AccValue*CorrVolume2));
							
							sum2+=w;
							sumw22+=(w*w);
							}
							
							if(qp2>4.){
							
							
							
							negtHighMass->Fill(-t,w/(AccValue*CorrVolume2));
							negXiHighMass->Fill(xi,w/(AccValue*CorrVolume2));
							negEgHighMass->Fill(Epho,w/(AccValue*CorrVolume2));
							tmeanHighMass->Fill(binningtPhysics,-t*(w/(AccValue*CorrVolume2)));
							tmeanNHighMass->Fill(binningtPhysics,(w/(AccValue*CorrVolume2)));
							EgmeanHighMass->Fill(binningEgPhysics,Epho*(w/(AccValue*CorrVolume2)));
							EgmeanNHighMass->Fill(binningEgPhysics,(w/(AccValue*CorrVolume2)));
							XimeanHighMass->Fill(binningXiPhysics,xi*(w/(AccValue*CorrVolume2)));
							XimeanNHighMass->Fill(binningXiPhysics,(w/(AccValue*CorrVolume2)));
							}
						}
						}
						
						if( /*(phi_>-50. && phi_<50. && theta_<70. && theta_>50.) || 
						
						( theta_<130. && theta_>110. && (phi_>130. || phi_<-130.) )*/ true  ){
						if(AccValue<=0.0)cout<<"hello"<<endl;//continue;
						if(binningtPhysics==0){test++;PhiVSThetaCM->Fill(cm.phi,cm.theta,1*(w/*/AccValue*/));}//*(L/L0));}
						if(binningtPhysics==1)PhiVSThetaCM1->Fill(cm.phi,cm.theta,1*(w/*/AccValue*/));//*(L/L0));
						if(binningtPhysics==2)PhiVSThetaCM2->Fill(cm.phi,cm.theta,1*(w/*/AccValue*/));//*(L/L0));
						if(binningtPhysics==3)PhiVSThetaCM3->Fill(cm.phi,cm.theta,1*(w/*/AccValue*/));//*(L/L0));
						}
					}
				}

		}

	}

	TCanvas *canU  = new TCanvas("canU","canU",1000,1000);
	canU->Divide(2,1);
	canU->cd(1);Udistribution->Draw();
	canU->cd(2);UdistributionVsTheta->Draw("colz");
	canU->SaveAs("canU.pdf");
	

	cout<<"test "<<test<<endl;

	cout<<"bin content"<<PhiVSThetaCM->GetBinContent(10,10)<<"bin error "<<PhiVSThetaCM->GetBinError(10,10)<<" accceptance "<<Acc[0]->GetBinContent(10,10)<<endl;

	TCanvas *canl  = new TCanvas("canl","tcs",4000,2000);
	canl->Divide(4,2);
	canl->cd(1);PhiVSThetaCM->Draw("colz");
	canl->cd(2);PhiVSThetaCM1->Draw("colz");
	canl->cd(3);PhiVSThetaCM2->Draw("colz");
	canl->cd(4);PhiVSThetaCM3->Draw("colz");
	canl->SaveAs("essai.root");

	EgVSQ2->SaveAs("EgVSQ2.root");

	//get t quartiles
	int nq=3;
	Double_t xq[nq];  // position where to compute the quantiles in [0,1]
	Double_t yq[nq];  // array to contain the quantiles
	for (Int_t i=0;i<nq;i++) xq[i] = Float_t(i+1)/(nq+1);
	FinalEventst->GetQuantiles(nq,yq,xq);
	cout<<"bint limits"<<yq[0]<<" "<<yq[1]<<" "<<yq[2]<<" "<<endl;


	TH1D stylet1("amplitude","amplitude;#phi;Y",binPhi,-180,180);
	TH1D *Yc[4];
	for(int bin=0;bin<4;bin++){
		Yc[bin]= new TH1D(stylet1);
		Yc[bin]->SetName("Y");
	}

	TH1D stylet1E("Eamplitude","Error amplitude;#phi;E(Y)",binPhi,-180,180);
	TH1D *EYc[4];
	for(int bin=0;bin<4;bin++){
		EYc[bin]= new TH1D(stylet1E);
		EYc[bin]->SetName("EY");
	}
/*
	//Calcul assymetry
	cout<<"Calcul"<<endl;
	double binsInT[4]={0.15+(tbins[0]-0.15)/2,tbins[0]+(tbins[1]-tbins[0])/2,tbins[1]+(tbins[2]-tbins[1])/2,tbins[2]+(-tbins[2]+0.8)/2};
	double widthT[4]={(tbins[0])/2,(tbins[1]-tbins[0])/2,(tbins[2]-tbins[1])/2,(-tbins[2]+1.)/2};	
	for(int i=0;i<tbin;i++){		
		cout<<"new T bin"<<tbin<<endl;
		double R0=0.0;
		double ER02=0.0;
		double Rc=0.0;
		double ERc2=0.0;
		double R=0.0;
		double ER=0.0; //not squared
		double widtht=0.0;
		double errorR=0.0;
		double valuet=0.0;

		for(int j=0;j<binPhi;j++){
			double phi1=((360./binPhi)*j)+(180./binPhi)-180;
			double ycVal=sumPhi[j][i];
			double ycE=sqrt(abs(ErrorSsumPhi[j][i]));
			Yc[i]->Fill(phi1,ycVal);
			EYc[i]->Fill(phi1,ycE);
			cout<<sumPhi[j][i]<<endl;
			R0+=sumPhi[j][i];
			Rc+=sumPhi[j][i]*cos(phi1*TMath::DegToRad());
			ER02+=abs(ErrorSsumPhi[j][i])+ abs(ErrorAcceptance[j][i]);
			ERc2+=(abs(ErrorSsumPhi[j][i])+ abs(ErrorAcceptance[j][i]))*cos(phi1*TMath::DegToRad())*cos(phi1*TMath::DegToRad());
		}	
		R=Rc/R0;
		ER=abs(R)*sqrt(abs(ER02*(1./R0)*(1./R0))+abs(ERc2*(1./Rc)*(1./Rc)));
		cout<<"R "<<R<<" R0  "<<R0<<" Rc  "<<Rc<<" ER "<<ER<<" ER0  "<<sqrt(ER02)<<" ERc  "<<sqrt(ERc2)<<endl;
		valuet=binsInT[i];
		widtht=widthT[i];

		Assym->SetPoint(i,valuet,R);
		//Assym->SetPointError(i,0.0,ER);
		
		cout<<"R "<<R<<" bin "<<i<<endl;
	}


	for(int i=0;i<tbin;i++){			
		Assym->SetPointError(i,0.0,smear(*Yc[i], *EYc[i],i+1));
	}
	cout<<smear(*Yc[0], *EYc[0],1)<<endl;
	cout<<smear(*Yc[1], *EYc[1],2)<<endl;
	cout<<smear(*Yc[2], *EYc[2],3)<<endl;
	cout<<smear(*Yc[3], *EYc[3],4)<<endl;
	
*/


	//highmass AFB
	/*
	TH1F *num= (TH1F*)postHighMass->Clone("postHighMass");
	TH1F *denom1= (TH1F*)postHighMass->Clone("postHighMass");
	num->Sumw2();
	denom1->Sumw2();
	
	num->Add(neg,-1);
	denom1->Add(neg);
	num->Divide(denom1);
	TGraphAsymmErrors* AFBt= new TGraphAsymmErrors(4);
	cout<<"going to do afbt"<<endl;
	for(int i=0;i<4;i++){
		double bincenter=(tmean->GetBinContent(i+1))/(tmeanN->GetBinContent(i+1));
		cout<<"AFBt "<<bincenter<<" "<<(num->GetBinContent(i+1))<<endl;
		AFBt->SetPoint(i,bincenter,(num->GetBinContent(i+1)));
		AFBt->SetPointError(i,(bincenter-t1[i]),(t1[i+1]-bincenter),(num->GetBinError(i+1)),(num->GetBinError(i+1)));
		cout<<t1[i]<<" "<<endl;
	}
	AFBt->SetMarkerSize(5);
   	AFBt->SetMarkerColor(kRed);
   	AFBt->SetLineColor(kRed);
   	AFBt->GetXaxis()->SetRangeUser(0.15,0.8);
	AFBt->SetTitle(";-t;A_{FB}");
	
	TH1F *numXi= (TH1F*)posXi->Clone("posXi");
	TH1F *denom1Xi= (TH1F*)posXi->Clone("posXi");
	numXi->Sumw2();
	denom1Xi->Sumw2();
	
	numXi->Add(negXi,-1);
	denom1Xi->Add(negXi);
	numXi->Divide(denom1Xi);
	TGraphAsymmErrors* AFBXi= new TGraphAsymmErrors(3);
	for(int i=0;i<3;i++){
		double bincenter=(Ximean->GetBinContent(i+1))/(XimeanN->GetBinContent(i+1));
		AFBXi->SetPoint(i,bincenter,(numXi->GetBinContent(i+1)));
		AFBXi->SetPointError(i,bincenter-Xi3[i],Xi3[i+1]-bincenter,(numXi->GetBinError(i+1)),(numXi->GetBinError(i+1)));
	}
	AFBXi->SetMarkerSize(5);
   	AFBXi->SetMarkerColor(kOrange);
   	AFBXi->SetLineColor(kOrange);
   	//AFBXi->GetXaxis()->SetRangeUser(4.,10.);
	AFBXi->SetTitle(";#xi;A_{FB}");
	
	
	
	
	TH1F *numEg= (TH1F*)posEg->Clone("posEg");
	TH1F *denom1Eg= (TH1F*)posEg->Clone("posEg");
	numEg->Sumw2();
	denom1Eg->Sumw2();
	
	numEg->Add(negEg,-1);
	denom1Eg->Add(negEg);
	numEg->Divide(denom1Eg);
	TGraphAsymmErrors* AFBEg= new TGraphAsymmErrors(3);
	for(int i=0;i<3;i++){
		double bincenter=(Egmean->GetBinContent(i+1))/(EgmeanN->GetBinContent(i+1));
		AFBEg->SetPoint(i,bincenter,(numEg->GetBinContent(i+1)));
		AFBEg->SetPointError(i,bincenter-Eg1[i],Eg1[i+1]-bincenter,(numEg->GetBinError(i+1)),(numEg->GetBinError(i+1)));
	}
	AFBEg->SetMarkerSize(5);
   	AFBEg->SetMarkerColor(kGreen);
   	AFBEg->SetLineColor(kGreen);
   	AFBEg->GetXaxis()->SetRangeUser(4.,10.);
	AFBEg->SetTitle(";E_{#gamma};A_{FB}");
	
	*/
	////
	
	TH1F *numQ2Phi= (TH1F*)posQ2Phi->Clone("posQ2Phi");
	TH1F *denom1Q2Phi= (TH1F*)posQ2Phi->Clone("posQ2Phi");
	numQ2Phi->Sumw2();
	denom1Q2Phi->Sumw2();
	
	numQ2Phi->Add(negQ2Phi,-1);
	denom1Q2Phi->Add(negQ2Phi);
	numQ2Phi->Divide(denom1Q2Phi);
	TGraphAsymmErrors* AFBQ2Phi= new TGraphAsymmErrors(15);
	cout<<"going to do afbQ2Phi"<<endl;
	for(int i=0;i<15;i++){
		double bincenter=posQ2Phi->GetBinCenter(i+1);
		AFBQ2Phi->SetPoint(i,bincenter,(numQ2Phi->GetBinContent(i+1)));
		AFBQ2Phi->SetPointError(i,0.0,0.0,(numQ2Phi->GetBinError(i+1)),(numQ2Phi->GetBinError(i+1)));
		
	}
	AFBQ2Phi->SetMarkerSize(5);
   	AFBQ2Phi->SetMarkerColor(kRed);
   	AFBQ2Phi->SetLineColor(kRed);
   	AFBQ2Phi->GetXaxis()->SetRangeUser(0.7,1.3);
	AFBQ2Phi->SetTitle(";M (GeV);A_{FB}");
	AFBQ2Phi->SaveAs("AFB_phi.root");
	
	
	TCanvas *canAFBphi  = new TCanvas("canAFBphi","canAFBphi",5000,1000);
	canAFBphi->Divide(4,1);
	canAFBphi->cd(1);posQ2Phi->Draw();
	canAFBphi->cd(2);negQ2Phi->Draw();
	canAFBphi->cd(3);AFBQ2Phi->Draw();
	canAFBphi->SaveAs("canAFBphi.pdf");
	
	//


	TH1F *num= (TH1F*)pos->Clone("pos");
	TH1F *denom1= (TH1F*)pos->Clone("pos");
	num->Sumw2();
	denom1->Sumw2();
	
	num->Add(neg,-1);
	denom1->Add(neg);
	num->Divide(denom1);
	TGraphAsymmErrors* AFBt= new TGraphAsymmErrors(4);
	cout<<"going to do afbt"<<endl;
	for(int i=0;i<4;i++){
		double bincenter=(tmean->GetBinContent(i+1))/(tmeanN->GetBinContent(i+1));
		cout<<"AFBt "<<bincenter<<" "<<(num->GetBinContent(i+1))<<endl;
		AFBt->SetPoint(i,bincenter,(num->GetBinContent(i+1)));
		AFBt->SetPointError(i,(bincenter-t1[i]),(t1[i+1]-bincenter),(num->GetBinError(i+1)),(num->GetBinError(i+1)));
		cout<<t1[i]<<" "<<endl;
	}
	AFBt->SetMarkerSize(5);
   	AFBt->SetMarkerColor(kRed);
   	AFBt->SetLineColor(kRed);
   	AFBt->GetXaxis()->SetRangeUser(0.15,0.8);
	AFBt->SetTitle(";-t;A_{FB}");
	
	
	//AFB_t exclusyst
	TH1F *numExclusyst= (TH1F*)posExclusyst->Clone("posExclusyst");
	TH1F *denom1Exclusyst= (TH1F*)posExclusyst->Clone("posExclusyst");
	numExclusyst->Sumw2();
	denom1Exclusyst->Sumw2();
	
	numExclusyst->Add(negExclusyst,-1);
	denom1Exclusyst->Add(negExclusyst);
	numExclusyst->Divide(denom1Exclusyst);
	TGraphAsymmErrors* AFBtExclusyst= new TGraphAsymmErrors(1);
	cout<<"going to do afbt"<<endl;
	
		double bincenterExclusyst=0.5;
		AFBtExclusyst->SetPoint(0,bincenterExclusyst,(numExclusyst->GetBinContent(1)));
		AFBtExclusyst->SetPointError(0,0.35,0.3,(numExclusyst->GetBinError(1)),(numExclusyst->GetBinError(1)));
		
	
	AFBtExclusyst->SetMarkerSize(5);
   	AFBtExclusyst->SetMarkerColor(kRed);
   	AFBtExclusyst->SetLineColor(kRed);
   	AFBtExclusyst->GetXaxis()->SetRangeUser(0.15,0.8);
	AFBtExclusyst->SetTitle(";-t;A_{FB}");
	AFBtExclusyst->SaveAs("AFBtExclusyst.root");
	
	
	//AFB t binné
	
	TH1F *num1= (TH1F*)pos1->Clone("pos1");
	TH1F *denom11= (TH1F*)pos1->Clone("pos1");
	num1->Sumw2();
	denom11->Sumw2();
	
	num1->Add(neg1,-1);
	denom11->Add(neg1);
	num1->Divide(denom11);
	TGraphAsymmErrors* AFBt1= new TGraphAsymmErrors(4);
	cout<<"going to do afbt1"<<endl;
	for(int i=0;i<4;i++){
		double bincenter=(tmean1->GetBinContent(i+1))/(tmeanN1->GetBinContent(i+1));
		cout<<"AFBt1 "<<bincenter<<" "<<(num1->GetBinContent(i+1))<<endl;
		AFBt1->SetPoint(i,bincenter,(num1->GetBinContent(i+1)));
		AFBt1->SetPointError(i,(bincenter-t1[i]),(t1[i+1]-bincenter),(num1->GetBinError(i+1)),(num1->GetBinError(i+1)));
		cout<<t1[i]<<" "<<endl;
	}
	AFBt1->SetMarkerSize(5);
   	AFBt1->SetMarkerColor(kRed);
   	AFBt1->SetLineColor(kRed);
   	AFBt1->GetXaxis()->SetRangeUser(0.15,0.8);
   	//AFBt1->GetYaxis()->SetRangeUser(0.,0.1);
	AFBt1->SetTitle("1.5 Gev<M<2 GeV;-t;A_{FB}");
	
	
	
	TH1F *num2= (TH1F*)pos2->Clone("pos2");
	TH1F *denom12= (TH1F*)pos2->Clone("pos2");
	num2->Sumw2();
	denom12->Sumw2();
	
	
	//get t quartiles
	int na=4;
	Double_t xa[na];  // position where to compute the quantiles in [0,1]
	Double_t ya[na];  // array to contain the quantiles
	for (Int_t i=0;i<na;i++) xa[i] = Float_t(i+1)/(na+1);
	CheckTbinning->GetQuantiles(na,ya,xa);
	cout<<"bint limits"<<ya[0]<<" "<<ya[1]<<" "<<ya[2]<<" "<<endl;
	
	
	cout<<pos2->GetEntries()<<endl;
	cout<<neg2->GetEntries()<<endl;
	num2->Add(neg2,-1);
	denom12->Add(neg2);
	num2->Divide(denom12);
	TGraphAsymmErrors* AFBt2= new TGraphAsymmErrors(4);
	cout<<"going to do afbt2"<<endl;
	for(int i=0;i<4;i++){
		double bincenter=(tmean2->GetBinContent(i+1))/(tmeanN2->GetBinContent(i+1));
		cout<<"AFBt2 "<<bincenter<<" "<<(num2->GetBinContent(i+1))<<endl;
		AFBt2->SetPoint(i,bincenter,(num2->GetBinContent(i+1)));
		AFBt2->SetPointError(i,(bincenter-t1b[i]),(t1b[i+1]-bincenter),(num2->GetBinError(i+1)),(num2->GetBinError(i+1)));
		cout<<t1b[i]<<" "<<endl;
	}
	AFBt2->SetMarkerSize(5);
   	AFBt2->SetMarkerColor(kRed);
   	AFBt2->SetLineColor(kRed);
   	AFBt2->GetXaxis()->SetRangeUser(0.15,0.8);
   	//AFBt2->GetYaxis()->SetRangeUser(0.,0.1);
	AFBt2->SetTitle("2 Gev<M<3 GeV;-t;A_{FB}");
	
	//AFB t2 exclu syst
	TH1F *num2Exclusyst= (TH1F*)pos2Exclusyst->Clone("pos2Exclusyst");
	TH1F *denom12Exclusyst= (TH1F*)pos2Exclusyst->Clone("pos2Exclusyst");
	num2Exclusyst->Sumw2();
	denom12Exclusyst->Sumw2();
	num2Exclusyst->Add(neg2Exclusyst,-1);
	denom12Exclusyst->Add(neg2Exclusyst);
	num2Exclusyst->Divide(denom12Exclusyst);
	TGraphAsymmErrors* AFBt2Exclusyst= new TGraphAsymmErrors(1);
		bincenterExclusyst=0.5;
		AFBt2Exclusyst->SetPoint(0,bincenterExclusyst,(num2Exclusyst->GetBinContent(1)));
		AFBt2Exclusyst->SetPointError(0,0.35,0.3,(num2Exclusyst->GetBinError(1)),(num2Exclusyst->GetBinError(1)));
	AFBt2Exclusyst->SetMarkerSize(5);
   	AFBt2Exclusyst->SetMarkerColor(kRed);
   	AFBt2Exclusyst->SetLineColor(kRed);
   	AFBt2Exclusyst->GetXaxis()->SetRangeUser(0.15,0.8);
   	//AFBt2->GetYaxis()->SetRangeUser(0.,0.1);
	AFBt2Exclusyst->SetTitle("2 Gev<M<3 GeV;-t;A_{FB}");
	AFBt2Exclusyst->SaveAs("AFBt2Exclusyst.root");
	
	TCanvas *canAFBtbinne  = new TCanvas("canAFBtbinne","canAFBtbinne",5000,3000);
	canAFBtbinne->Divide(3,2);
	canAFBtbinne->cd(1);AFBt->Draw("AP");
	canAFBtbinne->cd(2);AFBt1->Draw("AP");
	canAFBtbinne->cd(3);AFBt2->Draw("AP");
	canAFBtbinne->cd(4);pos2->Draw("");
	canAFBtbinne->cd(5);neg2->Draw("");
	canAFBtbinne->cd(6);CheckTbinning->Draw("");
	canAFBtbinne->SaveAs("canAFBtbinne.pdf");
	AFBt1->SaveAs("AFBt1.root");
	AFBt2->SaveAs("AFBt2.root");
	
	
	//////////////
	
	
	
	
	TH1F *numXi= (TH1F*)posXi->Clone("posXi");
	TH1F *denom1Xi= (TH1F*)posXi->Clone("posXi");
	numXi->Sumw2();
	denom1Xi->Sumw2();
	
	numXi->Add(negXi,-1);
	denom1Xi->Add(negXi);
	numXi->Divide(denom1Xi);
	TGraphAsymmErrors* AFBXi= new TGraphAsymmErrors(3);
	for(int i=0;i<3;i++){
		double bincenter=(Ximean->GetBinContent(i+1))/(XimeanN->GetBinContent(i+1));
		AFBXi->SetPoint(i,bincenter,(numXi->GetBinContent(i+1)));
		AFBXi->SetPointError(i,bincenter-Xi3[i],Xi3[i+1]-bincenter,(numXi->GetBinError(i+1)),(numXi->GetBinError(i+1)));
	}
	AFBXi->SetMarkerSize(5);
   	AFBXi->SetMarkerColor(kOrange);
   	AFBXi->SetLineColor(kOrange);
   	//AFBXi->GetXaxis()->SetRangeUser(4.,10.);
	AFBXi->SetTitle(";#xi;A_{FB}");
	
	
	
	
	TH1F *numEg= (TH1F*)posEg->Clone("posEg");
	TH1F *denom1Eg= (TH1F*)posEg->Clone("posEg");
	numEg->Sumw2();
	denom1Eg->Sumw2();
	
	numEg->Add(negEg,-1);
	denom1Eg->Add(negEg);
	numEg->Divide(denom1Eg);
	TGraphAsymmErrors* AFBEg= new TGraphAsymmErrors(3);
	for(int i=0;i<3;i++){
		double bincenter=(Egmean->GetBinContent(i+1))/(EgmeanN->GetBinContent(i+1));
		AFBEg->SetPoint(i,bincenter,(numEg->GetBinContent(i+1)));
		AFBEg->SetPointError(i,bincenter-Eg1[i],Eg1[i+1]-bincenter,(numEg->GetBinError(i+1)),(numEg->GetBinError(i+1)));
	}
	AFBEg->SetMarkerSize(5);
   	AFBEg->SetMarkerColor(kGreen);
   	AFBEg->SetLineColor(kGreen);
   	AFBEg->GetXaxis()->SetRangeUser(4.,10.);
	AFBEg->SetTitle(";E_{#gamma};A_{FB}");
	
	TH1F *numQ2= (TH1F*)posQ2->Clone("posQ2");
	TH1F *denom1Q2= (TH1F*)posQ2->Clone("posQ2");
	numQ2->Sumw2();
	denom1Q2->Sumw2();
	
	numQ2->Add(negQ2,-1);
	denom1Q2->Add(negQ2);
	numQ2->Divide(denom1Q2);
	TGraphAsymmErrors* AFBQ2= new TGraphAsymmErrors(4);
	for(int i=0;i<4;i++){
		double bincenterQ2=(Q2mean->GetBinContent(i+1))/(Q2meanN->GetBinContent(i+1));
		cout<<"bincenter Q2 "<<bincenterQ2<<" "<<(Q2mean->GetBinContent(i+1))<<" "<<(Q2meanN->GetBinContent(i+1))<<endl;
		AFBQ2->SetPoint(i,bincenterQ2,(numQ2->GetBinContent(i+1)));
		AFBQ2->SetPointError(i,(bincenterQ2-M1[i]),(M1[i+1]-bincenterQ2),(numQ2->GetBinError(i+1)),(numQ2->GetBinError(i+1)));
		cout<<M1[i]<<" "<<M1[i+1]<<endl;
	}
	AFBQ2->SetMarkerSize(5);
   	AFBQ2->SetMarkerColor(kBlue);
   	AFBQ2->SetLineColor(kBlue);
   	AFBQ2->GetXaxis()->SetRangeUser(1.5,3.);
	AFBQ2->SetTitle(";Q'^{2};A_{FB}");
	
	
	
	TCanvas *canAFB  = new TCanvas("canAFB","canAFB",5000,1000);
	canAFB->Divide(4,1);
	canAFB->cd(1);pos->Draw();
	canAFB->cd(2);neg->Draw();
	canAFB->cd(3);num->Draw();
	canAFB->SaveAs("canAFB.pdf");
	//num->SaveAs("AFB.root");
	
	TCanvas *canAFBEg  = new TCanvas("canAFBEg","canAFBEg",5000,1000);
	canAFBEg->Divide(4,1);
	canAFBEg->cd(1);posEg->Draw();
	canAFBEg->cd(2);negEg->Draw();
	canAFBEg->cd(3);numEg->Draw();
	canAFBEg->cd(4);AFBEg0_t->Draw();
	canAFBEg->SaveAs("canAFBEg.pdf");
	//numEg->SaveAs("AFBEg.root");
	
	TCanvas *canAFBQ2  = new TCanvas("canAFBQ2","canAFBQ2",5000,1000);
	canAFBQ2->Divide(4,1);
	canAFBQ2->cd(1);posQ2->Draw();
	canAFBQ2->cd(2);negQ2->Draw();
	canAFBQ2->cd(3);numQ2->Draw();
	canAFBQ2->SaveAs("canAFBQ2.pdf");
	//numQ2->SaveAs("AFBQ2.root");
	
	TCanvas *canAFBXi  = new TCanvas("canAFBXi","canAFBXi",5000,1000);
	canAFBXi->Divide(4,1);
	canAFBXi->cd(1);posXi->Draw();
	canAFBXi->cd(2);negXi->Draw();
	canAFBXi->cd(3);numXi->Draw();
	canAFBXi->SaveAs("canAFBXi.pdf");
	//numQ2->SaveAs("AFBQ2.root");
	
	TCanvas *canAFBtot  = new TCanvas("canAFBtot","canAFBtot",9000,3200);
	canAFBtot->Divide(4,2);
	num->SetStats(kFALSE);
	numEg->SetStats(kFALSE);
	numQ2->SetStats(kFALSE);
	numXi->SetStats(kFALSE);
	num->SetTitle(";-t;A_{FB}");
	numEg->SetTitle(";E_{#gamma};A_{FB}");
	numQ2->SetTitle(";Q'^{2};A_{FB}");
	num->SetMarkerSize(3);
   	num->SetMarkerColor(kRed);
   	num->SetLineColor(kRed);
   	
   	numEg->SetMarkerSize(3);
   	numEg->SetMarkerColor(kGreen);
   	numEg->SetLineColor(kGreen);
   	
   	numQ2->SetMarkerSize(3);
   	numQ2->SetMarkerColor(kBlue);
   	numQ2->SetLineColor(kBlue);
   	
	canAFBtot->cd(1);num->Draw();
	canAFBtot->cd(2);numEg->Draw();
	canAFBtot->cd(3);numQ2->Draw();
	canAFBtot->cd(5);AFBt->GetHistogram()->SetMaximum(0.5);AFBt->GetHistogram()->SetMinimum(-0.5);AFBt->Draw("A P");
	canAFBtot->cd(6);AFBEg->GetHistogram()->SetMaximum(0.5);AFBEg->GetHistogram()->SetMinimum(-0.5);AFBEg->Draw("A P");
	canAFBtot->cd(7);AFBQ2->GetHistogram()->SetMaximum(0.7);AFBQ2->GetHistogram()->SetMinimum(-0.7);AFBQ2->Draw("A P");
	canAFBtot->cd(8);AFBXi->GetHistogram()->SetMaximum(0.5);AFBXi->GetHistogram()->SetMinimum(-0.5);AFBXi->Draw("A P");
	canAFBtot->SaveAs("canAFBtot.pdf");
	AFBt->SaveAs("AFBt.root");
	AFBEg->SaveAs("AFBEg.root");
	AFBQ2->SaveAs("AFBQ2.root");
	AFBXi->SaveAs("AFBXi.root");
	
	
	double factorPola = 0.892;//0.863;//
	
	
	TH1F *numpola= (TH1F*)pospola->Clone("pospola");
	TH1F *denom1pola= (TH1F*)pospola->Clone("pospola");
	numpola->Sumw2();
	denom1pola->Sumw2();
	
	numpola->Add(negpola,-1);
	denom1pola->Add(negpola);
	numpola->Divide(denom1pola);
	numpola->Scale(1./factorPola);
	
	
	
	
	//BSA versus t
	TF1  *f2 = new TF1("f2","[0]*sin(x*[1]/[2])",0,360);
	f2->FixParameter(1,3.14159264);
	f2->FixParameter(2,180);
	//f2->SetParameter(3,0.0);
	TGraphAsymmErrors *BSAversust = new TGraphAsymmErrors(4);
	TGraphAsymmErrors *BSAversustmethod2 = new TGraphAsymmErrors(4);
	TCanvas *canpolat  = new TCanvas("canpolat","canpolat",4000,4750);
	canpolat->Divide(2,3);
	for(int tbin=0;tbin<4;tbin++){
	
		TH1F *numpolap= (TH1F*)BSApos[tbin]->Clone(Form("t bin %d pos",tbin));
		TH1F *denom1polap= (TH1F*)BSApos[tbin]->Clone(Form("t bin %d pos",tbin));
		numpolap->Sumw2();
		denom1polap->Sumw2();
		numpolap->SetTitle(Form("Bin %d;#phi;BSA",tbin));
		numpolap->Add(BSAneg[tbin],-1);
		denom1polap->Add(BSAneg[tbin]);
		numpolap->Divide(denom1polap);
		numpolap->Scale(1./factorPola);
		numpolap->GetXaxis()->SetTitle("#phi");
		gStyle->SetOptFit(1);
		gStyle->SetOptStat(0);
		canpolat->cd(tbin+1);numpolap->Fit("f2");
		
		numpolap->SaveAs(Form("FitT%d.root",tbin));
		//double tbincenter=(PhysicstM2forbinning->GetXaxis()->GetBinCenter(1+tbin));
		
		double fitresult = numpolap->GetFunction("f2")->GetParameter(0);
		double fiterror = numpolap->GetFunction("f2")->GetParError(0);

		//PairValue resultsFit= smearFit(*numpolap, tbin);
		
		//BSAversustmethod2->SetPoint(tbin,tbincenter,resultsFit.mean);
		//BSAversustmethod2->SetPointError(tbin,0.0,resultsFit.sigma);
		double tbincenter=(tmeanBSA->GetBinContent(tbin+1))/(tmeanBSAN->GetBinContent(tbin+1));
		cout<<" tbin center "<<tbincenter<<" "<<(tbincenter-t1[tbin])<<" "<<(t1[tbin+1]-tbincenter)<<" "<<fiterror<<endl;
		BSAversust->SetPoint(tbin,tbincenter,fitresult);
		BSAversust->SetPointError(tbin,(tbincenter-t1[tbin]),(t1[tbin+1]-tbincenter),fiterror,fiterror);
	}
	canpolat->cd(5);
	BSAversust->SetFillColor(kRed);
   	BSAversust->SetFillStyle(3005);
   	BSAversust->SetMarkerSize(3);
   	BSAversust->SetMarkerColor(kRed);
   	BSAversust->SetLineColor(kRed);
   	BSAversustmethod2->SetMarkerSize(3);
   	BSAversustmethod2->SetMarkerColor(kGreen);
   	BSAversustmethod2->SetLineColor(kGreen);
	BSAversust->SetTitle(";-t;BSA");BSAversust->Draw("A P");
	//BSAversustmethod2->Draw("P");
	canpolat->SaveAs("canpolat.pdf");
	canpolat->SaveAs("canpolat.root");
	BSAversust->SaveAs("BSAversust.root");
	
	
	//BSA versus t syst Exclucuts
	TGraphAsymmErrors *BSAversustExclusyst = new TGraphAsymmErrors(1);
	TCanvas *canpolatExclusyst  = new TCanvas("canpolatExclusyst","canpolatExclusyst",4000,4750);
	TH1F *numpolapExclusyst= (TH1F*)BSAposExclusyst[0]->Clone(Form("t bin %d pos Exclusyst",tbin));
	TH1F *denom1polapExclusyst= (TH1F*)BSAposExclusyst[0]->Clone(Form("t bin %d pos Exclusyst",tbin));
	numpolapExclusyst->Sumw2();
	denom1polapExclusyst->Sumw2();
	numpolapExclusyst->SetTitle(Form("Bin %d;#phi;BSA",tbin));
	numpolapExclusyst->Add(BSAnegExclusyst[0],-1);
	denom1polapExclusyst->Add(BSAnegExclusyst[0]);
	numpolapExclusyst->Divide(denom1polapExclusyst);
	numpolapExclusyst->Scale(1./factorPola);
	numpolapExclusyst->GetXaxis()->SetTitle("#phi");
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(0);
	numpolapExclusyst->Fit("f2");
		
	numpolapExclusyst->SaveAs(Form("FitT%dExclusyst.root",tbin));
		
	double fitresultExclusyst = numpolapExclusyst->GetFunction("f2")->GetParameter(0);
	double fiterrorExclusyst = numpolapExclusyst->GetFunction("f2")->GetParError(0);

	double tbincenterExclusyst=0.5;
	BSAversustExclusyst->SetPoint(0,tbincenterExclusyst,fitresultExclusyst);
	BSAversustExclusyst->SetPointError(0,0.35,0.3,0,0);
	canpolat->cd(5);
	BSAversustExclusyst->SetFillColor(kRed);
   	BSAversustExclusyst->SetFillStyle(3005);
   	BSAversustExclusyst->SetMarkerSize(3);
   	BSAversustExclusyst->SetMarkerColor(kRed);
   	BSAversustExclusyst->SetLineColor(kRed);
	BSAversustExclusyst->SetTitle(";-t;BSA");BSAversustExclusyst->Draw("A P");
	//BSAversustmethod2->Draw("P");
	canpolatExclusyst->SaveAs("canpolatExclusyst.pdf");
	canpolatExclusyst->SaveAs("canpolatExclusyst.root");
	BSAversustExclusyst->SaveAs("BSAversustExclusyst.root");
	
	
	// BSA versus xi
	TGraphAsymmErrors *BSAversusXi = new TGraphAsymmErrors(3);
	
	TGraphErrors *BSAversusXimethod2 = new TGraphErrors(3);
	TCanvas *canpolaXi  = new TCanvas("canpolaXi","canpolaXi",8000,3000);
	canpolaXi->Divide(4,2);
	for(int Xibin=0;Xibin<3;Xibin++){
	
		TH1F *numpolaXi= (TH1F*)BSAposXi[Xibin]->Clone(Form("Xi bin %d pos",Xibin));
		TH1F *denom1polaXi= (TH1F*)BSAposXi[Xibin]->Clone(Form("Xi bin %d pos",Xibin));
		numpolaXi->Sumw2();
		denom1polaXi->Sumw2();
	
		numpolaXi->Add(BSAnegXi[Xibin],-1);
		denom1polaXi->Add(BSAnegXi[Xibin]);
		numpolaXi->Divide(denom1polaXi);
		numpolaXi->Scale(1./factorPola);
	
		

		
		
		canpolaXi->cd(Xibin+1);numpolaXi->Fit("f2");//->Draw();//
		
		//double Xibincenter=(PhysicsXi->GetXaxis()->GetBinCenter(1+Xibin));
		
		double fitresultXi = numpolaXi->GetFunction("f2")->GetParameter(0);
		double fiterrorXi = numpolaXi->GetFunction("f2")->GetParError(0);
		
		
		//PairValue resultsFit= smearFit(*numpolaXi, Xibin);
		
		//BSAversusXimethod2->SetPoint(Xibin,Xibincenter,resultsFit.mean);
		//BSAversusXimethod2->SetPointError(Xibin,0.0,resultsFit.sigma);
		double Xibincenter=(XimeanBSA->GetBinContent(Xibin+1))/(XimeanBSAN->GetBinContent(Xibin+1));
		BSAversusXi->SetPoint(Xibin,Xibincenter,fitresultXi);
		BSAversusXi->SetPointError(Xibin,(Xibincenter-Xi3[Xibin]),(Xi3[Xibin+1]-Xibincenter),fiterrorXi,fiterrorXi);
		
	}
	canpolaXi->cd(5);
	BSAversusXi->SetFillColor(kRed);
   	BSAversusXi->SetFillStyle(3005);
   	BSAversusXi->SetMarkerSize(3);
   	BSAversusXi->SetMarkerColor(kRed);
   	BSAversusXi->SetLineColor(kRed);
   	BSAversusXimethod2->SetMarkerSize(3);
   	BSAversusXimethod2->SetMarkerColor(kGreen);
   	BSAversusXimethod2->SetLineColor(kGreen);
	BSAversusXi->SetTitle(";Xi;BSA");BSAversusXi->Draw("A P");
	//BSAversusXimethod2->Draw(" P");
	canpolaXi->SaveAs("canpolaXi.pdf");
	BSAversusXi->SaveAs("BSAversusXi.root");
	
	
	
	// BSA versus M
	TGraphAsymmErrors *BSAversusM = new TGraphAsymmErrors(3);
	
	TCanvas *canpolaM  = new TCanvas("canpolaM","canpolaM",8000,3000);
	canpolaM->Divide(4,2);
	for(int Mbin=0;Mbin<4;Mbin++){
	
		TH1F *numpolaM= (TH1F*)BSAposM[Mbin]->Clone(Form("M bin %d pos",Mbin));
		TH1F *denom1polaM= (TH1F*)BSAposM[Mbin]->Clone(Form("M bin %d pos",Mbin));
		numpolaM->Sumw2();
		denom1polaM->Sumw2();
	
		numpolaM->Add(BSAnegM[Mbin],-1);
		denom1polaM->Add(BSAnegM[Mbin]);
		numpolaM->Divide(denom1polaM);
		numpolaM->Scale(1./factorPola);
	
		

		
		
		canpolaM->cd(Mbin+1);numpolaM->Fit("f2");//->Draw();//
		
		//double Mbincenter=(PhysicstM2forbinning->GetYaxis()->GetBinCenter(1+Mbin));
		
		double fitresultM = numpolaM->GetFunction("f2")->GetParameter(0);
		double fiterrorM = numpolaM->GetFunction("f2")->GetParError(0);
		
		
		//PairValue resultsFit= smearFit(*numpolaXi, Xibin);
		
		//BSAversusXimethod2->SetPoint(Xibin,Xibincenter,resultsFit.mean);
		//BSAversusXimethod2->SetPointError(Xibin,0.0,resultsFit.sigma);
		double Mbincenter=(MmeanBSA->GetBinContent(Mbin+1))/(MmeanBSAN->GetBinContent(Mbin+1));
		cout<<" Mbincenter "<<Mbincenter<<" "<<(Mbincenter-M1[Mbin])<<" "<<(M1[Mbin+1]-Mbincenter)<<" "<<endl;
		BSAversusM->SetPoint(Mbin,Mbincenter,fitresultM);
		BSAversusM->SetPointError(Mbin,(Mbincenter-M1[Mbin]),(M1[Mbin+1]-Mbincenter),fiterrorM,fiterrorM);
		
	}
	canpolaM->cd(5);
	BSAversusM->SetFillColor(kRed);
   	BSAversusM->SetFillStyle(3005);
   	BSAversusM->SetMarkerSize(3);
   	BSAversusM->SetMarkerColor(kRed);
   	BSAversusM->SetLineColor(kRed);
	BSAversusM->SetTitle(";M;BSA");BSAversusM->Draw("A P");
	//BSAversusXimethod2->Draw(" P");
	canpolaM->SaveAs("canpolaM.pdf");
	BSAversusM->SaveAs("BSAversusM.root");
	
	
	//R ratio
	TF1  *f3 = new TF1("f2","cos(x*[1]/[2])",-180,180);
	f3->FixParameter(1,3.14159264);
	f3->FixParameter(2,180);
	TGraphAsymmErrors *Rratio = new TGraphAsymmErrors(4);
	TCanvas *canRatio  = new TCanvas("canRatio","canRatio",8000,4000);
	canRatio->Divide(4,3);
	for(int tbin=0;tbin<4;tbin++){
	
		TH1F *numR= (TH1F*)RratioNum[tbin]->Clone(Form("t bin %d num",tbin));
		TH1F *denomR= (TH1F*)RratioDenom[tbin]->Clone(Form("t bin %d denom",tbin));
		
	
		numR->Multiply(f3,1.);
		
		canRatio->cd(tbin+1);numR->Draw();
		canRatio->cd(tbin+1+4);RratioNum[tbin]->Draw();
		
		double tbincenter=(AcceptancetM2forbinning[0]->GetXaxis()->GetBinCenter(1+tbin));
		
		double RC = numR->Integral();
		double R0 = denomR->Integral();

		cout<<"RC R0 "<<RC<<" "<<R0<<endl;
		double R = RC/R0;
		
		
		double errorY=smear(*RratioNum[tbin],tbin);
		double bincenter=(tmeanRratio->GetBinContent(tbin+1))/(tmeanRratioN->GetBinContent(tbin+1));
		Rratio->SetPoint(tbin,bincenter,R);
		Rratio->SetPointError(tbin,(bincenter-t1[tbin]),(t1[tbin+1]-bincenter),errorY,errorY);
		//Rratio->SetPointError(tbin,0.0,fiterror);
	}
	canRatio->cd(9);
	Rratio->SetFillColor(kRed);
   	Rratio->SetFillStyle(3005);
   	Rratio->SetMarkerSize(3);
   	Rratio->SetMarkerColor(kRed);
   	Rratio->SetLineColor(kRed);
	Rratio->SetTitle(";-t;R'");
	Rratio->Draw("A P");
	canRatio->SaveAs("canRatio.pdf");
	Rratio->SaveAs("RRatiovst.root");
	
	//R ratio High Mass
	TGraphAsymmErrors *RratioHighMass = new TGraphAsymmErrors(4);
	TCanvas *canRatioHighMass  = new TCanvas("canRatioHighMass","canRatioHighMass",8000,4000);
	canRatioHighMass->Divide(4,3);
	for(int tbin=0;tbin<4;tbin++){
	
		TH1F *numRHighMass= (TH1F*)RratioNumHighMass[tbin]->Clone(Form("t bin %d numHighMass",tbin));
		TH1F *denomRHighMass= (TH1F*)RratioDenomHighMass[tbin]->Clone(Form("t bin %d denomHighMass",tbin));
		
	
		numRHighMass->Multiply(f3,1.);
		
		canRatioHighMass->cd(tbin+1);numRHighMass->Draw();
		canRatioHighMass->cd(tbin+1+4);RratioNumHighMass[tbin]->Draw();
		
		double tbincenterHighMass=(AcceptancetM2forbinning[0]->GetXaxis()->GetBinCenter(1+tbin));
		
		double RCHighMass = numRHighMass->Integral();
		double R0HighMass = denomRHighMass->Integral();

		cout<<"RC R0 "<<RCHighMass<<" "<<R0HighMass<<endl;
		double R = RCHighMass/R0HighMass;
		
		
		double errorY=smear(*RratioNumHighMass[tbin],tbin);
		double bincenter=(tmeanRratioHighMass->GetBinContent(tbin+1))/(tmeanRratioNHighMass->GetBinContent(tbin+1));
		RratioHighMass->SetPoint(tbin,bincenter,R);
		RratioHighMass->SetPointError(tbin,(bincenter-t1b[tbin]),(t1b[tbin+1]-bincenter),errorY,errorY);
		//Rratio->SetPointError(tbin,0.0,fiterror);
	}
	canRatioHighMass->cd(9);
	RratioHighMass->SetFillColor(kRed);
   	RratioHighMass->SetFillStyle(3005);
   	RratioHighMass->SetMarkerSize(3);
   	RratioHighMass->SetMarkerColor(kRed);
   	RratioHighMass->SetLineColor(kRed);
	RratioHighMass->SetTitle(";-t;R'");
	RratioHighMass->Draw("A P");
	canRatioHighMass->SaveAs("canRatioHighMass.pdf");
	RratioHighMass->SaveAs("RRatiovstHighMass.root");
	
	//R ratio Rafo
	TGraphAsymmErrors *RratioRafo = new TGraphAsymmErrors(4);
	TCanvas *canRatioRafo  = new TCanvas("canRatioRafo","canRatio",8000,4000);
	canRatioRafo->Divide(4,3);
	if(Rafo){
	for(int tbin=0;tbin<8;tbin++){
	
		TH1F *numRRafo= (TH1F*)RratioNumRafo[tbin]->Clone(Form("t bin %d num Rafo",tbin));
		TH1F *denomRRafo= (TH1F*)RratioDenomRafo[tbin]->Clone(Form("t bin %d denom Rafo",tbin));
		
	
		numRRafo->Multiply(f3,1.);
		
		canRatioRafo->cd(tbin+1);numRRafo->Draw();
		canRatioRafo->cd(tbin+1+4);RratioNumRafo[tbin]->Draw();
		
		double tbincenterRafo=(AcceptancetM2forbinning[0]->GetXaxis()->GetBinCenter(1+tbin));
		
		double RCRafo = numRRafo->Integral();
		double R0Rafo = denomRRafo->Integral();

		cout<<"RC R0 "<<RCRafo<<" "<<R0Rafo<<endl;
		double RRafo = RCRafo/R0Rafo;
		
		
		double errorYRafo=smear(*RratioNumRafo[tbin],tbin);
		double bincenter=(tmeanRratioRafo->GetBinContent(tbin+1))/(tmeanRratioNRafo->GetBinContent(tbin+1));
		RratioRafo->SetPoint(tbin,bincenter,RRafo);
		RratioRafo->SetPointError(tbin,(bincenter-t1Rafo[tbin]),(t1Rafo[tbin+1]-bincenter),errorYRafo,errorYRafo);
		//Rratio->SetPointError(tbin,0.0,fiterror);
	}
	}
	canRatioRafo->cd(9);
	RratioRafo->SetFillColor(kRed);
   	RratioRafo->SetFillStyle(3005);
   	RratioRafo->SetMarkerSize(3);
   	RratioRafo->SetMarkerColor(kRed);
   	RratioRafo->SetLineColor(kRed);
	RratioRafo->SetTitle(";-t;R'");
	RratioRafo->GetHistogram()->SetMaximum(1.);   // along          
   	RratioRafo->GetHistogram()->SetMinimum(-1.);
	RratioRafo->Draw("A P");
	Double_t x[4]  = {0.16,0.28,0.44,0.72};
   	Double_t y[4]  = {0.339593,0.534089,0.502156,0.331439};
   	Double_t ex[4] = {0.,0.,0.,0.};
   	Double_t ey[4] = {0.101213,0.112986,0.0969139,0.0962636};
   	TGraphErrors *RratioRafoThesis = new TGraphErrors(4,x,y,ex,ey);
   	RratioRafoThesis->SetMarkerColor(kBlue);
   	RratioRafoThesis->SetLineColor(kBlue);
   	RratioRafoThesis->Draw("P");
	canRatioRafo->SaveAs("canRatioRafo.pdf");
	RratioRafo->SaveAs("RRatiovstRafo.root");
	
	//R ratio Xi binning
	TGraphAsymmErrors *RratioXi = new TGraphAsymmErrors(3);
	TCanvas *canRatioXi  = new TCanvas("canRatioXi","canRatioXi",8000,4000);
	canRatioXi->Divide(4,3);
	for(int Xibin=0;Xibin<3;Xibin++){
	
		TH1F *numRXi= (TH1F*)RratioNumXi[Xibin]->Clone(Form("Xi bin %d num",Xibin));
		TH1F *denomRXi= (TH1F*)RratioDenomXi[Xibin]->Clone(Form("Xi bin %d denom",Xibin));
		
	
		numRXi->Multiply(f3,1.);
		
		canRatioXi->cd(Xibin+1);numRXi->Draw();
		canRatioXi->cd(Xibin+1+4);RratioNumXi[Xibin]->Draw();
		
		double Xibincenter=(PhysicsXi->GetXaxis()->GetBinCenter(1+Xibin));
		
		double RCXi = numRXi->Integral();
		double R0Xi = denomRXi->Integral();

		cout<<"Xi RC R0 "<<RCXi<<" "<<R0Xi<<endl;
		double RXi = RCXi/R0Xi;
		
		
		double errorY=smear(*RratioNumXi[Xibin],Xibin);
		double bincenter=(XimeanRratio->GetBinContent(Xibin+1))/(XimeanRratioN->GetBinContent(Xibin+1));
		RratioXi->SetPoint(Xibin,bincenter,RXi);
		RratioXi->SetPointError(Xibin,(bincenter-Xi3[Xibin]),(Xi3[Xibin+1]-bincenter),errorY,errorY);
		cout<<"in xi binning "<<bincenter<<" "<<Xi3[Xibin]<<" "<<Xi3[Xibin+1]<<endl;
		//Rratio->SetPointError(tbin,0.0,fiterror);
	}
	canRatioXi->cd(9);
	RratioXi->SetFillColor(kRed);
   	RratioXi->SetFillStyle(3005);
   	RratioXi->SetMarkerSize(3);
   	RratioXi->SetMarkerColor(kRed);
   	RratioXi->SetLineColor(kRed);
	RratioXi->SetTitle(";#xi;R'");
	RratioXi->Draw("A P");
	canRatioXi->SaveAs("canRatioXi.pdf");
	RratioXi->SaveAs("RratioXi.root");
	
	
	
	
	
	
	
	TCanvas *canpola  = new TCanvas("canpola","canpola",5000,1000);
	canpola->Divide(4,1);
	canpola->cd(1);pospola->Draw();
	canpola->cd(2);negpola->Draw();
	TF1  *f1 = new TF1("f1","[0]*sin(x*[1]/[2])",-180,180);
	f1->FixParameter(1,3.14159264);
	f1->FixParameter(2,180);
	canpola->cd(3);numpola->Fit("f1");//Draw();
	PolarizationTransfer->SetMaximum(1.2);
	canpola->cd(4);PolarizationTransfer->Draw();
	canpola->SaveAs("canpola.pdf");
	numpola->SaveAs("pola.root");

	TCanvas *canpolaTranfer  = new TCanvas("canpolaTranfer","canpolaTranfer",3000,2000);
	gStyle->SetOptStat(0);
	 PolarizationTransfer->GetYaxis()->SetRangeUser(0., 2.);
	PolarizationTransfer->SetTitle(";E_{#gamma};Pola. Transfer");
	PolarizationTransfer->Draw();
	canpolaTranfer->SaveAs("polarizationTransfer.pdf");
		gStyle->SetOptStat(111);
		
		
	Assym->SaveAs("Assym.root");
	cout<<"real number od equivalent tcs "<<(nEventTCS*nEventTCS)/denom<<endl;

	
	
	cout<<" Mean kinematics "<<endl;
	cout<<"AFBt0_Eg"<<endl;
	TCanvas *canAFBt0_Eg  = new TCanvas("canAFBt0_Eg","",500,500);      
	AFBt0_Q2->Draw();
	canAFBt0_Eg->SaveAs("canAFBt0_Eg.pdf");
	cout<<(AFBt0_Eg->GetMean(1))<<"  "<<(AFBt0_Eg->GetRMS(1))<<endl;
	cout<<AFBt0_Q2->GetMean(1)<<"  "<<(AFBt0_Q2->GetRMS(1))<<endl;
	cout<<	endl;
	cout<<"AFBEg0"	<<endl;							
	cout<<AFBEg0_t->GetMean(1)<<"  "<<(AFBEg0_t->GetRMS(1))<<endl;
	cout<<AFBEg0_Q2->GetMean(1)<<"  "<<(AFBEg0_Q2->GetRMS(1))<<endl;
	cout<<	endl;
	cout<<"AFBQ20"	<<endl;				
	cout<<AFBQ20_Eg->GetMean(1)<<"  "<<(AFBQ20_Eg->GetRMS(1))<<endl;
	cout<<AFBQ20_t->GetMean(1)<<"  "<<(AFBQ20_t->GetRMS(1))<<endl;
	cout<<	endl;
	cout<<"AFBXi0"		<<endl;			
	cout<<AFBXi0_Q2->GetMean(1)<<"  "<<(AFBXi0_Q2->GetRMS(1))<<endl;
	cout<<AFBXi0_Eg->GetMean(1)<<"  "<<(AFBXi0_Eg->GetRMS(1))<<endl;
	cout<<AFBXi0_t->GetMean(1)<<"  "<<(AFBXi0_t->GetRMS(1))<<endl;
	cout<<	endl;
	cout<<"AFBt1"		<<endl;				
	cout<<AFBt1_Eg->GetMean(1)<<"  "<<(AFBt1_Eg->GetRMS(1))<<endl;
	cout<<AFBt1_Q2->GetMean(1)<<"  "<<(AFBt1_Q2->GetRMS(1))<<endl;
	cout<<endl;
	cout<<"AFBt2"<<endl;
	cout<<AFBt2_Eg->GetMean(1)<<"  "<<(AFBt2_Eg->GetRMS(1))<<endl;
	cout<<AFBt2_Q2->GetMean(1)<<"  "<<(AFBt2_Q2->GetRMS(1))<<endl;
	cout<<	endl;
	cout<<"RratioHighMass"	<<endl;				
	cout<<RratioHighMass_Eg->GetMean(1)<<"  "<<(RratioHighMass_Eg->GetRMS(1))<<endl;
	cout<<RratioHighMass_Q2->GetMean(1)<<"  "<<(RratioHighMass_Q2->GetRMS(1))<<endl;
	cout<<	endl;
	cout<<"Rratio_Eg"<<endl;				
	cout<<Rratio_Eg->GetMean(1)<<"  "<<(Rratio_Eg->GetRMS(1))<<endl;
	cout<<Rratio_Q2->GetMean(1)<<"  "<<(Rratio_Q2->GetRMS(1))<<endl;
	cout<<Rratio_t->GetMean(1)<<"  "<<(Rratio_t->GetRMS(1))<<endl;
	cout<<endl;
	
	cout<<"EgmeanBSA method 1"<<endl;
	cout<<BSA_Eg->GetMean(1)<<(BSA_Eg->GetRMS(1))<<endl;
	cout<<"MmeanBSA method 1"<<endl;
	cout<<BSA_Q2->GetMean(1)<<(BSA_Q2->GetRMS(1))<<endl;
	cout<<"tmeanBSA method 1"<<endl;
	cout<<BSA_t->GetMean(1)<<(BSA_t->GetRMS(1))<<endl;
	
	cout<<"Q2meanBSA"<<endl;
	cout<<Q2meanBSA->GetMean(1)<<endl;
	cout<<endl;
	cout<<"tmeanBSA"<<endl;
	cout<<	(tmeanBSA->Integral()) /(tmeanBSAN->Integral())<<endl;
	cout<<"XimeanBSA"				<<endl;
	cout<<(XimeanBSA->Integral()) /(XimeanBSAN->Integral())<<endl;
	cout<<"MmeanBSA"					<<endl;
	cout<<(MmeanBSA->Integral())/(MmeanBSAN->Integral())<<endl;
	cout<<"BSA_theta"<<endl;
	cout<<BSA_theta->GetMean(1)<<endl;
	
	
	
	TCanvas *canyc  = new TCanvas("canyc","yc",4000,1000);
	canyc->Divide(4,2);
	canyc->cd(1);Yc[0]->Draw("hist");
	canyc->cd(2);Yc[1]->Draw("hist");
	canyc->cd(3);Yc[2]->Draw("hist");
	canyc->cd(4);Yc[3]->Draw("hist");
	canyc->cd(5);EYc[0]->Draw("hist");
	canyc->cd(6);EYc[1]->Draw("hist");
	canyc->cd(7);EYc[2]->Draw("hist");
	canyc->cd(8);EYc[3]->Draw("hist");
	canyc->SaveAs("yc.pdf");


	TCanvas *canBinning  = new TCanvas("canBinning","canBinning",4000,1000);
	canBinning->Divide(3,1);
	canBinning->cd(1);EgVST->Draw("col");
	canBinning->cd(2);MVST->Draw("col");
	canBinning->cd(3);MVSEg->Draw("col");
	canBinning->SaveAs("canBinning.pdf");
	EgVST->SaveAs("EgVST.root");
	MVST->SaveAs("MVST.root");
	MVSEg->SaveAs("MVSEg.root");
	
	/*TCanvas *canA  = new TCanvas("can1","tcs1",3000,3000);
	canA->Divide(4,3);
	Acc[0]->SetTitle("t bin 1");
	Acc[1]->SetTitle("t bin 2");
	Acc[2]->SetTitle("t bin 3");
	Acc[3]->SetTitle("t bin 4");
	canA->cd(1);Acc[0]->Draw("box");PhiVSThetaCM->Draw("same");MaxAcc0->Draw("same");MinAcc0->Draw("same");
	canA->cd(2);Acc[1]->Draw("box");PhiVSThetaCM1->Draw("same");MaxAcc1->Draw("same");MinAcc1->Draw("same");
	canA->cd(3);Acc[2]->Draw("box");PhiVSThetaCM2->Draw("same");MaxAcc2->Draw("same");MinAcc2->Draw("same");
	canA->cd(4);Acc[3]->Draw("box");PhiVSThetaCM3->Draw("same");MaxAcc3->Draw("same");MinAcc3->Draw("same");
	canA->cd(5);PhiVSThetaCM->Draw("colz");MaxAcc0->Draw("same");MinAcc0->Draw("same");
	canA->cd(6);PhiVSThetaCM1->Draw("colz");MaxAcc1->Draw("same");MinAcc1->Draw("same");
	canA->cd(7);PhiVSThetaCM2->Draw("colz");MaxAcc2->Draw("same");MinAcc2->Draw("same");
	canA->cd(8);PhiVSThetaCM3->Draw("colz");MaxAcc3->Draw("same");MinAcc3->Draw("same");
	canA->cd(9);Assym->Draw("AP");
	canA->cd(10);FinalEventst->Draw();
	canA->cd(11);FinalEventsEg->Draw();
	canA->cd(12);FinalEventsM->Draw();
	canA->SaveAs("tcs1Ana.pdf");*/
	
	/*gStyle->SetFrameLineWidth(4); 
	gStyle->SetLineWidth(5);
	gStyle->SetHistLineWidth(5);*/
	TCanvas *canAp  = new TCanvas("can1p","tcs1",6000,4000);
	canAp->Divide(3,2);
	FinalEventst->SetTitle("");
	FinalEventsEg->SetTitle("");
	FinalEventsM->SetTitle("");
	canAp->cd(1);FinalEventst->Draw();
	canAp->cd(2);FinalEventsEg->Draw();
	canAp->cd(3);FinalEventsM->Draw();
	canAp->cd(4);EgVSPeletron->Draw("col");
	canAp->cd(5);tVSPeletron->Draw("col");
	canAp->SaveAs("tcs3Ana.pdf");
	FinalEventst->SaveAs("FinalEventst.root");
	FinalEventsEg->SaveAs("FinalEventsEg.root");
	FinalEventsM->SaveAs("FinalEventsM.root");
	FinalEventsPhi->SaveAs("FinalEventsPhiProt.root");
	FinalEventsTheta->SaveAs("FinalEventsThetaProt.root");
	FinalEventsP->SaveAs("FinalEventsPProt.root");
	
	TCanvas *canAr  = new TCanvas("can1r","tcs1",6000,4000);
	canAr->Divide(2,2);
	PhiVSThetaCM->SetTitle("t bin 1;#phi;#theta");
	PhiVSThetaCM1->SetTitle("t bin 2;#phi;#theta");
	PhiVSThetaCM2->SetTitle("t bin 3;#phi;#theta");
	PhiVSThetaCM3->SetTitle("t bin 4;#phi;#theta");
	canAr->cd(1);PhiVSThetaCM->Draw("colz");//MaxAcc0->Draw("same");MinAcc0->Draw("same");
	canAr->cd(2);PhiVSThetaCM1->Draw("colz");//MaxAcc1->Draw("same");MinAcc1->Draw("same");
	canAr->cd(3);PhiVSThetaCM2->Draw("colz");//MaxAcc2->Draw("same");MinAcc2->Draw("same");
	canAr->cd(4);PhiVSThetaCM3->Draw("colz");//MaxAcc3->Draw("same");MinAcc3->Draw("same");
	canAr->SaveAs("tcs2Ana.pdf");
	
	TCanvas *canQ2  = new TCanvas("canQ2","canQ2",6000,4000);
	gPad->SetLogy();
	EhistQ2->Draw();
	EhistQ21->SetLineColor(kRed);
	//EhistQ21->SetHistColor(kRed);
	EhistQ21->Draw("same");
	canQ2->SaveAs("canQ2.pdf");



/*
	TCanvas *can  = new TCanvas("can","tcs",5000,23000);
	can->Divide(5,23);
	can->cd(1);histMMass->Draw();
	can->cd(2);histQ2->Draw();
	can->cd(3);gPad->SetLogz();histPMiss->Draw("col");
	can->cd(4);histMMPt->Draw("col");
	can->cd(5);histQP2->Draw();
	can->cd(6);histT->Draw();
	can->cd(7);histQ2t->Draw();
	can->cd(8);PhiVSPElectron->Draw("col");
	can->cd(9);PhiVSPPositron->Draw("col");
	can->cd(10);ThetaVSPElectron->Draw("col");
	can->cd(11);ThetaVSPPositron->Draw("col");
	can->cd(12);PhiVSThetaElectron->Draw("col");
	can->cd(13);PhiVSThetaPositron->Draw("col");
	can->cd(14);PhiVSThetaProton->Draw("col");
	can->cd(15);ThetaVSPProton->Draw("col");
	can->cd(16);EnergyConservation->Draw();
	can->cd(17);MMQ2->Draw("col");
	can->cd(18);EhistMMass->Draw();
	can->cd(19);EhistQ2->Draw();
	can->cd(20);gPad->SetLogz();EhistPMiss->Draw("col");
	can->cd(21);EhistMMPt->Draw("col");
	can->cd(22);EhistQP2->Draw();
	can->cd(23);EEnergyConservation->Draw();
	can->cd(24);EMMQ2->Draw("col");
	can->cd(25);EPhiVSPElectron->Draw("col");
	can->cd(26);EPhiVSPPositron->Draw("col");
	can->cd(27);EThetaVSPElectron->Draw("col");
	can->cd(28);EThetaVSPPositron->Draw("col");
	can->cd(29);EPhiVSThetaElectron->Draw("col");
	can->cd(30);EPhiVSThetaPositron->Draw("col");
	can->cd(31);EPhiVSThetaProton->Draw("col");
	can->cd(32);EThetaVSPProton->Draw("col");
	can->cd(33);EhistT->Draw();
	can->cd(34);SFelectron->Draw("col");
	can->cd(35);SFpositron->Draw("col");
	can->cd(36);CheElectron->Draw("col");
	can->cd(37);ChePositron->Draw("col");
	can->cd(38);ECelectron->Draw("col");
	can->cd(39);ECpositron->Draw("col");
	can->cd(40);BetaProton->Draw("col");pBetaProton->Draw("same");
	can->cd(41);EECelectron->Draw("col");
	can->cd(42);EECpositron->Draw("col");	
	can->cd(43);SFUelectron->Draw("col");	
	can->cd(44);SFVelectron->Draw("col");	
	can->cd(45);SFWelectron->Draw("col");
	can->cd(46);SFUpositron->Draw("col");
	can->cd(47);SFVpositron->Draw("col");
	can->cd(48);SFWpositron->Draw("col");
	can->cd(56);PhiVSThetaCM->Draw("colz");
	can->cd(57);PhiVSThetaCM1->Draw("colz");
	can->cd(58);PhiVSThetaCM2->Draw("colz");
	can->cd(59);PhiVSThetaCM3->Draw("colz");
	can->cd(64);vertexProtTheta->Draw("colz");
	can->cd(65);vertexTimePair->Draw();
	can->cd(66);vertexTimeP->Draw();
	can->cd(67);vertexPair->Draw();
	can->cd(68);vertexElecP->Draw();
	can->cd(69);vertexPosiP->Draw();
	can->cd(70);vertexElec->Draw();
	can->cd(71);vertexPosi->Draw();
	can->cd(72);vertexProt->Draw();
	can->cd(73);vertexTimePP->Draw("col");	
	can->cd(74);corrSFelectron->Draw("col");
	can->cd(75);corrSFpositron->Draw("col");			
	can->cd(76);histM->Draw();
	can->cd(77);EhistM->Draw();EhistMprim->Draw("same");
	can->cd(78);histME->Draw("col");
	can->cd(79);histMPElec->Draw("col");
	can->cd(80);histMPPosi->Draw("col");
	can->cd(81);EhistME->Draw("col");
	can->cd(82);EhistMPElec->Draw("col");
	can->cd(83);EhistMPPosi->Draw("col");			
	can->cd(84);EhistMPProton->Draw("col");			
	can->cd(85);histMME->Draw("col");
	can->cd(86);EhistMME->Draw("col");
	can->cd(87);EhistMXProton->Draw("col");
	can->cd(88);EBetaProton->Draw("col");
	can->cd(89);EvertexTimePair->Draw();
	can->cd(90);EvertexTimeP->Draw();
	can->cd(91);EvertexPair->Draw();
	can->cd(92);EvertexElecP->Draw();
	can->cd(93);EvertexPosiP->Draw();
	can->cd(94);EvertexTimePPCD->Draw("col");
	can->cd(95);EvertexTimePPFD->Draw("col");
	can->cd(96);EhistMPhiProton->Draw("col");
	can->cd(97);histMPhiProton->Draw("col");
	can->cd(98);histMPhiElectron->Draw("col");
	can->cd(99);histMPhiPositron->Draw("col");
	can->cd(100);histMThetaElectron->Draw("col");
	can->cd(101);histMThetaPositron->Draw("col");
	can->cd(102);ECheElectron->Draw("col");
	can->cd(103);EChePositron->Draw("col");
	can->cd(104);PhiVSPProton->Draw("col");
	can->cd(105);EPhiVSPProton->Draw("col");
	can->cd(106);Acc[0]->Draw("box");PhiVSThetaCM->Draw("same text");
	can->cd(107);Acc[1]->Draw("box");PhiVSThetaCM1->Draw("same text");
	can->cd(108);Acc[2]->Draw("box");PhiVSThetaCM2->Draw("same text");
	can->cd(109);Acc[3]->Draw("box");PhiVSThetaCM3->Draw("same text");
	can->cd(110);Assym->Draw("AP");
	can->cd(111);EhistMSectorDiff->Draw("col");
	can->cd(112);EhistEssai->Draw();
	can->cd(113);histMThetaProton->Draw("col");
	can->cd(114);gPad->SetLogy();Pproton->Draw();
	can->SaveAs(type+"tcs1.png");
*/

	gStyle->SetOptStat("e");
	TCanvas *canCheckT  = new TCanvas("canCheckT","canCheckT",4000,2000);
	canCheckT->Divide(3,2);
	canCheckT->cd(1);EnergyConservation->Draw();
	EnergyConservation->SetTitle("All events; #Delta E");
	canCheckT->cd(2);EEnergyConservation->Draw();
	EEnergyConservation->SetTitle("After exclu.; #Delta E");
	canCheckT->cd(3);EEnergyConservationEvents->Draw();
	EEnergyConservationEvents->SetTitle("Final events; #Delta E");
	canCheckT->cd(4);tCheckBefore->Draw();
	canCheckT->cd(5);tCheck->Draw();
	canCheckT->cd(6);tCheckEvents->Draw();
	canCheckT->SaveAs("canCheckT.pdf");
	canCheckT->SaveAs("canCheckT.root");
	
	TCanvas *canChi2 = new TCanvas("canChi2","canChi2",4000,2000);
	canChi2->Divide(4,2);
	canChi2->cd(1);Chi2ElectronAvant->Draw();
	canChi2->cd(2);Chi2PositronAvant->Draw();
	canChi2->cd(3);Chi2ProtonAvant->Draw();
	canChi2->cd(4);Chi2ProtonFDAvant->Draw();
	canChi2->cd(5);Chi2Electron->Draw();
	canChi2->cd(6);Chi2Positron->Draw();
	canChi2->cd(7);Chi2Proton->Draw();
	canChi2->cd(8);Chi2ProtonFD->Draw();
	canChi2->SaveAs("canChi2.pdf");
	canChi2->SaveAs("canChi2.root");
	Chi2ProtonAvant->SaveAs("Chi2ProtonAvant.root");
	Chi2ProtonFDAvant->SaveAs("Chi2ProtonFDAvant.root");
	
	TCanvas *canVT  = new TCanvas();
	canVT->Divide(2,2);
	canVT->cd(1);EvertexTimePair->Draw();
	canVT->cd(2);gPad->SetLogz();vertexTimePairMass->Draw("colz");
	canVT->cd(3);gPad->SetLogz();TimeCheDiffMass->Draw("colz");
	canVT->SaveAs("Canvt.pdf");

	TCanvas *can0  = new TCanvas("can0","tcs0",3000,3000);
	can0->Divide(2,2);
	can0->cd(1);HPosi->Draw();
	can0->cd(2);HNega->Draw();
	can0->cd(3);BAssym->Draw("AP");
	can0->cd(4);
	can0->SaveAs(type+"tcsBeamAssymetry.pdf");	



	gStyle->SetOptStat(0);
	TCanvas *canB  = new TCanvas("can1","tcs1",4000,2000);
	canB->Divide(3,2);
	ThetaVSPElectron->SetTitle("Electron");
	ThetaVSPPositron->SetTitle("Positron");
	ThetaVSPProton->SetTitle("Proton");
	EThetaVSPElectron->SetTitle("Electron");
	EThetaVSPPositron->SetTitle("Positron");
	EThetaVSPProton->SetTitle("Proton");
	canB->cd(1);ThetaVSPElectron->Draw("col");
	canB->cd(2);ThetaVSPPositron->Draw("col");
	canB->cd(3);ThetaVSPProton->Draw("col");
	canB->cd(4);EThetaVSPElectron->Draw("col");
	canB->cd(5);EThetaVSPPositron->Draw("col");
	canB->cd(6);EThetaVSPProton->Draw("col");
	canB->SaveAs(type+"tcs1AnaKine.pdf");
	
	TCanvas *canPositron  = new TCanvas("canPositron","canPositron",3000,2000);
	ThetaVSPPositron->Draw("col");
	canPositron->SaveAs("PositronKine.pdf");
	canPositron->SaveAs("PositronKine.root");

	TCanvas *canD  = new TCanvas("can2","tcs2",6000,1600);
	canD->Divide(3,1);
	histPMiss->SetTitle("Pt/P scattered electron");
	canD->cd(1);gPad->SetLogz();histPMiss->Draw("colz");
	histMMPt->SetTitle("Pt/P vs Mass^{2} scattered electron ");
	histMMPt->GetXaxis()->SetTitle("Mass^{2} e (GeV^{2})");
	histMMPt->GetYaxis()->SetTitle("Pt/P");
	histMMPt->GetYaxis()->SetTitleOffset(1.1);
	canD->cd(2);gPad->SetLogz();histMMPt->Draw("colz");
	MMQ2->SetTitle("Mass^{2} photon vs Mass^{2} scattered electron ");
	MMQ2->GetXaxis()->SetTitle("Mass^{2} e (GeV^{2})");
	MMQ2->GetYaxis()->SetTitle("Mass^{2} #gamma (GeV^{2})");
	canD->cd(3);gPad->SetLogz();MMQ2->Draw("colz");
	canD->SaveAs(type+"ExclusivityCuts.pdf");
	
	TCanvas *canDbis  = new TCanvas("can2","tcs2",12000,4000);
	canDbis->Divide(2,1);
	histMMPt->SetTitle("");
	
	histMMPt->SetStats(kFALSE);
	//gPad->SetMargin(0.1,0.1,0.1,0.1);
	gPad->SetRightMargin(0.5);
	canDbis->cd(1);gPad->SetLogz();histMMPt->Draw("col");
	
	gPad->Update();
	canDbis->SaveAs("ExclusivityCutsJolie.pdf");

	TCanvas *canB1  = new TCanvas("canB1","tcs1",4000,2000);
	canB1->Divide(3,2);
	ThetaVSPElectron->SetTitle("Electron");
	ThetaVSPPositron->SetTitle("Positron");
	ThetaVSPProton->SetTitle("Proton");
	TCSEThetaVSPElectron->SetTitle("Electron");
	TCSEThetaVSPPositron->SetTitle("Positron");
	TCSEThetaVSPProton->SetTitle("Proton");
	canB1->cd(1);ThetaVSPElectron->Draw("col");
	canB1->cd(2);gPad->SetLogz();ThetaVSPPositron->Draw("col");
	canB1->cd(3);ThetaVSPProton->Draw("col");
	canB1->cd(4);TCSEThetaVSPElectron->Draw("col");
	canB1->cd(5);TCSEThetaVSPPositron->Draw("col");
	canB1->cd(6);TCSEThetaVSPProton->Draw("col");
	canB1->SaveAs(type+"tcs1AnaKineTCS.pdf");


	TCanvas *canMass  = new TCanvas("canMass","Mass",6000,4000);
	canMass->Divide(1,1);
	canMass->cd(1);
	EhistM->SetLineColor(kRed);
	gPad->SetLogy();
	EhistM->SetStats(kFALSE);
	EhistM->SetTitle(";M (GeV); events");
	EhistM->Draw();
	
	//EhistMAfterFid->Draw("same");
	canMass->SaveAs("Mass.pdf");
	EhistM->SaveAs("Mass.root");

	TCanvas *canMassRatio  = new TCanvas("canMassRatio","MassRatio",4000,4000);
	canMassRatio->Divide(1,1);
	canMassRatio->cd(1);
	EhistM->Divide(EhistMAfterFid);
	EhistM->Draw();
	//EhistM->Draw("same");
	canMassRatio->SaveAs("MassRatio.pdf");

	TCanvas *canOther  = new TCanvas("canOther","Other",8000,4000);
	canOther->Divide(2,1);
	canOther->cd(1);
	VertexOtherPart->Draw();
	canOther->cd(2);
	Chi2OtherPart->Draw();
	canOther->SaveAs("Other.pdf");

	gStyle->SetOptFit(1);
	TCanvas *canCheckMass  = new TCanvas("canCheckMass","CheckMass",12000,12000);
	canCheckMass->Divide(3,3);
	canCheckMass->cd(1);EhistMPElec->Draw("col");
	canCheckMass->cd(2);EhistMPPosi->Draw("col");			
	canCheckMass->cd(3);EhistMPProton->Draw("col");	
	canCheckMass->cd(4);histMThetaElectron->Draw("col");
	canCheckMass->cd(5);histMThetaPositron->Draw("col");
	canCheckMass->cd(6);histMThetaProton->Draw("col");
	canCheckMass->cd(7);histCheCheElectron->Draw("col");
	TF1 *fvt = new TF1("fvt","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+x*[4]",0.1,0.2);//+[3]+x*[4]+x*x*[5]
	fvt->SetParameter(0,40.);
	fvt->SetParLimits(0,20,50);
	fvt->SetParameter(1,0.15);
	fvt->SetParLimits(1,0.1,0.2);
	fvt->SetParameter(2,0.1);
	fvt->SetParameter(3,0.0);
	fvt->SetParLimits(3,0,10);
	canCheckMass->cd(8);FuckingPi0->Fit("fvt","","",0.1,0.17);
	canCheckMass->SaveAs("CheckMass.pdf");


	TCanvas *canCheckVertex  = new TCanvas("canCheckVertex","ChecVertex",12000,8000);
	canCheckVertex->Divide(3,2);
	canCheckVertex->cd(1);histMVElectron->Draw("col");
	canCheckVertex->cd(2);histMVPositron->Draw("col");
	canCheckVertex->cd(3);histMV->Draw("col");
	canCheckVertex->cd(4);histMCheElectron->Draw("col");
	canCheckVertex->cd(5);histMChePositron->Draw("col");
	canCheckVertex->cd(6);EvertexPair->Draw("col");//->Fit("gaus");
	canCheckVertex->SaveAs("CheckVertex.pdf");

	TCanvas *canCheckVertex1  = new TCanvas("canCheckVertex1","Vertex",8000,16000);
	canCheckVertex1->Divide(2,4);
	canCheckVertex1->cd(1);vertexProtThetaFD->Draw("col");
	canCheckVertex1->cd(2);vertexProtThetaCD->Draw("col");
	canCheckVertex1->cd(3);vertexProtThetaFD->ProjectionY()->Draw();
	canCheckVertex1->cd(4);vertexProtThetaCD->ProjectionY()->Draw();
	canCheckVertex1->cd(6);vertexElecP->Draw();
	canCheckVertex1->cd(7);vertexPosiP->Draw();
	canCheckVertex1->cd(8);vertexElec->Draw();
	canCheckVertex1->SaveAs("CheckVertex1.pdf");

	TCanvas *canPCAL  = new TCanvas("A","A",3000,3000);
	canPCAL->Divide(1,1);
	canPCAL->cd(1);BeforeCuts->Draw("col");
	canPCAL->SaveAs("PCAL.pdf");

	TCanvas *canPCALPositron  = new TCanvas("Aff","A",3000,3000);
	canPCALPositron->Divide(1,1);
	canPCALPositron->cd(1);BeforeCutsPositron->Draw("col");
	canPCALPositron->SaveAs("PCALPositron.pdf");

	TCanvas *canCHEPositron  = new TCanvas("Ahhhff","A",6000,3000);
	canCHEPositron->Divide(2,1);
	canCHEPositron->cd(1);CheCooElec->Draw("col");
	canCHEPositron->cd(2);CheCooPosi->Draw("col");
	canCHEPositron->SaveAs("CHEposition.pdf");

	TCanvas *canCHEdiff  = new TCanvas("Ahhghff","A",9000,3000);
	canCHEdiff->Divide(3,1);
	canCHEdiff->cd(1);CheCooDiff->Draw("col");
	canCHEdiff->cd(2);CheCooDiffInf->Draw("col");
	canCHEdiff->cd(3);CheCooDiffSup->Draw("col");
	canCHEdiff->SaveAs("CHEpositiondiff.pdf");


	TCanvas *Cheelec = new TCanvas("Cheelec","Cheelec",2500,1500);
	CheElectron->SetTitle(";N_{PHE};counts;");
	CheElectron->SetStats(0);
	CheElectron->SetLabelSize(.04, "xyz");
	CheElectron->SetTitleSize(.04, "xyz");
	ECheElectron->Draw("hist");
	ECheElectron->SaveAs("Canche.root");
	Cheelec->SaveAs("Canche.pdf");

	TCanvas *canTestPio  = new TCanvas("TestPio","TestPio",15000,12000);
	canTestPio->Divide(4,5);
	canTestPio->cd(1);VertexPio->Draw("col");
	canTestPio->cd(2);MomentumPairPio->Draw("col");
	canTestPio->cd(3);ThetaPairPio->Draw("col");
	canTestPio->cd(4);histCheCheElectron->Draw("col");

	canTestPio->cd(5);VertexAutre->Draw("col");
	canTestPio->cd(6);MomentumPairAutre->Draw("col");
	canTestPio->cd(7);ThetaPairAutre->Draw("col");
	canTestPio->cd(8);histCheCheElectronSup->Draw("col");

	canTestPio->cd(9);VertexAutreInf->Draw("col");
	canTestPio->cd(10);MomentumPairAutreInf->Draw("col");
	canTestPio->cd(11);ThetaPairAutreInf->Draw("col");
	canTestPio->cd(12);histCheCheElectronInf->Draw("col");

	canTestPio->cd(13);histAngleMass->Draw("col");
	canTestPio->cd(14);TimeChePioInf->Draw("col");
	canTestPio->cd(15);TimeChePio->Draw("colz");
	canTestPio->cd(16);TimeChePioSup->Draw("colz");
	canTestPio->cd(17);gPad->SetLogy();TimeCheDiff->Draw("colz");
	canTestPio->cd(18);gPad->SetLogz();TimeCheDiffMass->Draw("colz");
	canTestPio->cd(19);SFMassElec->Draw("colz");
	canTestPio->cd(20);SFMassPosi->Draw("colz");
	canTestPio->SaveAs("TestPio.pdf");


	TCanvas *canSFcuts  = new TCanvas("SFcuts ","A",9000,12000);
	canSFcuts->Divide(2,3);
	canSFcuts->cd(1);SFPCALECALElecavant->Draw("col");
	canSFcuts->cd(2);SFPCALECALPosiavant->Draw("col");
	canSFcuts->cd(3);SFPCALECALElecapres->Draw("col");
	canSFcuts->cd(4);SFPCALECALPosiapres->Draw("col");
	canSFcuts->cd(5);SFPpcalzeroPosi->Draw("col");
	canSFcuts->cd(6);SFPpcalzeroElec->Draw("col");

	canSFcuts->SaveAs("SFcuts.pdf");

	TCanvas *canChePhi  = new TCanvas("ChePhi ","A",9000,9000);
	canChePhi->Divide(1,1);
	canChePhi->cd(1);ChePhi->Draw("colz");
	canChePhi->SaveAs("ChePhi.pdf");

	TCanvas *canPP = new TCanvas("ChePP ","A",7000,3000);
		canPP->Divide(2,1);
	PPosiPElec->SetTitleSize(.04, "xyz");
	PPosiPElec1->SetTitleSize(.04, "xyz");
	canPP->cd(1);gPad->SetLogz();PPosiPElec->Draw("col");
	canPP->cd(2);gPad->SetLogz();PPosiPElec1->Draw("col");
	canPP->SaveAs("canPP.pdf");
	
	cout<<"la"<<endl;
	TLine *line1 = new TLine(0.05,4.,0.35,4.);
	line1->SetLineColor(kRed);
	
	TLine *line2 = new TLine(0.11,2.,0.11,9.);
	line2->SetLineColor(kGreen);
	TLine *line3 = new TLine(0.2,2.,0.2,9.);
	line3->SetLineColor(kGreen);
	TCanvas *canxi  = new TCanvas("canxi","",500,500);                                  
	xihist->Draw("col");
	line1->Draw("same");	line2->Draw("same");	line3->Draw("same");
	canxi->SaveAs("canxi.pdf");
	
	TCanvas *canxi1  = new TCanvas("canxi1","",500,500);                                  
	xihist1->Draw("col");
	canxi1->SaveAs("canxi1.pdf");
	
		cout<<" ou la"<<endl;
	TLine *line2b = new TLine(0.11,0.15,0.11,0.8);
	line2b->SetLineColor(kGreen);
	TLine *line3b = new TLine(0.2,0.15,0.2,0.8);
	line3b->SetLineColor(kGreen);
	TCanvas *canti  = new TCanvas("canti","",500,500);                                  
	thist->Draw("col");
	line2b->Draw("same");	line3b->Draw("same");
	canti->SaveAs("canti.pdf");
	
	TCanvas *canti2  = new TCanvas("canti2","",500,500);                                  
	thist1->Draw("col");
	line2b->Draw("same");	line3b->Draw("same");
	canti2->SaveAs("canti2.pdf");
	
	TCanvas *canPhiMissingParticle  = new TCanvas("canMLPplot","",1700,1500);      
	gPad->SetLogz();PhiMissingParticle->Draw("col");
	canPhiMissingParticle->SaveAs("canPhiMissingParticle.pdf");
	
	TCanvas *canMLPplot  = new TCanvas("canMLPplot","",500,500);      
	gPad->SetLogy();MLPplot->Draw("");
	canMLPplot->SaveAs("MLPplot.pdf");
	
TCanvas *canti3  = new TCanvas("canti2","",500,500);      
	Q2t->Draw("colz");
	canti3->SaveAs("Q2t.pdf");

TCanvas *canPositronpt  = new TCanvas("canPositronpt","",500,500);      
	Positronpt->Draw();
	canPositronpt->SaveAs("canPositronpt.pdf");

	TCanvas *canCorrecRad  = new TCanvas("canCorrecRad","",750,500); 
	canCorrecRad->Divide(3,2);
	canCorrecRad->cd(1);diffPhidiffThetaRad->Draw();
	canCorrecRad->cd(2);diffPhidiffThetaRad1->Draw();
	canCorrecRad->cd(3);diffPhidiffThetaRad1vsTheta->Draw("colz");
	canCorrecRad->cd(4);diffPhidiffThetaRad1vsP->Draw("colz");
	canCorrecRad->cd(5);diffPhidiffThetaRad1vsPhi->Draw("colz");
	canCorrecRad->cd(6);DiffRadCor->Draw("colz");
	canCorrecRad->SaveAs("canCorrecRad.pdf");
	
	TCanvas *canExclu1D  = new TCanvas("canExclu1D","",1000,500); 
	canExclu1D->Divide(2,1);
	canExclu1D->cd(1);Pt1D->Draw("hist");
	canExclu1D->cd(2);MassScat->Draw("hist");
	canExclu1D->SaveAs("canExclu1D.pdf");
	Pt1D->SaveAs("Pt1D.root");
	MassScat->SaveAs("MassScat.root");
	Pt1D1->SaveAs("Pt1D1.root");
	MassScat1->SaveAs("MassScat1.root");
	
	RunCheck->SaveAs("RunCheck.root");

	outFile->cd();		
	outFile->Write();
	outFile->Close();


	cout<<"nb of file "<<nbf<<endl;
	cout<<"nb of events "<<nbrecEvent<<endl;
	cout<<"real number od equivalent tcs "<<(nEventTCS*nEventTCS)/denom<<endl;
	cout<<"nb event TCS "<<nEventTCS<<endl;
	cout<<"nb event CD "<<nCD<<endl;
	cout<<"nb event FD "<<nFD<<endl;
	cout<<"nb entries CheX "<<nentriesChe<<endl;

	cout<<"AfterCuts "<<AfterCuts<<endl;
	
	cout<<"nb JPSI "<<nbJPSI<<endl;
	cout<<" corr rad "<<corrrad<<endl;
	
	cout<<"first bin "<<((sum1*sum1)/sumw21)<<endl;;
	cout<<"second bin "<<((sum2*sum2)/sumw22)<<endl;;		
	
	//gROOT->ProcessLine(".q");

	gApplication->Terminate();

	return 0;
}
