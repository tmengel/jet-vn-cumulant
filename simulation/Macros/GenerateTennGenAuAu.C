#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TFile.h"	
#include "TRandom3.h"
#include "TTimeStamp.h"
#include "TMath.h"   
#include <fstream>
#include <iostream>
#include "TString.h"
#include "TTree.h"
#include "TSystem.h"
using namespace std;
using namespace ROOT;

const Int_t MAXPARTS = 2000;

Double_t dNdPhi(Double_t phiPart, Double_t pT, Double_t Psi1 , Double_t Psi2, Double_t Psi3, Double_t Psi4, Double_t v1,Double_t v2,Double_t v3,Double_t v4){
    return (1.0+2.0*(v1*TMath::Cos(phiPart-Psi1) + v2*TMath::Cos(2.0*(phiPart-Psi2)) + v3*TMath::Cos(3.0*(phiPart-Psi3))+ v4*TMath::Cos(4.0*(phiPart-Psi4)) ) );
}

Double_t MyIntegrandBG(Double_t *x, Double_t *p){

  Double_t x0 = x[0]; 
  
  Double_t mass     = p[0];
  Double_t pT       = p[1];
  Double_t beta_max = p[2];
  Double_t temp     = p[3];
  Double_t n      = p[4];

  // Keep beta within reasonable limits
  Double_t beta = beta_max * TMath::Power(x0, n);
  if (beta > 0.9999999999999999) beta = 0.9999999999999999;

  Double_t mT      = TMath::Sqrt(mass*mass+pT*pT);

  Double_t rho0   = TMath::ATanH(beta);  
  Double_t arg00 = pT*TMath::SinH(rho0)/temp;
  if (arg00 > 700.) arg00 = 700.; // avoid FPE
  Double_t arg01 = mT*TMath::CosH(rho0)/temp;
  Double_t f0 = x0*mT*TMath::BesselI0(arg00)*TMath::BesselK1(arg01);


  return f0;
}

Double_t MyStaticBGdNdPtTimesPt(Double_t *x, Double_t *p) {
  // implementation of BGBW (1/pT dNdpt)

  Double_t pT = x[0];;
  

  Double_t mass    = p[0];
  Double_t beta    = p[1];
  Double_t temp    = p[2];
  Double_t n       = p[3];
  Double_t norm    = p[4];

  static TF1 * fIntBG = 0;
  if(!fIntBG)
    fIntBG = new TF1 ("fIntBG", MyIntegrandBG, 0, 1, 5);

  fIntBG->SetParameters(mass, pT, beta, temp,n);
  Double_t result = fIntBG->Integral(0,1);
  return result*norm*pT;//*1e30;;
}

Double_t BlastWavedNdptTimesPt(Double_t *x, Double_t *p) {
  Double_t pT = x[0];
  
  Double_t mass    = p[0];
  Double_t A       = p[5];

  return MyStaticBGdNdPtTimesPt(x,p);
}

Double_t HarmonicFunction(Double_t pT, Int_t Cent, Int_t Vn, Int_t KF){
	

	Int_t part;
	if(pT > 4.0){  pT = 4.0;  }
	else if(pT < 0.5){ return 0.0; }

	if(KF ==211||KF==-211){ part = 0; }
	else if(KF ==321||KF==-321){ part = 1; }
	else if(KF ==2212||KF==-2212){ part = 2; }
	else{ part =0; }


	if(Cent==0){
		if(part==0){
			if(Vn==0){ return -0.0143392+0.0759773*pT - 0.0154024*pT*pT - 0.000915513*pT*pT*pT+0.000264875*pT*pT*pT*pT;}
			else if(Vn==1){ return -0.00558048+0.0136122*pT+0.0356842*pT*pT - 0.0164884*pT*pT*pT+0.00195171*pT*pT*pT*pT;}
			else if(Vn==2){ return -0.00256332 - 0.00483424*pT+0.0390473*pT*pT - 0.016438*pT*pT*pT+0.00205353*pT*pT*pT*pT;}
		}
		else if(part==1){
			if(Vn==0){ return -0.0287711+0.0777838*pT - 0.0110908*pT*pT - 0.00251265*pT*pT*pT+0.000435634*pT*pT*pT*pT;}
			else if(Vn==1){ return -0.00843333+0.00413273*pT+0.03708*pT*pT - 0.0140296*pT*pT*pT+0.00138853*pT*pT*pT*pT;}
			else if(Vn==2){ return -0.0147847+0.0238129*pT+0.00149144*pT*pT+0.00122299*pT*pT*pT - 0.000583605*pT*pT*pT*pT;}
		}
		else if(part==2){
			if(Vn==0){ return 0.0139405 - 0.0778337*pT+0.122313*pT*pT - 0.0441081*pT*pT*pT+0.00502149*pT*pT*pT*pT;}
			else if(Vn==1){ return 0.0204243 - 0.0952528*pT+0.118049*pT*pT - 0.0375978*pT*pT*pT+0.00388916*pT*pT*pT*pT;}
			else if(Vn==2){ return 0.0368641 - 0.132059*pT+0.140577*pT*pT - 0.0465704*pT*pT*pT+0.00527664*pT*pT*pT*pT;}
		}
	}
	else if(Cent==1){
		if(part==0){
			if(Vn==0){ return -0.0164604+0.119236*pT - 0.0140501*pT*pT - 0.00715798*pT*pT*pT+0.00137041*pT*pT*pT*pT;}
			else if(Vn==1){ return -0.00615225+0.0205719*pT+0.0373663*pT*pT - 0.0176272*pT*pT*pT+0.00201875*pT*pT*pT*pT;}
			else if(Vn==2){ return 0.00343526 - 0.0156758*pT+0.058631*pT*pT - 0.0234438*pT*pT*pT+0.00258536*pT*pT*pT*pT;}
		}
		else if(part==1){
			if(Vn==0){ return -0.0335949+0.108473*pT+0.00189956*pT*pT - 0.0119077*pT*pT*pT+0.00177362*pT*pT*pT*pT;}
			else if(Vn==1){ return -0.0141866+0.0239085*pT+0.0233996*pT*pT - 0.0081926*pT*pT*pT+0.000430152*pT*pT*pT*pT;}
			else if(Vn==2){ return 0.0178734 - 0.0725294*pT+0.105302*pT*pT - 0.0391775*pT*pT*pT+0.00466704*pT*pT*pT*pT;}
		}
		else if(part==2){
			if(Vn==0){ return 0.0147481 - 0.0885341*pT+0.16892*pT*pT - 0.0604128*pT*pT*pT+0.00661366*pT*pT*pT*pT;}
			else if(Vn==1){ return 0.020801 - 0.0910493*pT+0.118184*pT*pT - 0.0365487*pT*pT*pT+0.0035618*pT*pT*pT*pT;}
			else if(Vn==2){ return 0.0300511 - 0.108966*pT+0.122315*pT*pT - 0.0365423*pT*pT*pT+0.00350489*pT*pT*pT*pT;}
		}
	}
	else if(Cent==2){
		if(part==0){
			if(Vn==0){ return -0.0220529+0.172125*pT - 0.0353618*pT*pT - 0.003559*pT*pT*pT+0.00113968*pT*pT*pT*pT;}
			else if(Vn==1){ return -0.0066372+0.0262161*pT+0.0372216*pT*pT - 0.0187145*pT*pT*pT+0.00228567*pT*pT*pT*pT;}
			else if(Vn==2){ return 0.000152642+0.00135534*pT+0.0523496*pT*pT - 0.0225954*pT*pT*pT+0.0025451*pT*pT*pT*pT;}
		}
		else if(part==1){
			if(Vn==0){ return -0.0424241+0.152629*pT - 0.00506494*pT*pT - 0.0151633*pT*pT*pT+0.00254353*pT*pT*pT*pT;}
			else if(Vn==1){ return -0.0149325+0.0253627*pT+0.0329371*pT*pT - 0.0153877*pT*pT*pT+0.00170996*pT*pT*pT*pT;}
			else if(Vn==2){ return -0.0171898+0.0261749*pT+0.032913*pT*pT - 0.0180592*pT*pT*pT+0.00240376*pT*pT*pT*pT;}
		}
		else if(part==2){
			if(Vn==0){ return 0.0128407 - 0.0812974*pT+0.196424*pT*pT - 0.0729275*pT*pT*pT+0.0081403*pT*pT*pT*pT;}
			else if(Vn==1){ return 0.0216277 - 0.0905268*pT+0.125852*pT*pT - 0.0410326*pT*pT*pT+0.00433817*pT*pT*pT*pT;}
			else if(Vn==2){ return 0.0296393 - 0.113592*pT+0.137947*pT*pT - 0.0424535*pT*pT*pT+0.00422479*pT*pT*pT*pT;}
		}
	}
	else if(Cent==3){
		if(part==0){
			if(Vn==0){ return -0.0273469+0.215291*pT - 0.0580156*pT*pT+0.0015503*pT*pT*pT+0.00068957*pT*pT*pT*pT;}
			else if(Vn==1){ return -0.00634738+0.0244379*pT+0.0472794*pT*pT - 0.0265474*pT*pT*pT+0.00383202*pT*pT*pT*pT;}
			else if(Vn==2){ return 0.00529299 - 0.0155944*pT+0.0851034*pT*pT - 0.0399046*pT*pT*pT+0.00537977*pT*pT*pT*pT;}
		}
		else if(part==1){
			if(Vn==0){ return -0.0457415+0.184931*pT - 0.0184578*pT*pT - 0.0130774*pT*pT*pT+0.00241422*pT*pT*pT*pT;}
			else if(Vn==1){ return 0.00059818 - 0.0174573*pT+0.0752039*pT*pT - 0.0318181*pT*pT*pT+0.00386052*pT*pT*pT*pT;}
			else if(Vn==2){ return 0.00319935 - 0.0357498*pT+0.0956003*pT*pT - 0.0389201*pT*pT*pT+0.0046787*pT*pT*pT*pT;}
		}
		else if(part==2){
			if(Vn==0){ return 0.00914554 - 0.0597874*pT+0.203465*pT*pT - 0.0797661*pT*pT*pT+0.00929514*pT*pT*pT*pT;}
			else if(Vn==1){ return 0.0304227 - 0.111558*pT+0.150866*pT*pT - 0.0511995*pT*pT*pT+0.00556649*pT*pT*pT*pT;}
			else if(Vn==2){ return -0.0025491 - 0.0227755*pT+0.0628781*pT*pT - 0.0165041*pT*pT*pT+0.00111185*pT*pT*pT*pT;}
		}
	}
	else if(Cent==4){
		if(part==0){
			if(Vn==0){ return -0.0300557+0.23959*pT - 0.0712208*pT*pT+0.004233*pT*pT*pT+0.000504197*pT*pT*pT*pT;}
			else if(Vn==1){ return -0.0047109+0.0195728*pT+0.0522525*pT*pT - 0.0282469*pT*pT*pT+0.00377098*pT*pT*pT*pT;}
			else if(Vn==2){ return -0.0132215+0.0468*pT+0.0341852*pT*pT - 0.0206421*pT*pT*pT+0.00294137*pT*pT*pT*pT;}
		}
		else if(part==1){
			if(Vn==0){ return -0.0425067+0.191418*pT - 0.0147714*pT*pT - 0.0177701*pT*pT*pT+0.00341417*pT*pT*pT*pT;}
			else if(Vn==1){ return -0.000136675 - 0.0175618*pT+0.0863983*pT*pT - 0.0430817*pT*pT*pT+0.00620464*pT*pT*pT*pT;}
			else if(Vn==2){ return -0.0570229+0.17675*pT - 0.123802*pT*pT+0.0478088*pT*pT*pT - 0.00638515*pT*pT*pT*pT;}
		}
		else if(part==2){
			if(Vn==0){ return 0.0054852 - 0.0327023*pT+0.19693*pT*pT - 0.0815048*pT*pT*pT+0.0098101*pT*pT*pT*pT;}
			else if(Vn==1){ return 0.0109575 - 0.0600514*pT+0.115052*pT*pT - 0.0418587*pT*pT*pT+0.00470501*pT*pT*pT*pT;}
			else if(Vn==2){ return 0.04566 - 0.148386*pT+0.193706*pT*pT - 0.0675996*pT*pT*pT+0.00792379*pT*pT*pT*pT;}
		}
	}
	else if(Cent==5){
		if(part==0){
			if(Vn==0){ return -0.0327924+0.250176*pT - 0.0765101*pT*pT+0.00390845*pT*pT*pT+0.000819225*pT*pT*pT*pT;}
			else if(Vn==1){ return -0.003819+0.0112537*pT+0.0588299*pT*pT - 0.0333377*pT*pT*pT+0.0046983*pT*pT*pT*pT;}
			else if(Vn==2){ return -0.00266408+0.00640717*pT+0.023585*pT*pT - 0.0121294*pT*pT*pT+0.00172664*pT*pT*pT*pT;}
		}
		else if(part==1){
			if(Vn==0){ return -0.0131631+0.0325158*pT - 0.00707803*pT*pT+0.000728541*pT*pT*pT - 8.91768e-05*pT*pT*pT*pT;}
			else if(Vn==1){ return -0.0571195+0.249406*pT - 0.074045*pT*pT+0.00463722*pT*pT*pT+0.000540562*pT*pT*pT*pT;}
			else if(Vn==2){ return -0.0143528+0.0402737*pT+0.0106742*pT*pT - 0.00873702*pT*pT*pT+0.000978765*pT*pT*pT*pT;}
		}
		else if(part==2){
			if(Vn==0){ return -0.00897794+0.0202506*pT+0.159824*pT*pT - 0.0719297*pT*pT*pT+0.00894275*pT*pT*pT*pT;}
			else if(Vn==1){ return -0.00854808+0.0237419*pT+0.00737678*pT*pT+0.00711372*pT*pT*pT - 0.00275382*pT*pT*pT*pT;}
			else if(Vn==2){ return 0.0; } //return 6.95289e-310+4.67364e-310*pT+4.67364e-310*pT*pT+6.95289e-310*pT*pT*pT+5.05923e-321*pT*pT*pT*pT;}
		}
	}
}

Int_t SpectraCentBin(Int_t multiplicity){

        if(multiplicity>=500){return 0;}
        else if(multiplicity>=345){return 1;}
        else if(multiplicity>=154){return 2;}
        else if(multiplicity>=58){return 3;}
        else{return 6;}
}

Int_t HarmonicCentBin(Int_t multiplicity){

        if(multiplicity>=500){return 0;}
        else if(multiplicity>=345){return 1;}
        else if(multiplicity>=234 ){return 2;}
        else if(multiplicity>=154){return 3;}
        else if(multiplicity>=99){return 4;}
        else if(multiplicity>=58){return 5;}
        else{return 6;}
}

const Double_t Multiplicity_Params[10] = {0.0569104,-0.00129485,1.45861e-05,-8.83828e-08,3.10488e-10,-6.54566e-13,8.25612e-16,-5.9263e-19,2.10337e-22,-2.40376e-26};

const Double_t MASS[4][6] ={{ 0.13957 , 0.13957 , 0.49368 , 0.49368 , 0.93827 , 0.93827 }  , { 0.13957 , 0.13957 , 0.49368 , 0.49368 , 0.93827 , 0.93827 }  , { 0.13957 , 0.13957 , 0.49368 , 0.49368 , 0.93827 , 0.93827 }  , { 0.13957 , 0.13957 , 0.49368 , 0.49368 , 0.93827 , 0.93827 }  };
const Double_t BETA[4][6] ={{ 0.759913 , 0.760621 , 0.7319 , 0.732047 , 0.648344 , 0.680417 }  , { 0.765472 , 0.763144 , 0.739863 , 0.739893 , 0.641706 , 0.684773 }  , { 0.772956 , 0.773497 , 0.754914 , 0.756607 , 0.661272 , 0.692553 }  , { 0.913187 , 0.788666 , 0.780603 , 0.780076 , 0.681797 , 0.696818 }  };
const Double_t TEMP[4][6] ={{ 0.175254 , 0.17647 , 0.163823 , 0.167908 , 0.25008 , 0.256333 }  , { 0.175317 , 0.17444 , 0.166239 , 0.170082 , 0.247217 , 0.256199 }  , { 0.173122 , 0.17437 , 0.164637 , 0.168712 , 0.240833 , 0.249546 }  , { 0.0835104 , 0.169124 , 0.16093 , 0.163796 , 0.218115 , 0.22822 }  };
const Double_t NPAR[4][6] ={{ 32.371 , 33.1419 , 11.3702 , 12.1617 , 44.6847 , 74.0065 }  , { 30.7 , 28.7665 , 11.8083 , 12.3836 , 35.0342 , 65.2486 }  , { 28.408 , 29.0054 , 12.4268 , 13.1964 , 34.7224 , 54.8703 }  , { 27.5124 , 27.7284 , 14.806 , 14.6286 , 27.4966 , 36.6379 }  };
const Double_t NORM[4][6] ={{ 10211.8 , 9733.9 , 7669.71 , 6958.68 , 714.188 , 833.707 }  , { 6923.66 , 7024.51 , 4679.75 , 4275.71 , 521.372 , 559.871 }  , { 4030.29 , 3828.84 , 2667.37 , 2413.47 , 341.159 , 360.574 }  , { 109206 , 1663.87 , 1117.2 , 1030.07 , 254.529 , 241.505 }  };

const Double_t piplusRatios[4] = {0.396738386,0.400412797,0.400958581,0.40230616};
const Double_t piminusRatios[4] = {0.396738386,0.400412797,0.400958581,0.40230616};
const Double_t kplusRatios[4] = {0.063108127,0.061919505,0.060984334,0.059731231};
const Double_t kminusRatios[4] = {0.06118953,0.059236326, 0.058838257,0.057484421};
const Double_t proRatios[4] = {0.043099904,0.041486068,0.042385006,0.042604604};
const Double_t pbarRatios[4] = {0.03295875,0.032404541,0.033371486,0.034295646};

Int_t GenerateTennGenAuAu(Int_t nevent, Int_t centbin, Double_t etaRange, Int_t prefix){

    Int_t nparts;
    Double_t particle_px[MAXPARTS], particle_py[MAXPARTS], particle_pz[MAXPARTS], particle_e[MAXPARTS];
    Int_t particle_id[MAXPARTS];
    Double_t event_psi1, event_psi2, event_psi3, event_psi4;

    TH2D *v2_pi = new TH2D("v2_pi","v2_pi",50,0,5,50,0,1);
    TH2D *v2_k = new TH2D("v2_k","v2_k",50,0,5,50,0,1);
    TH2D *v2_p = new TH2D("v2_p","v2_p",50,0,5,50,0,1);
    TH2D *v3_pi = new TH2D("v3_pi","v3_pi",50,0,5,50,0,1);
    TH2D *v3_k = new TH2D("v3_k","v3_k",50,0,5,50,0,1);
    TH2D *v3_p = new TH2D("v3_p","v3_p",50,0,5,50,0,1);
    TH2D *v4_pi = new TH2D("v4_pi","v4_pi",50,0,5,50,0,1);
    TH2D *v4_k = new TH2D("v4_k","v4_k",50,0,5,50,0,1);
    TH2D *v4_p = new TH2D("v4_p","v4_p",50,0,5,50,0,1);

    TF1 *pTdist_piPlus;
    TF1 *pTdist_piMinus;
    TF1 *pTdist_kPlus;
    TF1 *pTdist_kMinus;
    TF1 *pTdist_p;
    TF1 *pTdist_pbar;
    TF1 *multiplicty_distro = new TF1("multiplicity_distro","pol9",0,700);

    multiplicty_distro->SetNpx (10000);
    for(Int_t i=0;i<10;i++){
        multiplicty_distro->SetParameter(i,Multiplicity_Params[i]);
    }

    pTdist_piPlus = new TF1("pTdist_piPlus",BlastWavedNdptTimesPt ,0,100,5); 
    pTdist_piPlus->SetNpx (1000);
    pTdist_piMinus = new TF1("pTdist_piMinus",BlastWavedNdptTimesPt ,0,100,5); 
    pTdist_piMinus->SetNpx (1000);  
    pTdist_kPlus = new TF1("pTdist_kPlus",BlastWavedNdptTimesPt ,0,100,5); 
    pTdist_kPlus->SetNpx (1000);  
    pTdist_kMinus = new TF1("pTdist_kMinus",BlastWavedNdptTimesPt ,0,100,5); 
    pTdist_kMinus->SetNpx (1000);  
    pTdist_p = new TF1("pTdist_p",BlastWavedNdptTimesPt ,0,100,5); 
    pTdist_p->SetNpx (1000);  
    pTdist_pbar = new TF1("pTdist_pbar",BlastWavedNdptTimesPt ,0,100,5); 
    pTdist_pbar->SetNpx (1000);  
    
    TRandom3 *etadist;
    TRandom3 *Harmonics_Phi_Dist_Rand;
    TRandom3 *Psi_1;
    TRandom3 *Psi_3;
    TRandom3 *Multi_Rand;
    TRandom3 *pT_Rand;

    TTimeStamp *timestamp = new TTimeStamp();
    UInt_t Seed = UInt_t(timestamp->GetSec());
    UInt_t seed_eta = Seed * 7;
    etadist = new TRandom3(seed_eta);
    UInt_t seed_psi_1 = (Seed*2)+3;
    Psi_1 = new TRandom3(seed_psi_1);
    UInt_t seed_psi_3 = (Seed*3)+2;
    Psi_3 = new TRandom3(seed_psi_3);
    UInt_t seed_uniform_random_phi = (Seed*3)*4;
    Harmonics_Phi_Dist_Rand = new TRandom3(seed_uniform_random_phi);
    UInt_t seed_multi = (Seed*2)*3;
    Multi_Rand = new TRandom3(seed_multi);
    UInt_t seed_pT = (Seed*2)*3;
    pT_Rand = new TRandom3(seed_pT);    


    // Output file
    TString datadir = "../root-files/AuAu/";
    if(gSystem->AccessPathName(datadir)) gSystem->mkdir(datadir);

    if(centbin == 0) datadir += "0-10/";
    else if(centbin == 1) datadir += "10-20/";
    else if(centbin == 2) datadir += "20-40/";
    else if(centbin == 3) datadir += "40-60/";
    else {
        cout << "Invalid centrality bin" << endl;
        return 0;
    }
    if(gSystem->AccessPathName(datadir)) gSystem->mkdir(datadir);

    pTdist_piPlus->SetParameter(0,MASS[centbin][1]);
    pTdist_piPlus->SetParameter(1,BETA[centbin][1]);
    pTdist_piPlus->SetParameter(2,TEMP[centbin][1]);
    pTdist_piPlus->SetParameter(3,NPAR[centbin][1]);
    pTdist_piPlus->SetParameter(4,NORM[centbin][1]);

    pTdist_piMinus->SetParameter(0,MASS[centbin][0]);
    pTdist_piMinus->SetParameter(1,BETA[centbin][0]);
    pTdist_piMinus->SetParameter(2,TEMP[centbin][0]);
    pTdist_piMinus->SetParameter(3,NPAR[centbin][0]);
    pTdist_piMinus->SetParameter(4,NORM[centbin][0]);

    pTdist_kPlus->SetParameter(0,MASS[centbin][3]);
    pTdist_kPlus->SetParameter(1,BETA[centbin][3]);
    pTdist_kPlus->SetParameter(2,TEMP[centbin][3]);
    pTdist_kPlus->SetParameter(3,NPAR[centbin][3]);
    pTdist_kPlus->SetParameter(4,NORM[centbin][3]);

    pTdist_kMinus->SetParameter(0,MASS[centbin][2]);
    pTdist_kMinus->SetParameter(1,BETA[centbin][2]);
    pTdist_kMinus->SetParameter(2,TEMP[centbin][2]);
    pTdist_kMinus->SetParameter(3,NPAR[centbin][2]);
    pTdist_kMinus->SetParameter(4,NORM[centbin][2]);

    pTdist_p->SetParameter(0,MASS[centbin][5]);
    pTdist_p->SetParameter(1,BETA[centbin][5]);
    pTdist_p->SetParameter(2,TEMP[centbin][5]);
    pTdist_p->SetParameter(3,NPAR[centbin][5]);
    pTdist_p->SetParameter(4,NORM[centbin][5]);

    pTdist_pbar->SetParameter(0,MASS[centbin][4]);
    pTdist_pbar->SetParameter(1,BETA[centbin][4]);
    pTdist_pbar->SetParameter(2,TEMP[centbin][4]);
    pTdist_pbar->SetParameter(3,NPAR[centbin][4]);
    pTdist_pbar->SetParameter(4,NORM[centbin][4]);

    TString filename = Form("%s200GeV_AuAu_ptbin%d.root",datadir.Data(),prefix);
    TFile *outFile = new TFile(filename.Data(),"RECREATE");

    TTree* outTree = new TTree("tree", "tree");
    TTree* eventInfo = new TTree("eventInfo","eventInfo");

    eventInfo->Branch("nevent",&nevent,"nevent/I");
    eventInfo->Branch("centbin",&centbin,"centbin/I");
    eventInfo->Branch("etarange",&etaRange,"etarange/D");


    outTree->Branch("nparts",&nparts,"nparts/I");
    outTree->Branch("particle_px", particle_px, "particle_px[nparts]/D");
    outTree->Branch("particle_py", particle_py, "particle_py[nparts]/D");
    outTree->Branch("particle_pz", particle_pz, "particle_pz[nparts]/D");
    outTree->Branch("particle_e", particle_e, "particle_e[nparts]/D");
    outTree->Branch("particle_id", particle_id, "particle_id[nparts]/I");
    outTree->Branch("psi1", &event_psi1, "psi1/D");
    outTree->Branch("psi2", &event_psi2, "psi2/D");
    outTree->Branch("psi3", &event_psi3, "psi3/D");
    outTree->Branch("psi4", &event_psi4, "psi4/D");

    cout << "Starting event loop" << endl;
    for (Int_t ievent=0; ievent<nevent; ievent++){
    
    
        Double_t Psi_1_event = Psi_1->Uniform( 0 , 2.0*TMath::Pi());
        Double_t Psi_2_event = 0.0;
        Double_t Psi_3_event = Psi_3->Uniform( 0 , 2.0*TMath::Pi());           
        Double_t Psi_4_event = 0.0;

      
        Int_t particle_multiplicity = -1;
        // cout << "Starting multiplicity loop" << endl;
        while (1){
            particle_multiplicity = Int_t((etaRange/0.5)*multiplicty_distro->GetRandom(Multi_Rand));
            if(SpectraCentBin(particle_multiplicity) == centbin) break;
        }
        Int_t harmonic_centbin = HarmonicCentBin(particle_multiplicity);

        Int_t particle_yeilds[6] = {0,0,0,0,0,0};
        particle_yeilds[0] = Int_t(piminusRatios[centbin]*particle_multiplicity);
        particle_yeilds[1] = Int_t(piplusRatios[centbin]*particle_multiplicity);
        particle_yeilds[2] = Int_t(kminusRatios[centbin]*particle_multiplicity);
        particle_yeilds[3] = Int_t(kplusRatios[centbin]*particle_multiplicity);
        particle_yeilds[4] = Int_t(pbarRatios[centbin]*particle_multiplicity);
        particle_yeilds[5] = Int_t(proRatios[centbin]*particle_multiplicity);
        
        nparts = 0;
        // if(ievent%10 == 0) cout << "Event " << ievent << endl;
        // cout << "Starting particle loop" << endl;
        for(Int_t n = 0; n<6; n++)
        {
            for(Int_t ipart = 0; ipart<particle_yeilds[n];ipart++)
            {
                
                Int_t KF = 0;
                if (n==0) KF = 211;
                if (n==1) KF = -211;
                if (n==2) KF = 321;
                if (n==3) KF = -321;
                if (n==4) KF = 2212;
                if (n==5) KF = -2212;

                Double_t particle_mass = 0.0;
                if(n == 0) particle_mass = 0.139570;
                else if(n==1) particle_mass = 0.139570;
                else if(n==2) particle_mass = 0.493677;
                else if(n==3) particle_mass = 0.493677;
                else if(n==4) particle_mass = 0.938272;
                else if(n==5) particle_mass = 0.938272;

                Double_t pT = 0.0;
                if (n==0) pT = pTdist_piMinus->GetRandom(pT_Rand);
                if (n==1) pT = pTdist_piPlus->GetRandom(pT_Rand);
                if (n==2) pT = pTdist_kMinus->GetRandom(pT_Rand);
                if (n==3) pT = pTdist_kPlus->GetRandom(pT_Rand);
                if (n==4) pT = pTdist_pbar->GetRandom(pT_Rand);
                if (n==5) pT = pTdist_p->GetRandom(pT_Rand);
                

                Double_t eta = etadist->Uniform(-1.0*etaRange , 1.0*etaRange);


                Double_t phi = 0.0;
                Double_t v2_particle = HarmonicFunction(pT,harmonic_centbin,0,KF);
                Double_t v3_particle = HarmonicFunction(pT,harmonic_centbin,1,KF);
                Double_t v4_particle = HarmonicFunction(pT,harmonic_centbin,2,KF);
                Double_t v1_particle = 0.02*v2_particle;

                if(n ==0 || n==1){
                    v2_pi->Fill(pT,v2_particle);
                    v3_pi->Fill(pT,v3_particle);
                    v4_pi->Fill(pT,v4_particle);
                } 
                if(n ==2 || n==3){
                    v2_k->Fill(pT,v2_particle);
                    v3_k->Fill(pT,v3_particle);
                    v4_k->Fill(pT,v4_particle);
                }
                if(n ==4 || n==5){
                    v2_p->Fill(pT,v2_particle);
                    v3_p->Fill(pT,v3_particle);
                    v4_p->Fill(pT,v4_particle);
                }


                Double_t Max_dNdPhi = (1.0+2.0*(TMath::Abs(v1_particle)+TMath::Abs(v2_particle)+TMath::Abs(v3_particle)+TMath::Abs(v4_particle)));
                Int_t CHECK = 0.0;
                while(CHECK!=1){
                    Double_t dndphi = Harmonics_Phi_Dist_Rand->Uniform(0.0,Max_dNdPhi);
                    Double_t test_phi = Harmonics_Phi_Dist_Rand->Uniform(0.0,2.0*TMath::Pi());
                    Double_t  Value = dNdPhi(test_phi,pT,Psi_1_event,Psi_2_event,Psi_3_event,Psi_4_event,v1_particle,v2_particle,v3_particle,v4_particle);
                    if(dndphi < Value) {
                        CHECK = 1;
                        phi = test_phi;
                    }
                }

                Double_t px = pT*TMath::Cos(phi);
                Double_t py = pT*TMath::Sin(phi);
                Double_t pz = pT*TMath::SinH(eta);
                Double_t E = TMath::Sqrt(px*px+py*py+pz*pz+particle_mass*particle_mass);

                particle_px[nparts] = px;
                particle_py[nparts] = py;
                particle_pz[nparts] = pz;
                particle_e[nparts] = E;
                particle_id[nparts] = KF;

                nparts++;



            } // end species loop

        } // end particle loop

        event_psi1 = Psi_1_event;
        event_psi2 = Psi_2_event;
        event_psi3 = Psi_3_event;
        event_psi4 = Psi_4_event;

        outTree->Fill();

        for (Int_t j=0; j< nparts; j++){
            particle_px[j] = 0.0;
            particle_py[j] = 0.0;
            particle_pz[j] = 0.0;
            particle_e[j] = 0.0;
            particle_id[j] = 0;  
        }
    
 
    }
    eventInfo->Fill();
    v2_pi->Write();
    v3_pi->Write();
    v4_pi->Write();
    v2_k->Write();
    v3_k->Write();
    v4_k->Write();
    v2_p->Write();
    v3_p->Write();
    v4_p->Write();

    outFile->Write();
    outFile->Close();
    delete outFile;   
    cout << "Done" << endl;
    return 0;

}

Int_t main(Int_t argc, char** argv){
    
    Int_t nevent = 100;
    Int_t centbin = 0;
    Double_t etaRange = 1.1;
    Int_t prefix = 0;
    if(argc <= 4){
        cout << "Usage: ./GenerateTennGenAuAu <nevents> <centrality> <etarange> <prefix>" << endl;
        return 1;
    }
    nevent = atoi(argv[1]);
    centbin = atoi(argv[2]);
    etaRange = atof(argv[3]);
    prefix = atoi(argv[4]);
    cout << "Generating " << nevent << " events with " << centbin << " cent bin and " << etaRange << " eta range" << endl;
    cout << "Prefix: " << prefix << endl;
    return GenerateTennGenAuAu(nevent, centbin, etaRange, prefix); 
}
