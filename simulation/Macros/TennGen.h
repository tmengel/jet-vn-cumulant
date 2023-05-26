#include <iostream>
#include <vector>

#include <TString.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TRandom3.h>
#include <TTimeStamp.h>
#include <TMath.h>

using namespace std;
using namespace ROOT;

struct TGParticle {
    Float_t px;
    Float_t py;
    Float_t pz;
    Float_t E;
    Float_t mass;
    Float_t pT;
    Float_t eta;
    Float_t phi;
    Int_t charge;
    Int_t pdg;
};

struct TGEvent {
    std::vector<TGParticle> particles;
    Int_t nparts;
    Double_t etaRange;
    Int_t centBin;
    Int_t collEn;
    Float_t psi1;
    Float_t psi2;
    Float_t psi3;
    Float_t psi4;
    Float_t psi5;
};

class TennGen
{
    private:
  
    const Int_t availible_collision_energies[2] = {200, 2760};
    const vector<TString> availible_collision_species = {"AuAu", "PbPb"};
    const Int_t availible_centraility_bins_auau[4] ={0,1,2,3};
    const Int_t availible_centraility_bins_pbpb[6] ={0,1,2,3,4,5};
    const Double_t max_eta_range = 1.1;
    const vector<TString> availible_centraility_strings_auau = {"0-10%", "10-20%", "20-40%", "40-60%"};
    const vector<TString> availible_centraility_strings_pbpb = {"0-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%"};
    const Int_t availible_harmonics_auau[4] = {1,2,3,4};
    const Int_t availible_harmonics_pbpb[5] = {1,2,3,4,5};


    Int_t fCollEn;
    Int_t fCentBin;
    Double_t fEtaRange;

    Int_t doVn[5] = {0,0,0,0,0};
    Int_t fixPsi[5] = {0,1,0,1,0};
    Float_t fixed_psi[5] = {0.0,0.0,0.0,0.0,0.0};
    

    TF1 *v1_pi;
    TF1 *v2_pi;
    TF1 *v3_pi;
    TF1 *v4_pi;
    TF1 *v5_pi;
    TF1 *v1_K;
    TF1 *v2_K;
    TF1 *v3_K;
    TF1 *v4_K;
    TF1 *v5_K;
    TF1 *v1_P;
    TF1 *v2_P;
    TF1 *v3_P;
    TF1 *v4_P;
    TF1 *v5_P;

    TF1 *pTdist_piPlus;
    TF1 *pTdist_piMinus;
    TF1 *pTdist_kPlus;
    TF1 *pTdist_kMinus;
    TF1 *pTdist_p;
    TF1 *pTdist_pbar;
    TF1 *pTdist_piZero;
    
    TF1 *multiplicty_distro;

    TTimeStamp *timestamp;

    TRandom3 *etadist;
    TRandom3 *Harmonics_Phi_Dist_Rand;
    TRandom3 *Psi_1;
    TRandom3 *Psi_2;
    TRandom3 *Psi_3;
    TRandom3 *Psi_4;
    TRandom3 *Psi_5;    
    TRandom3 *Multi_Rand;
    TRandom3 *pT_Rand;

    TFile *fout;

public:

    TennGen(Int_t collen, Int_t centbin, Double_t etaRange);

    ~TennGen();

    static Double_t MyIntegrandBG(Double_t *x, Double_t *p){

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

    static Double_t MyStaticBGdNdPtTimesPt(Double_t *x, Double_t *p) {
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

    static Double_t BlastWavedNdptTimesPt(Double_t *x, Double_t *p) {
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

    void CheckInputs(){
         
        Int_t valid_collen_flag = 0;
        Int_t valid_centbin_flag = 0;
        Int_t valid_etarange_flag = 0;
        for (Int_t i = 0; i < 2; i++) {
            if (fCollEn == availible_collision_energies[i]) {
                valid_collen_flag = 1;
            }
        }

        if (!valid_collen_flag) {
            cout << "Invalid collision energy. Please choose from 200 or 2760 GeV." << endl;
            exit(1);
        }

        if(fCollEn == 200){
            for (Int_t i = 0; i < 4; i++) {
                if (fCentBin == availible_centraility_bins_auau[i]) {
                    valid_centbin_flag = 1;
                }
            }
        }
        
        else if(fCollEn == 2760){
            for (Int_t i = 0; i < 6; i++) {
                if (fCentBin == availible_centraility_bins_pbpb[i]) {
                    valid_centbin_flag = 1;
                }
            }
        }
    
        if(!valid_centbin_flag){
            cout << "Invalid centrality bin. Please choose from the following:" << endl;
            if(fCollEn == 200){
                for (Int_t i = 0; i < 4; i++) {
                    cout << "Cent " << availible_centraility_strings_auau[i] << " : " << availible_centraility_bins_auau[i] << endl;
                }
            }
            else if(fCollEn == 2760){
                for (Int_t i = 0; i < 6; i++) {
                    cout << "Cent " << availible_centraility_strings_pbpb[i] << " : " << availible_centraility_bins_pbpb[i] << endl;
                }
            }
            exit(1);
        }

        if(fEtaRange < 0 || fEtaRange > max_eta_range){
            cout << "Invalid eta range. Please choose from 0 to 1.1" << endl;
            exit(1);
        }
        
    }

    void InitAuAu(){

        const Double_t Multiplicity_Params[10] = {0.0569104,-0.00129485,1.45861e-05,-8.83828e-08,3.10488e-10,-6.54566e-13,8.25612e-16,-5.9263e-19,2.10337e-22,-2.40376e-26};
        const Double_t MASS[4][6] ={{ 0.13957 , 0.13957 , 0.49368 , 0.49368 , 0.93827 , 0.93827 }  , { 0.13957 , 0.13957 , 0.49368 , 0.49368 , 0.93827 , 0.93827 }  , { 0.13957 , 0.13957 , 0.49368 , 0.49368 , 0.93827 , 0.93827 }  , { 0.13957 , 0.13957 , 0.49368 , 0.49368 , 0.93827 , 0.93827 }  };
        const Double_t BETA[4][6] ={{ 0.759913 , 0.760621 , 0.7319 , 0.732047 , 0.648344 , 0.680417 }  , { 0.765472 , 0.763144 , 0.739863 , 0.739893 , 0.641706 , 0.684773 }  , { 0.772956 , 0.773497 , 0.754914 , 0.756607 , 0.661272 , 0.692553 }  , { 0.913187 , 0.788666 , 0.780603 , 0.780076 , 0.681797 , 0.696818 }  };
        const Double_t TEMP[4][6] ={{ 0.175254 , 0.17647 , 0.163823 , 0.167908 , 0.25008 , 0.256333 }  , { 0.175317 , 0.17444 , 0.166239 , 0.170082 , 0.247217 , 0.256199 }  , { 0.173122 , 0.17437 , 0.164637 , 0.168712 , 0.240833 , 0.249546 }  , { 0.0835104 , 0.169124 , 0.16093 , 0.163796 , 0.218115 , 0.22822 }  };
        const Double_t NPAR[4][6] ={{ 32.371 , 33.1419 , 11.3702 , 12.1617 , 44.6847 , 74.0065 }  , { 30.7 , 28.7665 , 11.8083 , 12.3836 , 35.0342 , 65.2486 }  , { 28.408 , 29.0054 , 12.4268 , 13.1964 , 34.7224 , 54.8703 }  , { 27.5124 , 27.7284 , 14.806 , 14.6286 , 27.4966 , 36.6379 }  };
        const Double_t NORM[4][6] ={{ 10211.8 , 9733.9 , 7669.71 , 6958.68 , 714.188 , 833.707 }  , { 6923.66 , 7024.51 , 4679.75 , 4275.71 , 521.372 , 559.871 }  , { 4030.29 , 3828.84 , 2667.37 , 2413.47 , 341.159 , 360.574 }  , { 109206 , 1663.87 , 1117.2 , 1030.07 , 254.529 , 241.505 }  };

        multiplicty_distro = new TF1("multiplicity_distro","pol9",0,670);

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

        timestamp = new TTimeStamp();
        UInt_t Seed = UInt_t(timestamp->GetSec());
        UInt_t seed_eta = Seed * 7;
        etadist = new TRandom3(seed_eta);
        UInt_t seed_psi_1 = (Seed*2)+3;
        Psi_1 = new TRandom3(seed_psi_1);
        UInt_t seed_psi_2 = (Seed*2)+1;
        Psi_2 = new TRandom3(seed_psi_2);
        UInt_t seed_psi_3 = (Seed*3)+2;
        Psi_3 = new TRandom3(seed_psi_3);
        UInt_t seed_psi_4 = (Seed*3)+4;
        Psi_4 = new TRandom3(seed_psi_4);
        UInt_t seed_uniform_random_phi = (Seed*3)*4;
        Harmonics_Phi_Dist_Rand = new TRandom3(seed_uniform_random_phi);
        UInt_t seed_multi = (Seed*2)*3;
        Multi_Rand = new TRandom3(seed_multi);
        UInt_t seed_pT = (Seed*2)*3;
        pT_Rand = new TRandom3(seed_pT);   
        
        pTdist_piPlus->SetParameter(0,MASS[fCentBin][1]);
        pTdist_piPlus->SetParameter(1,BETA[fCentBin][1]);
        pTdist_piPlus->SetParameter(2,TEMP[fCentBin][1]);
        pTdist_piPlus->SetParameter(3,NPAR[fCentBin][1]);
        pTdist_piPlus->SetParameter(4,NORM[fCentBin][1]);

        pTdist_piMinus->SetParameter(0,MASS[fCentBin][0]);
        pTdist_piMinus->SetParameter(1,BETA[fCentBin][0]);
        pTdist_piMinus->SetParameter(2,TEMP[fCentBin][0]);
        pTdist_piMinus->SetParameter(3,NPAR[fCentBin][0]);
        pTdist_piMinus->SetParameter(4,NORM[fCentBin][0]);

        pTdist_kPlus->SetParameter(0,MASS[fCentBin][3]);
        pTdist_kPlus->SetParameter(1,BETA[fCentBin][3]);
        pTdist_kPlus->SetParameter(2,TEMP[fCentBin][3]);
        pTdist_kPlus->SetParameter(3,NPAR[fCentBin][3]);
        pTdist_kPlus->SetParameter(4,NORM[fCentBin][3]);

        pTdist_kMinus->SetParameter(0,MASS[fCentBin][2]);
        pTdist_kMinus->SetParameter(1,BETA[fCentBin][2]);
        pTdist_kMinus->SetParameter(2,TEMP[fCentBin][2]);
        pTdist_kMinus->SetParameter(3,NPAR[fCentBin][2]);
        pTdist_kMinus->SetParameter(4,NORM[fCentBin][2]);

        pTdist_p->SetParameter(0,MASS[fCentBin][5]);
        pTdist_p->SetParameter(1,BETA[fCentBin][5]);
        pTdist_p->SetParameter(2,TEMP[fCentBin][5]);
        pTdist_p->SetParameter(3,NPAR[fCentBin][5]);
        pTdist_p->SetParameter(4,NORM[fCentBin][5]);

        pTdist_pbar->SetParameter(0,MASS[fCentBin][4]);
        pTdist_pbar->SetParameter(1,BETA[fCentBin][4]);
        pTdist_pbar->SetParameter(2,TEMP[fCentBin][4]);
        pTdist_pbar->SetParameter(3,NPAR[fCentBin][4]);
        pTdist_pbar->SetParameter(4,NORM[fCentBin][4]);

        doVn[0] = 1;
        doVn[1] = 1;
        doVn[2] = 1;
        doVn[3] = 1;
        doVn[4] = 0;

        fixPsi[0] = 0;
        fixPsi[1] = 1;
        fixPsi[2] = 0;
        fixPsi[3] = 1;
        fixPsi[4] = 1;
        
        fixed_psi[1] = 0.0;
        fixed_psi[3] = 0.0;
        fixed_psi[4] = 0.0;

    }

    void InitPbPb(){
        
        const Double_t v2_pi_params[6][6] = { {-6.29928e-03 , 6.29213e-02 , -1.69603e-02 , 1.24408e-03 , 0.0 , 5.6864 } , { -8.15353e-03 ,  9.46349e-02 , -2.27846e-02 , 1.44930e-03 , 0.0 , 5.68688} , { -1.27599e-02 , 1.37523e-01 , -3.44225e-02 , 2.40203e-03 , 0.0 , 5.6862 } , { -1.66651e-02 , 1.85760e-01 , -4.85808e-02 , 3.52947e-03 , 0.0 ,  5.68643 } , {  -2.05442e-02 , 2.19814e-01 , -6.16465e-02 , 4.91569e-03 , 0.0 ,  5.68753 }  , { -2.08654e-02 , 2.34196e-01 , -6.95590e-02 , 5.88573e-03 , 0.0 , 5.68962 } }; 
        const Double_t v2_K_params[6][6] = { { -2.28987e-02 , 6.91034e-02 , -1.48611e-02  ,  5.67949e-04 , 0.0 , 3.79629 } , { -3.28047e-02 ,  1.02630e-01  , -2.03472e-02 , 7.34138e-04 , 0.0 , 3.79652 } , { -4.44470e-02 , 1.43257e-01 , -2.79470e-02  , 8.88236e-04 , 0.0 , 3.79639 } , { -5.50614e-02 , 1.89560e-01 , -3.80298e-02 , 1.13837e-03 , 0.0 , 3.79629 } , { -6.46262e-02 ,  2.32903e-01 , -5.62942e-02 , 3.39183e-03 , 0.0 ,  3.79624 }  , { -6.68228e-02 , 2.58028e-01 , -7.21831e-02 ,  5.84160e-03 , 0.0 , 3.79576 } };
        const Double_t v2_P_params[6][6] = { {1.34354e-03 ,   -2.69847e-02 , 4.68928e-02 , -1.30276e-02 , 1.03297e-03 , 5.74363 } , { 1.68188e-02 , -5.83743e-02 , 7.93535e-02 , -2.04909e-02 , 1.53971e-03  , 5.74345 } , { 2.33248e-02 , -7.58826e-02 , 1.08507e-01 , -2.81839e-02 , 2.11781e-03 , 5.74522 } , { 1.17865e-02 , -5.86004e-02 , 1.24604e-01 , -3.45561e-02 , 2.68779e-03 ,  5.74372 } , { -1.40410e-03 , -3.01182e-02 , 1.26896e-01 , -3.80857e-02 , 3.12556e-03 ,  5.74289 }  , { -2.61424e-02 , 3.52815e-02 , 9.56710e-02  , -3.30046e-02 , 2.85671e-03 ,  5.74307 } };
        const Double_t v3_pi_params[6][6] = { { -2.53410e-02 , 7.65854e-02  , -1.64087e-02 , 9.27325e-04 , 0.0 , 5.6864  } , { -2.47587e-02  ,  8.05937e-02 , -1.60442e-02 , 7.45643e-04 , 0.0 , 5.68688 } , {  -2.40356e-02 , 8.39185e-02 ,  -1.60066e-02 , 6.23750e-04 , 0.0, 5.6862 } , { -2.75520e-02  ,  9.74475e-02 , -2.11485e-02 , 1.13159e-03 , 0.0 , 5.68643 } , { -2.89726e-02  , 1.05326e-01 , -2.52450e-02 ,  1.65368e-03 , 0.0 , 5.68753 }  , { -2.29852e-02  , 9.41400e-02  , -1.86962e-02 , -5.14131e-04  , 2.43303e-04 , 5.68962 } };
        const Double_t v3_K_params[6][6] = { { -2.00467e-02 , 4.04763e-02 , 4.87640e-03 , -2.25452e-03 , 0.0 , 3.79629 } , {  -2.12149e-02 , 4.44332e-02 , 5.71953e-03 , -2.55308e-03 , 0.0 , 3.79652  } , { -2.35466e-02 , 5.22392e-02 , 4.10174e-03  ,  -2.57112e-03 , 0.0 , 3.79639 } , { -2.76820e-02 , 6.67519e-02 , -1.63168e-03  , -2.03702e-03 , 0.0 , 3.79629  } , { -3.02956e-02 , 7.58570e-02  , -5.81264e-03 , -1.57544e-03  , 0.0 , 3.79624 }  , { -2.46987e-02 , 6.96050e-02 , -4.33787e-03 , -1.83356e-03 , 0.0 ,  3.79576 } };
        const Double_t v3_P_params[6][6] = { { -4.69818e-03 , -1.24501e-02 , 2.67730e-02 ,  -4.02833e-03 , 0.0 , 5.74363 } , { -9.67248e-03  , -4.64947e-03 , 2.65129e-02 , -4.13887e-03 , 0.0 , 5.74345 } , { -8.40764e-03 , -9.29454e-03  , 3.21732e-02 ,  -5.16580e-03 , 0.0 , 5.74522 } , { -2.37472e-02 , 1.73042e-02 , 2.41911e-02 , -4.40537e-03 , 0.0 , 5.74372 } , { -2.26647e-02 , 2.18590e-02 , 2.34411e-02 , -4.53662e-03  , 0.0 , 5.74289 }  , { -4.06477e-02 , 5.58970e-02 ,  9.23320e-03 , -2.90313e-03 , 0.0 ,  5.74307 } };
        const Double_t v4_pi_params[6][6] = { { -2.61775e-02 , 5.41027e-02 , -8.11629e-03 , 1.07636e-04 , 0.0 , 5.56759 } , { -2.84327e-02  , 6.15741e-02 , -1.08541e-02 , 3.74468e-04 , 0.0 , 5.56826 } , { -2.50069e-02 ,  5.64826e-02  , -7.82525e-03 , -5.60331e-05 , 0.0 , 5.56784 } , { -3.66736e-02 , 8.52490e-02 , -2.18799e-02 , 1.76828e-03 , 0.0 ,  5.56903 } , { -2.67396e-02 , 6.76440e-02 , -1.25867e-02 , 3.54348e-04 , 0.0 , 5.57063 }  , { -2.83923e-02 ,  7.37708e-02 , -1.73843e-02  , 1.07088e-03 , 0.0 , 5.57323 } };
        const Double_t v4_K_params[6][6] = { { -1.96656e-02 , 2.61057e-02 , 7.38473e-03  , -2.14579e-03 , 0.0 , 3.88195 } , { -2.09796e-02  , 2.93941e-02 ,  6.44392e-03 , -2.01671e-03 , 0.0 , 3.88229 } , { -2.32744e-02 , 3.62767e-02  ,  3.93014e-03  , -1.89347e-03 , 0.0 , 3.88197 } , { -2.02513e-02 , 2.93811e-02 ,  1.03741e-02 , -3.26730e-03 , 0.0 , 3.8817 } , { -1.45235e-02 , 1.84089e-02 , 1.63240e-02 , -4.16862e-03 , 0.0 , 3.88098 }  , { -1.96178e-02 , 3.49596e-02 , 5.84558e-03 , -2.54327e-03 , 0.0 , 3.8814 } };
        const Double_t v4_P_params[6][6] = { {  1.68565e-02 , -4.84043e-02 , 3.57252e-02 , -4.69793e-03 , 0.0 , 5.62931 } , { 1.84244e-02  , -5.34851e-02 , 4.00556e-02 , -5.30060e-03 , 0.0 , 5.62921 } , { 1.53098e-02 , -4.93400e-02 , 3.99301e-02 , -5.45571e-03 , 0.0 , 5.63246 } , { 1.31460e-02 , -4.65074e-02 , 4.21324e-02 , -5.94997e-03 , 0.0 ,  5.6322 } , { 9.26040e-03 , -4.29033e-02 , 4.25759e-02 , -6.31299e-03 , 0.0 , 5.6311 }  , { -2.62813e-02 , 2.86358e-02  , 8.47114e-03 , -1.90998e-03 , 0.0 , 5.63034 } };
        const Double_t v5_pi_params[6][6] = { {  -2.27633e-02 , 3.62812e-02 , -5.09327e-03 , 3.04512e-06 , 0.0 , 5.3046 } , { -2.11594e-02 , 3.40405e-02 , -4.51201e-03 ,  1.65056e-04 , 0.0 , 5.30516 } , { -1.79688e-02 , 2.73746e-02 , 1.11420e-03  , -1.03606e-03 , 0.0 , 5.30769 } , { -1.64322e-02 , 2.88779e-02 , 3.15650e-05 , -6.48281e-04  , 0.0 , 5.31082 } , { -2.00310e-02 , 3.98351e-02 , -7.27710e-03  , 5.18693e-04, 0.0 ,  5.31308 }  , {  -2.14945e-02 , 4.30007e-02  , -8.88504e-03 , 4.45375e-04 , 0.0 , 5.31678 } };
        const Double_t v5_K_params[6][6] = { { -2.93359e-02 , 4.03216e-02 , -6.72350e-03 , 2.49656e-04 , 0.0 , 3.68094 } , { -2.45592e-02 ,  3.08562e-02 , -1.20602e-03 , -5.88864e-04 , 0.0 , 3.68144 } , { -2.23697e-02 , 2.31996e-02 , 5.23514e-03, -1.79606e-03 , 0.0 , 3.68252 } , { -2.38504e-02 , 3.07280e-02 , 3.01792e-04 , -1.04716e-03 , 0.0 , 3.68158 } , { -3.33063e-02 , 5.52523e-02 , -1.62386e-02 , 1.98994e-03 , 0.0 , 3.68023 }  , { -3.34578e-02 , 5.08559e-02 , -1.28618e-02 , 1.32781e-03 , 0.0 , 3.68031 } };
        const Double_t v5_P_params[6][6] = { { 3.20157e-02 ,  -8.18906e-02 , 5.85698e-02 ,  -1.34328e-02  ,  1.01128e-03  , 5.37351 } , { 1.08684e-02 , -3.07064e-02 , 2.07656e-02 , -2.71231e-03  , 4.65723e-05 , 5.37475 } , { 3.93613e-02 , -1.01226e-01 , 7.84549e-02 , -1.98730e-02 , 1.67558e-03  , 5.37907 } , { 5.87097e-02 ,  -1.42471e-01 , 1.07133e-01 , -2.72022e-02 ,  2.30510e-03 , 5.3852 } , {  5.54421e-02 , -1.37351e-01 , 1.02622e-01 , -2.52028e-02  , 1.99010e-03 ,  5.3828 }  , { 7.55687e-02 , -2.06518e-01 , 1.57785e-01 , -4.01594e-02 , 3.28859e-03 ,  5.38206 } };
        const Double_t piPlus_params[8][5] = { { 0.139570 , 0.947908 , 0.0737422 , 0.968789 , 2207060 } , { 0.139570 , 0.943360 , 0.078973 , 1.02539 , 1500940 } , { 0.139570 , 0.95381 , 0.0708164 , 0.987509 , 1603690 } , { 0.139570 , 0.95409 , 0.0673556 , 1.00806 , 1319880 } , { 0.139570 , 0.959811 , 0.0676647 , 1.08715 , 861341 } , { 0.139570 , 0.957266 , 0.08703763 , 1.21600, 468819 } , { 0.139570 , 0.949426 , 0.078257 , 1.46891 , 190742 } , { 0.139570 , 0.954145 , 0.0744700 , 1.58732, 119812 }};
        const Double_t piMinus_params[8][5] = { { 0.139570 , 0.942100 , 0.072976 , 0.967010 , 1856920 } , { 0.139570 , 0.944275 , 0.0766907 , 0.986230 , 1583560 } , { 0.139570 , 0.945046 , 0.072708 , 1.01734 , 1160050 } , { 0.139570 , 0.952089 , 0.0729955 , 1.03405 , 970911 } , { 0.139570 , 0.952325 , 0.0737461 , 1.12383 , 620421 } , { 0.139570 , 0.940441 , 0.0844949 , 1.39199, 240963 } , { 0.139570 , 0.933892 , 0.0902015 , 1.63995 , 112371 } , { 0.139570 , 0.927508 , 0.0966857 , 2.00, 47185.9 }};
        const Double_t kPlus_params[8][5] = { { 0.493677 , 0.826760 , 0.129893 , 0.703341 , 230444 } , { 0.493677 , 0.820456 , 0.135416 , 0.754479 , 154330 } , { 0.493677 , 0.822368 , 0.136030 , 0.790558 , 113729 } , { 0.493677 , 0.812882 , 0.14257 , 0.868716 , 60544.6 } , { 0.493677 , 0.808265 , 0.148363 , 1.01543 , 32570.8 } , { 0.493677 , 0.741293 , 0.185352 , 1.36225, 7164.93 } , { 0.493677 , 0.771795 , 0.180722, 1.77221 , 4536.61 } , { 0.493677 , 0.823896 , 0.156262 , 1.94719, 4555.18 }};
        const Double_t kMinus_params[8][5] = { { 0.493677 , 0.821301 , 0.133049 , 0.748879 , 203970 } , { 0.493677 , 0.821564 , 0.134495 , 0.781053 , 160695 } , { 0.493677 , 0.816101 , 0.138573 , 0.816978 , 103440 } , { 0.493677 , 0.807483 , 0.146798 , 0.936718 , 52652.1 } , { 0.493677 , 0.802314 , 0.152261 , 1.08377 , 28720.5 } , { 0.493677 , 0.777818 , 0.171125 , 1.41012, 10213.3 } , { 0.493677 , 0.805967 , 0.162440 , 1.76754 , 7498.70 } , { 0.493677 , 0.830498 , 0.151675 , 1.97103, 5279.26 }};
        const Double_t p_params[8][5] = { { 0.938272 , 0.861255 , 0.113053 , 0.646603 , 3540660 } , { 0.938272 , 0.829143 , 0.134133 , 0.583048 , 608721 } , { 0.93827 , 0.784519 , 0.158416 , 0.489724 , 121635 } , { 0.938272 , 0.755180 , 0.170715 , 0.472479 , 49103.6 } , { 0.938272 , 0.738171 , 0.178862 , 0.549364 , 24070.7 } , { 0.938272 , 0.718834 , 0.188403 , 0.698602 , 11050.3 } , { 0.938272 , 0.685109 , 0.204819 , 0.914109 , 3943.15 } , { 0.938272 , 0.638029 , 0.244535 , 1.78914 , 768.023 }};
        const Double_t pbar_params[8][5] = { { 0.938272 , 0.863601 , 0.11476 , 0.658438 , 4073300 } , { 0.938272 , 0.855439 , 0.116446 , 0.647756 , 2212850 } , { 0.938272 , 0.832520 , 0.131665 , 0.622898 , 552092 } , { 0.938272 , 0.778317 , 0.160757 , 0.539450 , 77034.0 } , { 0.938272 , 0.739188 , 0.174989 , 0.532042 , 28702.3 } , { 0.938272 , 0.718917 , 0.186137 , 0.680732 , 12351.2 } , { 0.938272 , 0.688116 , 0.202278 , 0.939983 , 4397.05 } , { 0.938272 , 0.637255 , 0.244434 , 1.83999 , 799.620 }};
        const Double_t piZero_params[8][5] = { { 0.139977 , 0.942100 , 0.072976 , 0.967010 , 1856920 } , { 0.139977 , 0.944275 , 0.0766907 , 0.986230 , 1583560 } , { 0.139977, 0.945046 , 0.072708 , 1.01734 , 1160050 } , { 0.139977 , 0.952089 , 0.0729955 , 1.03405 , 970911 } , { 0.139977 , 0.952325 , 0.0737461 , 1.12383 , 620421 } , { 0.139977 , 0.940441 , 0.0844949 , 1.39199, 240963 } , { 0.139977 , 0.933892 , 0.0902015 , 1.63995 , 112371 } , { 0.139977 , 0.927508 , 0.0966857 , 2.00, 47185.9 }};        
        const Double_t v1_pi_params[6][6] = { {-6.29928e-03 , 6.29213e-02 , -1.69603e-02 , 1.24408e-03 , 0.0 , 5.6864 } , { -8.15353e-03 ,  9.46349e-02 , -2.27846e-02 , 1.44930e-03 , 0.0 , 5.68688} , { -1.27599e-02 , 1.37523e-01 , -3.44225e-02 , 2.40203e-03 , 0.0 , 5.6862 } , { -1.66651e-02 , 1.85760e-01 , -4.85808e-02 , 3.52947e-03 , 0.0 ,  5.68643 } , {  -2.05442e-02 , 2.19814e-01 , -6.16465e-02 , 4.91569e-03 , 0.0 ,  5.68753 }  , { -2.08654e-02 , 2.34196e-01 , -6.95590e-02 , 5.88573e-03 , 0.0 , 5.68962 } };
        const Double_t v1_K_params[6][6] = {{ -2.28987e-02 , 6.91034e-02 , -1.48611e-02  ,  5.67949e-04 , 0.0 , 3.79629 } , { -3.28047e-02 ,  1.02630e-01  , -2.03472e-02 , 7.34138e-04 , 0.0 , 3.79652 } , { -4.44470e-02 , 1.43257e-01 , -2.79470e-02  , 8.88236e-04 , 0.0 , 3.79639 } , { -5.50614e-02 , 1.89560e-01 , -3.80298e-02 , 1.13837e-03 , 0.0 , 3.79629 } , { -6.46262e-02 ,  2.32903e-01 , -5.62942e-02 , 3.39183e-03 , 0.0 ,  3.79624 }  , { -6.68228e-02 , 2.58028e-01 , -7.21831e-02 ,  5.84160e-03 , 0.0 , 3.79576 } };
        const Double_t v1_P_params[6][6] = {{1.34354e-03 ,   -2.69847e-02 , 4.68928e-02 , -1.30276e-02 , 1.03297e-03 , 5.74363 } , { 1.68188e-02 , -5.83743e-02 , 7.93535e-02 , -2.04909e-02 , 1.53971e-03  , 5.74345 } , { 2.33248e-02 , -7.58826e-02 , 1.08507e-01 , -2.81839e-02 , 2.11781e-03 , 5.74522 } , { 1.17865e-02 , -5.86004e-02 , 1.24604e-01 , -3.45561e-02 , 2.68779e-03 ,  5.74372 } , { -1.40410e-03 , -3.01182e-02 , 1.26896e-01 , -3.80857e-02 , 3.12556e-03 ,  5.74289 }  , { -2.61424e-02 , 3.52815e-02 , 9.56710e-02  , -3.30046e-02 , 2.85671e-03 ,  5.74307 } };

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
        pTdist_piZero = new TF1("pTdist_piZero",BlastWavedNdptTimesPt ,0,100,5); 
        pTdist_piZero->SetNpx (1000);  

        v1_pi = new TF1("v1_pi", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4) - 0.02) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4) - 0.02)", 0, 100 );
        v1_pi->SetNpx(1000);
        v2_pi = new TF1("v2_pi", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4) ) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
        v2_pi->SetNpx(1000);
        v3_pi = new TF1("v3_pi", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
        v3_pi->SetNpx(1000);
        v4_pi = new TF1("v4_pi", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
        v4_pi->SetNpx(1000);
        v5_pi = new TF1("v5_pi", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
        v5_pi->SetNpx(1000);
        v1_K = new TF1("v1_K", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4) - 0.02) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4) - 0.02)", 0, 100 );
        v1_K->SetNpx(1000);
        v2_K = new TF1("v2_K", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
        v2_K->SetNpx(1000);
        v3_K = new TF1("v3_K", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
        v3_K->SetNpx(1000);
        v4_K = new TF1("v4_K", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
        v4_K->SetNpx(1000);
        v5_K = new TF1("v5_K", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
        v5_K->SetNpx(1000);
        v1_P = new TF1("v1_P", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4) - 0.02) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4) - 0.02)", 0, 100 );
        v1_P->SetNpx(1000);
        v2_P = new TF1("v2_P", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
        v2_P->SetNpx(1000);
        v3_P = new TF1("v3_P", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
        v3_P->SetNpx(1000);
        v4_P = new TF1("v4_P", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
        v4_P->SetNpx(1000);
        v5_P = new TF1("v5_P", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
        v5_P->SetNpx(1000);

        v1_pi->SetParameter(0,v1_pi_params[fCentBin][0]);
        v1_pi->SetParameter(1,v1_pi_params[fCentBin][1]);
        v1_pi->SetParameter(2,v1_pi_params[fCentBin][2]);
        v1_pi->SetParameter(3,v1_pi_params[fCentBin][3]);
        v1_pi->SetParameter(4,v1_pi_params[fCentBin][4]);
        v1_pi->SetParameter(5,v1_pi_params[fCentBin][5]);
        
        v1_K->SetParameter(0,v1_K_params[fCentBin][0]);
        v1_K->SetParameter(1,v1_K_params[fCentBin][1]);
        v1_K->SetParameter(2,v1_K_params[fCentBin][2]);
        v1_K->SetParameter(3,v1_K_params[fCentBin][3]);
        v1_K->SetParameter(4,v1_K_params[fCentBin][4]);
        v1_K->SetParameter(5,v1_K_params[fCentBin][5]);
        
        v1_P->SetParameter(0,v1_P_params[fCentBin][0]);
        v1_P->SetParameter(1,v1_P_params[fCentBin][1]);
        v1_P->SetParameter(2,v1_P_params[fCentBin][2]);
        v1_P->SetParameter(3,v1_P_params[fCentBin][3]);
        v1_P->SetParameter(4,v1_P_params[fCentBin][4]);
        v1_P->SetParameter(5,v1_P_params[fCentBin][5]);
        
        v2_pi->SetParameter(0,v2_pi_params[fCentBin][0]);
        v2_pi->SetParameter(1,v2_pi_params[fCentBin][1]);
        v2_pi->SetParameter(2,v2_pi_params[fCentBin][2]);
        v2_pi->SetParameter(3,v2_pi_params[fCentBin][3]);
        v2_pi->SetParameter(4,v2_pi_params[fCentBin][4]);
        v2_pi->SetParameter(5,v2_pi_params[fCentBin][5]);
        
        v2_K->SetParameter(0,v2_K_params[fCentBin][0]);
        v2_K->SetParameter(1,v2_K_params[fCentBin][1]);
        v2_K->SetParameter(2,v2_K_params[fCentBin][2]);
        v2_K->SetParameter(3,v2_K_params[fCentBin][3]);
        v2_K->SetParameter(4,v2_K_params[fCentBin][4]);
        v2_K->SetParameter(5,v2_K_params[fCentBin][5]);
        
        v2_P->SetParameter(0,v2_P_params[fCentBin][0]);
        v2_P->SetParameter(1,v2_P_params[fCentBin][1]);
        v2_P->SetParameter(2,v2_P_params[fCentBin][2]);
        v2_P->SetParameter(3,v2_P_params[fCentBin][3]);
        v2_P->SetParameter(4,v2_P_params[fCentBin][4]);
        v2_P->SetParameter(5,v2_P_params[fCentBin][5]);

    
        v3_pi->SetParameter(0,v3_pi_params[fCentBin][0]);
        v3_pi->SetParameter(1,v3_pi_params[fCentBin][1]);
        v3_pi->SetParameter(2,v3_pi_params[fCentBin][2]);
        v3_pi->SetParameter(3,v3_pi_params[fCentBin][3]);
        v3_pi->SetParameter(4,v3_pi_params[fCentBin][4]);
        v3_pi->SetParameter(5,v3_pi_params[fCentBin][5]);
        
        v3_K->SetParameter(0,v3_K_params[fCentBin][0]);
        v3_K->SetParameter(1,v3_K_params[fCentBin][1]);
        v3_K->SetParameter(2,v3_K_params[fCentBin][2]);
        v3_K->SetParameter(3,v3_K_params[fCentBin][3]);
        v3_K->SetParameter(4,v3_K_params[fCentBin][4]);
        v3_K->SetParameter(5,v3_K_params[fCentBin][5]);
        
        v3_P->SetParameter(0,v3_P_params[fCentBin][0]);
        v3_P->SetParameter(1,v3_P_params[fCentBin][1]);
        v3_P->SetParameter(2,v3_P_params[fCentBin][2]);
        v3_P->SetParameter(3,v3_P_params[fCentBin][3]);
        v3_P->SetParameter(4,v3_P_params[fCentBin][4]);
        v3_P->SetParameter(5,v3_P_params[fCentBin][5]);

        
        v4_pi->SetParameter(0,v4_pi_params[fCentBin][0]);
        v4_pi->SetParameter(1,v4_pi_params[fCentBin][1]);
        v4_pi->SetParameter(2,v4_pi_params[fCentBin][2]);
        v4_pi->SetParameter(3,v4_pi_params[fCentBin][3]);
        v4_pi->SetParameter(4,v4_pi_params[fCentBin][4]);
        v4_pi->SetParameter(5,v4_pi_params[fCentBin][5]);
        
        v4_K->SetParameter(0,v4_K_params[fCentBin][0]);
        v4_K->SetParameter(1,v4_K_params[fCentBin][1]);
        v4_K->SetParameter(2,v4_K_params[fCentBin][2]);
        v4_K->SetParameter(3,v4_K_params[fCentBin][3]);
        v4_K->SetParameter(4,v4_K_params[fCentBin][4]);
        v4_K->SetParameter(5,v4_K_params[fCentBin][5]);
        
        v4_P->SetParameter(0,v4_P_params[fCentBin][0]);
        v4_P->SetParameter(1,v4_P_params[fCentBin][1]);
        v4_P->SetParameter(2,v4_P_params[fCentBin][2]);
        v4_P->SetParameter(3,v4_P_params[fCentBin][3]);
        v4_P->SetParameter(4,v4_P_params[fCentBin][4]);
        v4_P->SetParameter(5,v4_P_params[fCentBin][5]);

        v5_pi->SetParameter(0,v5_pi_params[fCentBin][0]);
        v5_pi->SetParameter(1,v5_pi_params[fCentBin][1]);
        v5_pi->SetParameter(2,v5_pi_params[fCentBin][2]);
        v5_pi->SetParameter(3,v5_pi_params[fCentBin][3]);
        v5_pi->SetParameter(4,v5_pi_params[fCentBin][4]);
        v5_pi->SetParameter(5,v5_pi_params[fCentBin][5]);
        
        v5_K->SetParameter(0,v5_K_params[fCentBin][0]);
        v5_K->SetParameter(1,v5_K_params[fCentBin][1]);
        v5_K->SetParameter(2,v5_K_params[fCentBin][2]);
        v5_K->SetParameter(3,v5_K_params[fCentBin][3]);
        v5_K->SetParameter(4,v5_K_params[fCentBin][4]);
        v5_K->SetParameter(5,v5_K_params[fCentBin][5]);
        
        v5_P->SetParameter(0,v5_P_params[fCentBin][0]);
        v5_P->SetParameter(1,v5_P_params[fCentBin][1]);
        v5_P->SetParameter(2,v5_P_params[fCentBin][2]);
        v5_P->SetParameter(3,v5_P_params[fCentBin][3]);
        v5_P->SetParameter(4,v5_P_params[fCentBin][4]);
        v5_P->SetParameter(5,v5_P_params[fCentBin][5]);


        pTdist_piPlus->SetParameter(0,piPlus_params[fCentBin][0]);
        pTdist_piPlus->SetParameter(1,piPlus_params[fCentBin][1]);
        pTdist_piPlus->SetParameter(2,piPlus_params[fCentBin][2]);
        pTdist_piPlus->SetParameter(3,piPlus_params[fCentBin][3]);
        pTdist_piPlus->SetParameter(4,piPlus_params[fCentBin][4]);

    
        pTdist_piMinus->SetParameter(0,piMinus_params[fCentBin][0]);
        pTdist_piMinus->SetParameter(1,piMinus_params[fCentBin][1]);
        pTdist_piMinus->SetParameter(2,piMinus_params[fCentBin][2]);
        pTdist_piMinus->SetParameter(3,piMinus_params[fCentBin][3]);
        pTdist_piMinus->SetParameter(4,piMinus_params[fCentBin][4]);

        pTdist_kPlus->SetParameter(0,kPlus_params[fCentBin][0]);
        pTdist_kPlus->SetParameter(1,kPlus_params[fCentBin][1]);
        pTdist_kPlus->SetParameter(2,kPlus_params[fCentBin][2]);
        pTdist_kPlus->SetParameter(3,kPlus_params[fCentBin][3]);
        pTdist_kPlus->SetParameter(4,kPlus_params[fCentBin][4]);

    
        pTdist_kMinus->SetParameter(0,kMinus_params[fCentBin][0]);
        pTdist_kMinus->SetParameter(1,kMinus_params[fCentBin][1]);
        pTdist_kMinus->SetParameter(2,kMinus_params[fCentBin][2]);
        pTdist_kMinus->SetParameter(3,kMinus_params[fCentBin][3]);
        pTdist_kMinus->SetParameter(4,kMinus_params[fCentBin][4]);

        pTdist_p->SetParameter(0,p_params[fCentBin][0]);
        pTdist_p->SetParameter(1,p_params[fCentBin][1]);
        pTdist_p->SetParameter(2,p_params[fCentBin][2]);
        pTdist_p->SetParameter(3,p_params[fCentBin][3]);
        pTdist_p->SetParameter(4,p_params[fCentBin][4]);


        pTdist_pbar->SetParameter(0,pbar_params[fCentBin][0]);
        pTdist_pbar->SetParameter(1,pbar_params[fCentBin][1]);
        pTdist_pbar->SetParameter(2,pbar_params[fCentBin][2]);
        pTdist_pbar->SetParameter(3,pbar_params[fCentBin][3]);
        pTdist_pbar->SetParameter(4,pbar_params[fCentBin][4]);


        pTdist_piZero->SetParameter(0,piZero_params[fCentBin][0]);
        pTdist_piZero->SetParameter(1,piZero_params[fCentBin][1]);
        pTdist_piZero->SetParameter(2,piZero_params[fCentBin][2]);
        pTdist_piZero->SetParameter(3,piZero_params[fCentBin][3]);
        pTdist_piZero->SetParameter(4,piZero_params[fCentBin][4]);

        timestamp = new TTimeStamp();
        UInt_t Seed = UInt_t(timestamp->GetSec());
        UInt_t seed_eta = Seed * 7;
        etadist = new TRandom3(seed_eta);
        UInt_t seed_psi_1 = (Seed*2)+3;
        Psi_1 = new TRandom3(seed_psi_1);
        UInt_t seed_psi_2 = (Seed*2)+1;
        Psi_2 = new TRandom3(seed_psi_2);
        UInt_t seed_psi_3 = (Seed*3)+2;
        Psi_3 = new TRandom3(seed_psi_3);
        UInt_t seed_psi_4 = (Seed*3)+4;
        Psi_4 = new TRandom3(seed_psi_4);
        UInt_t seed_psi_5 = (Seed*3)+5;
        Psi_5 = new TRandom3(seed_psi_5);
        UInt_t seed_uniform_random_phi = (Seed*3)*4;
        Harmonics_Phi_Dist_Rand = new TRandom3(seed_uniform_random_phi);
        UInt_t seed_pT = (Seed*2)*3;
        pT_Rand = new TRandom3(seed_pT);  

        doVn[0] = 1;
        doVn[1] = 1;
        doVn[2] = 1;
        doVn[3] = 1;
        doVn[4] = 1;

        fixPsi[0] = 0;
        fixPsi[1] = 1;
        fixPsi[2] = 0;
        fixPsi[3] = 1;
        fixPsi[4] = 0;
        
        fixed_psi[1] = 0.0;
        fixed_psi[3] = 0.0;
    }

    // setters
    void SetDoVn(Int_t n, Int_t dovn) {doVn[n-1] = dovn;};
    
    void SetFixedPsi(Int_t n, Double_t fixedpsi) {
        fixed_psi[n-1] = fixedpsi;
        fixPsi[n-1] = 1;
    };

    void doQA(){
        cout << "TODO: QA" << endl;
    }

    TGEvent Next();

    TGEvent NextPbPb();

    TGEvent NextAuAu();

};

TennGen::TennGen(Int_t collen, Int_t centbin, Double_t etarange) 
{
    fCollEn = collen;
    fCentBin = centbin;
    fEtaRange = etarange;
    
    CheckInputs();

    if(fCollEn == 200){
        InitAuAu();
    }
    else if(fCollEn == 2760){
        InitPbPb();
    }
    
}

TennGen::~TennGen()
{
    // delete pTdist_piPlus;
    // delete pTdist_piMinus;
    // delete pTdist_kPlus;
    // delete pTdist_kMinus;
    // delete pTdist_p;
    // delete pTdist_pbar;
    // delete etadist;
    // delete Psi_1;
    // delete Psi_2;
    // delete Psi_3;
    // delete Psi_4;
    // delete Harmonics_Phi_Dist_Rand;
    // delete Multi_Rand;
    // delete pT_Rand;
    // delete timestamp;
    // delete fout;

    // delete v2_pi;
    // delete v2_k;
    // delete v2_p;
    // delete v3_pi;
    // delete v3_k;
    // delete v3_p;
    // delete v4_pi;
    // delete v4_k;
    // delete v4_p;

}

TGEvent TennGen::Next(){

    // check collEn
    if(fCollEn == 200){
        return NextAuAu();
    }
    else if(fCollEn == 2760){
        return NextPbPb();
    }
}

TGEvent TennGen::NextPbPb(){

    TGEvent event;
    event.particles.clear();

    if(fixPsi[0]) event.psi1 = fixed_psi[0];
    else event.psi1 = Psi_1->Uniform( 0 , 2.0*TMath::Pi());
    if(fixPsi[1]) event.psi2 = fixed_psi[1];
    else event.psi2 = 0.0;
    if(fixPsi[2]) event.psi3 = fixed_psi[2];
    else event.psi3 = Psi_3->Uniform( 0 , 2.0*TMath::Pi());
    if(fixPsi[3]) event.psi4 = fixed_psi[3];
    else event.psi4 = 0.0;
    if(fixPsi[4]) event.psi5 = fixed_psi[4];
    else event.psi5 = Psi_5->Uniform( 0 , 2.0*TMath::Pi());

    
    event.nparts = -1;
    event.centBin = fCentBin;
    event.etaRange = fEtaRange;
    event.collEn = fCollEn;

    const Double_t yield_arr[8][7] = { { (TMath::Ceil((fEtaRange * 654)/0.5)) , (TMath::Ceil((fEtaRange * 654)/0.5)) , (TMath::Ceil((fEtaRange * 97)/0.5)) , (TMath::Ceil((fEtaRange * 97)/0.5)) , (TMath::Ceil((fEtaRange * 30)/0.5)) , (TMath::Ceil((fEtaRange * 29)/0.5)) , (TMath::Ceil((fEtaRange * 654)/0.5)) } , { (TMath::Ceil((fEtaRange * 541)/0.5)) , (TMath::Ceil((fEtaRange * 539)/0.5)) , (TMath::Ceil((fEtaRange * 81)/0.5)) , (TMath::Ceil((fEtaRange * 80)/0.5)) , (TMath::Ceil((fEtaRange * 25)/0.5)) , (TMath::Ceil((fEtaRange * 25)/0.5)) , (TMath::Ceil((fEtaRange * 539)/0.5)) } , { (TMath::Ceil((fEtaRange * 406)/0.5)) , (TMath::Ceil((fEtaRange * 404)/0.5)) , (TMath::Ceil((fEtaRange * 61)/0.5)) , (TMath::Ceil((fEtaRange * 61)/0.5)) , (TMath::Ceil((fEtaRange * 19)/0.5)) , (TMath::Ceil((fEtaRange * 19)/0.5)) , (TMath::Ceil((fEtaRange * 404)/0.5)) } , { (TMath::Ceil((fEtaRange * 274)/0.5)) , (TMath::Ceil((fEtaRange * 273 )/0.5)), (TMath::Ceil((fEtaRange * 41)/0.5)) , (TMath::Ceil((fEtaRange * 41)/0.5)) , (TMath::Ceil((fEtaRange * 13)/0.5)) , (TMath::Ceil((fEtaRange * 13)/0.5)) , (TMath::Ceil((fEtaRange * 273)/0.5)) } , { (TMath::Ceil((fEtaRange * 179)/0.5)) , (TMath::Ceil((fEtaRange * 179)/0.5)) , (TMath::Ceil((fEtaRange * 27)/0.5)) , (TMath::Ceil((fEtaRange * 27)/0.5)) , (TMath::Ceil((fEtaRange * 9)/0.5)) , (TMath::Ceil((fEtaRange * 9)/0.5)) , (TMath::Ceil((fEtaRange * 179)/0.5)) } , { (TMath::Ceil((fEtaRange * 111)/0.5)) , (TMath::Ceil((fEtaRange * 110)/0.5)) , (TMath::Ceil((fEtaRange * 16)/0.5)) , (TMath::Ceil((fEtaRange * 16)/0.5)) , (TMath::Ceil((fEtaRange * 5)/0.5)) , (TMath::Ceil((fEtaRange * 5)/0.5)) , 	(TMath::Ceil((fEtaRange * 110)/0.5)) } , { (TMath::Ceil((fEtaRange * 63)/0.5)) , (TMath::Ceil((fEtaRange * 63)/0.5)) , (TMath::Ceil((fEtaRange * 9)/0.5)) , (TMath::Ceil((fEtaRange * 9)/0.5)) , (TMath::Ceil((fEtaRange * 4)/0.5)) , (TMath::Ceil((fEtaRange * 4)/0.5)) , (TMath::Ceil((fEtaRange * 63)/0.5)) } , { (TMath::Ceil((fEtaRange * 33)/0.5)) , (TMath::Ceil((fEtaRange * 33)/0.5)) , ((TMath::Ceil(fEtaRange * 4)/0.5)) , (TMath::Ceil((fEtaRange * 4)/0.5)) , (TMath::Ceil((fEtaRange * 2)/0.5)) , (TMath::Ceil((fEtaRange * 2)/0.5)) , (TMath::Ceil((fEtaRange * 33)/0.5)) } }; //setting multiplicities for particles

    for(Int_t n = 0; n<7; n++)
    {
        for(Int_t ipart = 0; ipart< yield_arr[fCentBin][n]; ipart++)
        {
            TGParticle particle;    
            if (n==0){
                particle.pdg = -211;
                particle.mass = 0.139570;
                particle.charge = -1;
                particle.pT = pTdist_piMinus->GetRandom(pT_Rand);
            }
            if (n==1){
                particle.pdg = 211;
                particle.mass = 0.139570;
                particle.charge = 1;
                particle.pT = pTdist_piPlus->GetRandom(pT_Rand);
            } 
            if (n==2){
                particle.pdg = -321;
                particle.mass = 0.493677;
                particle.charge = -1;
                particle.pT = pTdist_kMinus->GetRandom(pT_Rand);
            }
            if (n==3){
                particle.pdg = 321;
                particle.mass = 0.493677;
                particle.charge = 1;
                particle.pT = pTdist_kPlus->GetRandom(pT_Rand);
            }
            if (n==4){
                particle.pdg = -2212;
                particle.mass = 0.938272;
                particle.charge = -1;
                particle.pT = pTdist_pbar->GetRandom(pT_Rand);
            }
            if (n==5){
                particle.pdg = 2212;
                particle.mass = 0.938272;
                particle.charge = 1;
                particle.pT = pTdist_p->GetRandom(pT_Rand);
            }
            //pi0
            if(n==6){
                particle.pdg = 111;
                particle.mass = 0.1349766;
                particle.charge = 0;
                particle.pT = pTdist_piZero->GetRandom(pT_Rand);
            }
            
            particle.eta = etadist->Uniform(-1.0*event.etaRange , 1.0*event.etaRange);

            Double_t v1_eval, v2_eval, v3_eval, v4_eval, v5_eval;
            if(n==0||n==1||n==6){ // get pion vN
                v1_eval = v1_pi->Eval( particle.pT , 0.0 , 0.0 ,  0.0 );
                v2_eval = v2_pi->Eval( particle.pT , 0.0 , 0.0 ,  0.0 );
                v3_eval = v3_pi->Eval( particle.pT , 0.0 , 0.0 ,  0.0 );
                v4_eval = v4_pi->Eval( particle.pT , 0.0 , 0.0 ,  0.0 );
                v5_eval = v5_pi->Eval( particle.pT , 0.0 , 0.0 ,  0.0 );
            }
            else if(n==2||n==3){ // get kaon vN
                v1_eval = v1_K->Eval( particle.pT , 0.0 , 0.0 ,  0.0 );
                v2_eval = v2_K->Eval( particle.pT , 0.0 , 0.0 ,  0.0 );
                v3_eval = v3_K->Eval( particle.pT , 0.0 , 0.0 ,  0.0 );
                v4_eval = v4_K->Eval( particle.pT , 0.0 , 0.0 ,  0.0 );
                v5_eval = v5_K->Eval( particle.pT , 0.0 , 0.0 ,  0.0 );
            }
            else if(n==4||n==5){ // get proton vN
                v1_eval = v1_P->Eval( particle.pT , 0.0 , 0.0 ,  0.0 );
                v2_eval = v2_P->Eval( particle.pT , 0.0 , 0.0 ,  0.0 );
                v3_eval = v3_P->Eval( particle.pT , 0.0 , 0.0 ,  0.0 );
                v4_eval = v4_P->Eval( particle.pT , 0.0 , 0.0 ,  0.0 );
                v5_eval = v5_P->Eval( particle.pT , 0.0 , 0.0 ,  0.0 );
            }

            if(v2_eval<0.0) v2_eval =0;
            if(v3_eval<0.0) v3_eval =0;
            if(v4_eval<0.0) v4_eval =0;
            if(v5_eval<0.0) v5_eval =0;

           
            // turn off harmonics which are not desired 
            Double_t v1 = doVn[0]*v1_eval;
            Double_t v2 = doVn[1]*v2_eval;
            Double_t v3 = doVn[2]*v3_eval;
            Double_t v4 = doVn[3]*v4_eval;
            Double_t v5 = doVn[4]*v5_eval;

            Double_t Max_dNdPhi = (1.0+2.0*(TMath::Abs(v1)+TMath::Abs(v2)+TMath::Abs(v3)+TMath::Abs(v4)+TMath::Abs(v5)));
            while(1){
                particle.phi = Harmonics_Phi_Dist_Rand->Uniform(0.0,2.0*TMath::Pi());
                Double_t von_nuemann_rejection = Harmonics_Phi_Dist_Rand->Uniform(0.0,Max_dNdPhi);
                Double_t dNdPhi = (1.0+2.0*(v1*TMath::Cos(particle.phi-event.psi1)
                                            +v2*TMath::Cos(2.0*(particle.phi-event.psi2))
                                            +v3*TMath::Cos(3.0*(particle.phi-event.psi3))
                                            +v4*TMath::Cos(4.0*particle.phi-event.psi4))
                                            +v5*TMath::Cos(5.0*particle.phi-event.psi5));
                if (von_nuemann_rejection<dNdPhi) break;
            }

            particle.px = particle.pT*TMath::Cos(particle.phi);
            particle.py = particle.pT*TMath::Sin(particle.phi);
            particle.pz = particle.pT*TMath::SinH(particle.eta);
            particle.E = TMath::Sqrt(particle.px*particle.px+particle.py*particle.py+particle.pz*particle.pz+particle.mass*particle.mass);

            event.particles.push_back(particle);

        } // end species loop

    } // end particle loop

    return event;
   

}

TGEvent TennGen::NextAuAu(){

    TGEvent event;
    event.particles.clear();

    if(fixPsi[0]) event.psi1 = fixed_psi[0];
    else event.psi1 = Psi_1->Uniform( 0 , 2.0*TMath::Pi());
    if(fixPsi[1]) event.psi2 = fixed_psi[1];
    else event.psi2 = 0.0;
    if(fixPsi[2]) event.psi3 = fixed_psi[2];
    else event.psi3 = Psi_3->Uniform( 0 , 2.0*TMath::Pi());
    if(fixPsi[3]) event.psi4 = fixed_psi[3];
    else event.psi4 = 0.0;

    event.nparts = -1;
    event.centBin = fCentBin;
    event.etaRange = fEtaRange;
    event.collEn = fCollEn;

    while (1){
            event.nparts = Int_t((event.etaRange/0.5)*multiplicty_distro->GetRandom(Multi_Rand));
            if(SpectraCentBin(event.nparts) == event.centBin) break;
        }
    Int_t harmonic_centbin = HarmonicCentBin(event.nparts);

    
    const Double_t piplusRatios[4] = {0.396738386,0.400412797,0.400958581,0.40230616};
    const Double_t piminusRatios[4] = {0.396738386,0.400412797,0.400958581,0.40230616};
    const Double_t kplusRatios[4] = {0.063108127,0.061919505,0.060984334,0.059731231};
    const Double_t kminusRatios[4] = {0.06118953,0.059236326, 0.058838257,0.057484421};
    const Double_t proRatios[4] = {0.043099904,0.041486068,0.042385006,0.042604604};
    const Double_t pbarRatios[4] = {0.03295875,0.032404541,0.033371486,0.034295646};


    Int_t particle_yeilds[6] = {0,0,0,0,0,0};
    particle_yeilds[0] = Int_t(piminusRatios[event.centBin]*event.nparts);
    particle_yeilds[1] = Int_t(piplusRatios[event.centBin]*event.nparts);
    particle_yeilds[2] = Int_t(kminusRatios[event.centBin]*event.nparts);
    particle_yeilds[3] = Int_t(kplusRatios[event.centBin]*event.nparts);
    particle_yeilds[4] = Int_t(pbarRatios[event.centBin]*event.nparts);
    particle_yeilds[5] = Int_t(proRatios[event.centBin]*event.nparts);

    for(Int_t n = 0; n<6; n++)
    {
        for(Int_t ipart = 0; ipart<particle_yeilds[n];ipart++)
        {
            TGParticle particle;    
            if (n==0){
                particle.pdg = -211;
                particle.mass = 0.139570;
                particle.charge = -1;
                particle.pT = pTdist_piMinus->GetRandom(pT_Rand);
            }
            if (n==1){
                particle.pdg = 211;
                particle.mass = 0.139570;
                particle.charge = 1;
                particle.pT = pTdist_piPlus->GetRandom(pT_Rand);
            } 
            if (n==2){
                particle.pdg = -321;
                particle.mass = 0.493677;
                particle.charge = -1;
                particle.pT = pTdist_kMinus->GetRandom(pT_Rand);
            }
            if (n==3){
                particle.pdg = 321;
                particle.mass = 0.493677;
                particle.charge = 1;
                particle.pT = pTdist_kPlus->GetRandom(pT_Rand);
            }
            if (n==4){
                particle.pdg = -2212;
                particle.mass = 0.938272;
                particle.charge = -1;
                particle.pT = pTdist_pbar->GetRandom(pT_Rand);
            }
            if (n==5){
                particle.pdg = 2212;
                particle.mass = 0.938272;
                particle.charge = 1;
                particle.pT = pTdist_p->GetRandom(pT_Rand);
            }
            
            particle.eta = etadist->Uniform(-1.0*event.etaRange , 1.0*event.etaRange);
            Double_t v2 = doVn[1]*HarmonicFunction(particle.pT,harmonic_centbin,0, particle.pdg);
            Double_t v3 = doVn[2]*HarmonicFunction(particle.pT,harmonic_centbin,1, particle.pdg);
            Double_t v4 = doVn[3]*HarmonicFunction(particle.pT,harmonic_centbin,2, particle.pdg);
            Double_t v1 = doVn[0]*0.02*v2;

            Double_t Max_dNdPhi = (1.0+2.0*(TMath::Abs(v1)+TMath::Abs(v2)+TMath::Abs(v3)+TMath::Abs(v4)));
            while(1){
                particle.phi = Harmonics_Phi_Dist_Rand->Uniform(0.0,2.0*TMath::Pi());
                Double_t von_nuemann_rejection = Harmonics_Phi_Dist_Rand->Uniform(0.0,Max_dNdPhi);
                Double_t dNdPhi = (1.0+2.0*(v1*TMath::Cos(particle.phi-event.psi1)
                                            +v2*TMath::Cos(2.0*(particle.phi-event.psi2))
                                            +v3*TMath::Cos(3.0*(particle.phi-event.psi3))
                                            +v4*TMath::Cos(4.0*particle.phi-event.psi4)));
                if (von_nuemann_rejection<dNdPhi) break;
            }

            particle.px = particle.pT*TMath::Cos(particle.phi);
            particle.py = particle.pT*TMath::Sin(particle.phi);
            particle.pz = particle.pT*TMath::SinH(particle.eta);
            particle.E = TMath::Sqrt(particle.px*particle.px+particle.py*particle.py+particle.pz*particle.pz+particle.mass*particle.mass);

            event.particles.push_back(particle);

        } // end species loop

    } // end particle loop

    return event;
        

}

