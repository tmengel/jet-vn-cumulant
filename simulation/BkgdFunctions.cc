#include "BkgdFunctions.h"


#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>

namespace toymodel  
{

double BWintegrand(double *x, double *p)
{

  double x0 = x[0]; 
  double mass     = p[0];
  double pT       = p[1];
  double beta_max = p[2];
  double temp     = p[3];
  double n      = p[4];

  // Keep beta within reasonable limits
  double beta = beta_max * TMath::Power(x0, n);
  if (beta > 0.9999999999999999) beta = 0.9999999999999999;

  double mT      = TMath::Sqrt(mass*mass+pT*pT);

  double rho0   = TMath::ATanH(beta);  
  double arg00 = pT*TMath::SinH(rho0)/temp;
  if (arg00 > 700.) arg00 = 700.; // avoid FPE
  double arg01 = mT*TMath::CosH(rho0)/temp;
  double f0 = x0*mT*TMath::BesselI0(arg00)*TMath::BesselK1(arg01);
  return f0;
}

double MyStaticBGdNdPtTimesPt(double *x, double *p)
{

  double pT = x[0];;
  double mass    = p[0];
  double beta    = p[1];
  double temp    = p[2];
  double n       = p[3];
  double norm    = p[4];

  static TF1 * fIntBG = 0;
  if(!fIntBG)
  {
    fIntBG = new TF1 ("fIntBG", BWintegrand, 0, 1, 5);
  }

  fIntBG->SetParameters(mass, pT, beta, temp,n);
  double result = fIntBG->Integral(0,1);
  return result*norm*pT;

}

double BGBW(double *x, double *p)
{
  return MyStaticBGdNdPtTimesPt(x,p);
}

double VnFunction(double *x, double *p)
{
  // TF1 * vn_pi = new TF1("vn_pi", "(x<[5])*((x>[6])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) +(x<[6])*(0)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );

  double pT = x[0];

  double p0 = p[0];
  double p1 = p[1];
  double p2 = p[2];
  double p3 = p[3];
  double p4 = p[4];
  double upper_lim = p[5];
  double lower_lim = p[6];

  double vn = 0.0;
  if (pT > upper_lim) pT = upper_lim;
  if (pT > lower_lim)
  {
    vn = p0 + (pT*p1) + (pT*pT*p2) + (pT*pT*pT*p3) + (pT*pT*pT*pT*p4);
  }
  return vn; 
}

double TruthRefVn(int harmonic, double min_particle_pt)
{

    // ================================================================
    // ===================== Function Parmeters ========================
    // ================================================================

    // Vn functions
    TF1 * vn_pi =  new TF1("vn_pi", VnFunction, min_particle_pt, 100, 7);
    TF1 * vn_k =  new TF1("vn_k", VnFunction, min_particle_pt, 100, 7);
    TF1 * vn_pro =  new TF1("vn_pro", VnFunction, min_particle_pt, 100, 7);

    std::vector<TF1*> vn_functions = {vn_pi, vn_k, vn_pro};
    for (auto tf1 : vn_functions){ tf1->SetNpx(10000); }
    for (int i =0; i < 7; i++)
    {
        vn_pi->SetParameter(i,pi_vn_params[harmonic][i]);
        vn_k->SetParameter(i,k_vn_params[harmonic][i]);
        vn_pro->SetParameter(i,pro_vn_params[harmonic][i]);
    }

    // pT distributions
    TF1 * f1_piPluspT = new TF1("f1_piPluspT_tmp", BGBW, min_particle_pt, 100, 5);
    TF1 * f1_piMinuspT = new TF1("f1_piMinuspT_tmp", BGBW, min_particle_pt, 100, 5);
    TF1 * f1_KPluspT = new TF1("f1_KPluspT_tmp", BGBW, min_particle_pt, 100, 5);
    TF1 * f1_KMinuspT = new TF1("f1_KMinuspT_tmp", BGBW, min_particle_pt, 100, 5);
    TF1 * f1_protonpT = new TF1("f1_protonpT_tmp", BGBW, min_particle_pt, 100, 5);
    TF1 * f1_pbarpT = new TF1("f1_pbarpT_tmp", BGBW, min_particle_pt, 100, 5);
    std::vector<TF1*> pT_distos = {f1_piPluspT, f1_piMinuspT, f1_KPluspT, f1_KMinuspT, f1_protonpT, f1_pbarpT};
    for (auto tf1 : pT_distos){ tf1->SetNpx(10000); }
    for (int i=0;i<6;i++)
    {
        for (int j=0;j<5;j++)
        {
            pT_distos[i]->SetParameter(j,blastwave_params[j][i]);
        }
    }
    
    // ================================================================
    // ===================== Initialize Histos ========================
    // ================================================================
    std::vector<TH1D*> h1_dNdpT;
    std::vector<TH1D*> h1_vn;
    std::vector<int> vn_to_dndpt_map = {0, 0, 1, 1, 2, 2};

    for (auto f1: pT_distos)
    {
        TH1D * h1 = (TH1D*)f1->CreateHistogram();
        h1->SetName(f1->GetName());
        h1_dNdpT.push_back(h1);
    }
    for (auto tf1 : vn_functions)
    {
        TH1D * h1 = (TH1D*)tf1->CreateHistogram();
        h1->SetName(tf1->GetName());
        h1_vn.push_back(h1);
    }

    // Scale by particle yeilds
    for (unsigned int i = 0; i < n_species; i++)
    {
        h1_dNdpT.at(i)->Scale(1.0/h1_dNdpT.at(i)->Integral());
        h1_dNdpT.at(i)->Scale(particle_yeilds.at(i));
    }

    double dNdpT_total = 0.0;
    for (auto h1 : h1_dNdpT)
    {
        dNdpT_total += h1->Integral();
    }

    double vn_dNdpT_total = 0.0;
    for (unsigned int i = 0; i < n_species; i++)
    {
        TH1D * h1 = (TH1D*)h1_vn.at(vn_to_dndpt_map.at(i))->Clone();
        TH1D * dNdpT = (TH1D*)h1_dNdpT.at(i)->Clone();
        h1->Multiply(dNdpT);
        vn_dNdpT_total += h1->Integral();
    }

    double reference_vn = vn_dNdpT_total / dNdpT_total;
    std::cout << "truth reference v" << harmonic+2 << ": " << reference_vn << std::endl;

    return reference_vn;


}

} // namespace toymodel