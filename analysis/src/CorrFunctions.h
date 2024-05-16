#ifndef CORRFUNCTIONS_H
#define CORRFUNCTIONS_H

#include <TComplex.h>
#include <TMath.h>



namespace CorrFunctions 
{


static double TwoPartRefWeight(int M)
{
    double m = 1.0*M;
    double w = m*(m-1);
    return  w;
}

static double FourPartRefWeight(int M)
{
    double m = 1.0*M;
    double w = m*(m-1)*(m-2)*(m-3);
    return w;
}

static double TwoPartDiffWeight(int M, int mp, int mq=0)
{
    double m = 1.0*M;
    double p = 1.0*mp;
    double q = 1.0*mq;
    double w = (p*m) - q;
    return w;
}

static double FourPartDiffWeight(int M, int mp, int mq=0)
{
    double m = 1.0*M;
    double p = 1.0*mp;
    double d = 1.0*mq;
    double w = (p*m - 3.0*d)*(m-1)*(m-2);
    return w;
}

static double TwoPartRef(TComplex Q, int M)
{
    // return (Q.Rho2() - M) / (M * (M - 1));
    double m = 1.0*M;
    double A = Q.Rho2();
    double w = CorrFunctions::TwoPartRefWeight(M);
    double val = A - m;
    return val/w;
}

static double FourPartRef(TComplex Q, TComplex Q2, int M)
{
    // double A = (Q.Rho2()*Q.Rho2() + Q2.Rho2() - 2.0*(Q2*TComplex::Conjugate(Q)*TComplex::Conjugate(Q)).Re())/(M*(M-1)*(M-2)*(M-3));
    // double B = (2.0*(M-2)*Q.Rho2() - 1.0*M*(M-3))/(1.0*M*(M-1)*(M-2)*(M-3));

    TComplex qstar = TComplex::Conjugate(Q);
    TComplex q2star = TComplex::Conjugate(Q2);

    TComplex qmagfour = Q*Q*qstar*qstar;
    TComplex qmagtwo = Q*qstar;
    TComplex q2magtwo = Q2*q2star;
    TComplex q2qstarqstar = Q2*qstar*qstar;

    double qmagfour_real = qmagfour.Re();
    double qmagtwo_real = qmagtwo.Re();
    double q2magtwo_real = q2magtwo.Re();
    double q2qstarqstar_real = q2qstarqstar.Re();
    double m = 1.0*M;
    double A = (qmagfour_real + q2magtwo_real - 2.0*q2qstarqstar_real);
    double B = -2.0*(2.0*(m-2)*qmagtwo_real - m*(m-3));
    double w = CorrFunctions::FourPartRefWeight(M);
    double val = A + B;
    return val/w;
}

static double TwoPartDiff(TComplex Q, TComplex p, int mp, int M, int mq =0)
{
    TComplex qstar = TComplex::Conjugate(Q);
    TComplex pqstar = p*qstar;
    double w = CorrFunctions::TwoPartDiffWeight(M, mp, mq);

    double A = pqstar.Re();
    double MQ = 1.0*mq;
    double val = A - MQ;
    return val / w;

    // return (p*TComplex::Conjugate(Q)) / (1.0*mp*M);
}

static double FourPartDiff(TComplex Q, TComplex Q2, TComplex p, int mp, int M, TComplex u=0, TComplex u2=0, int mq = 0)
{
    TComplex Qstar = TComplex::Conjugate(Q);
    TComplex Q2star = TComplex::Conjugate(Q2);
    

    TComplex A = p*Q*Qstar*Qstar;
    TComplex B = p*Q*Q2star;
    TComplex C = p*Qstar;
   

    double A_real = A.Re();
    double B_real = B.Re();
    double C_real = C.Re();
 
    double m = 1.0*M;
    double m_p = 1.0*mp;

    double val = A_real - B_real - 2.0*m*C_real + 2*C_real;
    double w = CorrFunctions::FourPartDiffWeight(M, mp, mq);

    if(mq != 0)
    {

        TComplex ustar = TComplex::Conjugate(u);
        TComplex Au = u2*Qstar*Qstar;
        TComplex Bu = u*Qstar;
        TComplex mag_Q_sqrd = Q*Qstar;
        TComplex Cu = Q*ustar;
        TComplex Du = u2*Q2star;

        double Au_real = Au.Re();
        double Bu_real = Bu.Re();
        double Cu_real = Cu.Re();
        double Du_real = Du.Re();
        double mag_Q_sqrd_real = mag_Q_sqrd.Re();

        double m_q = 1.0*mq;

        double val_u = -1.0*Au_real - 2.0*m_q*mag_Q_sqrd_real + 7.0*Bu_real - Cu_real + Du_real + 2.0*m_q*m - 6.0*m_q;
        val += val_u;
    }

    return val / w;
    // return (p*Q*TComplex::Conjugate(Q)*TComplex::Conjugate(Q) 
    // - p*Q*TComplex::Conjugate(Q2) 
    // - 2.0*M*p*TComplex::Conjugate(Q) 
    // + 2.0*p*TComplex::Conjugate(Q)) / (1.0*mp*M*(M-1)*(M-2));
}

static TComplex pvec(double phi, int n)
{
    double N = 1.0*n;
    TComplex p(TMath::Cos(N*phi), TMath::Sin(N*phi));
    return p;
}

}

#endif // CORRFUNCTIONS_H