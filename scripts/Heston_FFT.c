#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "fftw3.h"

complex double Hestoncf(complex double phi, double S, double T, double kappa,
                        double rho, double sigma, double theta, double v0,
                        double r, double div){
    double xx = log(S);
    complex double gamma = kappa - rho*sigma*phi*I;
    complex double z = -0.5 * phi*(phi + I);
    complex double psi  = csqrt(gamma*gamma - 2*sigma*sigma*z);


    complex double CCaux = (2.0 * psi - (psi - gamma) * (1.0 - cexp(-psi * T))) / (2.0 * psi);
    complex double CC = (-(kappa * theta) / (sigma*sigma)) * (2 * log(CCaux) + (psi - gamma) * T);

    complex double BB = (2.0 * z * (1.0 - cexp(-psi * T)) * v0) / (2.0 * psi - (psi - gamma) * (1.0 - cexp(-psi * T)));
    complex double AA = I * phi * (xx + r * T);

    return cexp(AA + BB + CC);
}

complex double HestoncfMod(double phi, double S, double T, double kappa,
                        double rho, double sigma, double theta, double v0,
                        double r, double div, double alpha){
    complex double Newphi = phi - (alpha + 1.0)*I;
    complex double CF = Hestoncf(Newphi, S, T,  kappa,rho, sigma,  theta,  v0, r, div);

    complex double NumCF = cexp(-r * T) * CF;
    complex double DenCF = alpha*(alpha + 1.0) - phi*phi + I*phi * (2.0*alpha + 1.0);


    return NumCF/DenCF;
}

__declspec(dllexport)
void CallHestonFFT(double Spot, double Maturity, double kappa, double rho, double volvol,
                   double theta, double var0, double rate, double div,
                   double *Strikes, int numStrikes, double *CallPrices) {

    int NN = 4096;
    double alpha = 1.25;
    double cc = 1000.0;
    double eta = cc / NN;
    double Lambda = (2.0 * M_PI) / (NN * eta);
    double bb = (NN * Lambda) / 2.0;

    // Allocate arrays
    double *ku = malloc(NN * sizeof(double));
    double *NewPhi = malloc(NN * sizeof(double));
    fftw_complex *ModCF = fftw_malloc(sizeof(fftw_complex) * NN);
    fftw_complex *Simpson = fftw_malloc(sizeof(fftw_complex) * NN);
    fftw_complex *FuncFFT = fftw_malloc(sizeof(fftw_complex) * NN);
    fftw_complex *Payoff = fftw_malloc(sizeof(fftw_complex) * NN);

    // Fill ku and NewPhi
    for (int j = 0; j < NN; j++) {
        ku[j] = -bb + Lambda * j;
        NewPhi[j] = eta * j;
    }

    // Compute ModCF
    for(int i=0; i<NN; i++)
        ModCF[i] = HestoncfMod( NewPhi[i], Spot, Maturity, kappa, rho, volvol, theta, var0, rate, div, alpha);

    // Simpson weights
    for (int j = 0; j < NN; j++) {
        double delta = (j == 0) ? 1.0 : 0.0;
        Simpson[j] = (eta / 3.0) * (3.0 + cpow(-I, j) - delta);
    }

    // FuncFFT
    for (int j = 0; j < NN; j++) {
        FuncFFT[j] = cexp(I * bb * NewPhi[j]) * ModCF[j] * Simpson[j];
    }

    // FFT
    fftw_plan plan = fftw_plan_dft_1d(NN, FuncFFT, Payoff, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    // Final call values
    double *CallValues = malloc(NN * sizeof(double));
    for (int j = 0; j < NN; j++) {
        CallValues[j] = creal(Payoff[j]) * exp(-alpha * ku[j]) / M_PI;
    }

    // Interpolate for requested strikes
    for (int s = 0; s < numStrikes; s++) {
        double logK = log(Strikes[s]);
        int pos = (int)(((logK + bb) / Lambda) + 1.0);
        if (pos >= 0 && pos < NN) {
            CallPrices[s] = CallValues[pos];
        } else {
            CallPrices[s] = 0.0; // Out of bounds
        }
    }

    // Cleanup
    fftw_destroy_plan(plan);
    fftw_free(ModCF);
    fftw_free(Simpson);
    fftw_free(FuncFFT);
    fftw_free(Payoff);
    free(CallValues);
    free(ku);
    free(NewPhi);
}


int main(){
    int NNStrikes = 10;
    double Strikes[NNStrikes], prices[NNStrikes];
    for(int i=0; i< NNStrikes; i++)
        Strikes[i] = 100 + (i-5)*5;
    CallHestonFFT(100.0, 1.0, 1.5,-0.4, 0.4, 0.03, 0.014,0.02,  0, Strikes, NNStrikes, prices);
    for(int j=0; j< NNStrikes; j++)
        printf("%lf\t%lf\n", Strikes[j],prices[j]);
    return 0;
}
