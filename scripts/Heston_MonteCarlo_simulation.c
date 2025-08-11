#include <math.h>
#include <stdlib.h>
#include <time.h>


double normal_std() {
    double u1 = rand() / (RAND_MAX + 1.0);
    double u2 = rand() / (RAND_MAX + 1.0);
    if (u1 < 1e-10)
    u1 = 1e-10;
    return sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);
}

__declspec(dllexport)
void Heston(double S0, double v0, double mu, double kappa, double theta, double sigma, double rho, double T, double *vol, double *precios, int N, int seed) {
    srand(seed);
    double dt = T / N;
    double S = S0;
    double v = v0 * v0;

    precios[0] = S;
    vol[0] = v0;

    for (int i = 1; i <= N; ++i) {
        double z1 = normal_std();
        double z2 = normal_std();
        double w1 = z1;
        double w2 = rho * z1 + sqrt(1 - rho * rho) * z2;

        double sqrt_v = sqrt(fmax(v, 0.0));

        S = S * exp((mu - 0.5 * v) * dt + sqrt_v * sqrt(dt) * w1);
        v = v + kappa * (theta - v) * dt + sigma * sqrt_v * sqrt(dt) * w2;
        v = fmax(v, 0.0);

        precios[i] = S;
        vol[i] = sqrt(v);
    }
}

