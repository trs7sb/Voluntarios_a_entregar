#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>

// Definimos constantes
#define h 0.01 // Espaciado espacial
#define PI 3.141592653589793
#define N 5000 // Número de puntos espaciales
#define I _Complex_I // Definición de la unidad imaginaria
int T = 3000; // Número de pasos temporales

// Inicializa la función de onda Phi
void inicializar_phi(double complex *Phi, double n_ciclos) {
    double k_tilde = (2.0 * PI * n_ciclos) / N;
    double x0 = (N * h) / 4.0;
    double sigma = (N * h) / 16.0;

    for (int j = 1; j < N; j++) {
        double x = j * h;
        Phi[j] = cexp(I * k_tilde * j) * cexp(-pow((x - x0), 2) / (2 * sigma * sigma));
    }
    Phi[0] = 0.0;
    Phi[N] = 0.0;

    // Normalización de Phi
    double norma = 0.0;
    for (int j = 0; j <= N; j++) {
        norma += pow(cabs(Phi[j]), 2) * h;
    }
    norma = sqrt(norma);
    for (int j = 1; j < N; j++) {
        Phi[j] /= norma;
    }
}

// Inicializa el potencial V con n_barreras de ancho 200 y separación 200
void inicializar_potencial_multibarrera(double *V, double lambda, double k_tilde, int n_barreras, int *inicio_detector) {
    double V_j_h2 = lambda * k_tilde * k_tilde;
    int ancho = 200;
    int separacion = 200;
    int inicio, fin;

    for (int j = 0; j <= N; j++) {
        V[j] = 0.0;
    }

    int offset = (int)(N * h / 4.0 / h) + N/16;
    for (int b = 0; b < n_barreras; b++) {
        inicio = offset + b * (ancho + separacion);
        fin = inicio + ancho - 1;
        if (fin > N) break;
        for (int j = inicio; j <= fin && j <= N; j++) {
            V[j] = V_j_h2;
        }
    }
    // El detector empieza justo después de la última barrera
    *inicio_detector = offset + n_barreras * (ancho + separacion);
}

int maximo_global(double *PD, int T) { 
    int max_index = 0;
    double max_value = PD[0];
    for (int n = 1; n < T; n++) {
        if (PD[n] > max_value) {
            max_value = PD[n];
            max_index = n;
        }
    }
    return max_index;
}

void calcula_alpha(double complex *Phi, double complex *a, double complex *b, double complex *c, double complex *gamma, double complex *alpha) {
    // Calculamos alpha fuera del bucle de tiempo
    alpha[N-1] = 0.0;

    for (int j = N-1; j > 0; j--) {
        gamma[j] = b[j] + c[j] * alpha[j];
        alpha[j-1] = -a[j] / gamma[j];
    }
}

void calcular_coeficientes_tridiagonal(double complex *Phi, double complex *c, double complex *gamma, double complex *beta, double complex *chi, double tilde_s, double complex *alpha) {
    // COEFICIENTES TRIDIAGONAL
   double complex d[N + 1];

    for (int j = 1; j < N; j++) {
       d[j] = 4.0 * I / tilde_s * Phi[j];
   }
   beta[N - 1] = 0.0 + 0.0*I;    
   // Orden descendente
   for (int j = N - 1; j > 0; j--) {
       beta[j - 1] = (d[j] - c[j] * beta[j]) / gamma[j];
   }

   // Orden ascendente
   chi[0] = 0.0;
   chi[N] = 0.0;
   for (int j = 0; j < N - 1; j++) {
       chi[j + 1] = alpha[j]*chi[j] + beta[j];
   }
}


int main() {
    srand((unsigned)time(NULL));

    int nexp = 0; 
    int m_T = 0; 
    double PD[T]; 
    int n_D = 0;  
    double PD_nD = 0.0;

    double complex Phi[N + 1], Phi_next[N + 1], chi[N + 1];
    double V[N + 1];
    complex double a[N + 1], b[N + 1], c[N + 1];
    complex double gamma[N - 1];

    // Parámetros iniciales
    int n_ciclos = N / 16;
    double lambda = 0.3;
    double k_tilde = (2.0 * PI * n_ciclos) / N;
    double tilde_s = 1.0 / (4.0 * k_tilde * k_tilde);
    complex double alpha[N];
    complex double beta[N]; 

  
    int separacion = 200;
    int ancho_detector = 200;
    int ancho = N/5;
    int gauss = (int)(N * h / 4.0 / h) + N/16; // o el valor que uses para el paquete gaussiano

    int n_barreras = (N - gauss - ancho_detector) / (ancho + separacion);
    if (n_barreras < 1) n_barreras = 1;

    FILE *f_pd = fopen("probabilidades_pd.txt", "w");
    if (f_pd == NULL) {
        fprintf(stderr, "Error al abrir el fichero para guardar PD.\n");
        return 1;
    }

    for (int experimento = 0; experimento <= nexp; experimento++) {
        a[0] = 0.0;
        b[0] = 0.0;
        c[0] = 0.0;

        inicializar_phi(Phi, n_ciclos);
        int inicio_detector = 0;
        inicializar_potencial_multibarrera(V, lambda, k_tilde, n_barreras, &inicio_detector);

        FILE *f_normas = fopen("normas.txt", "w");
        if (f_normas == NULL) {
            fprintf(stderr, "Error al abrir el fichero para guardar las normas.\n");
            return 1;
        }

        FILE *f_datos = fopen("schrodinger.dat", "w");
        if (f_datos == NULL) {
            fprintf(stderr, "Error al abrir el fichero para guardar los datos.\n");
            return 1;
        }

        for (int j = 1; j < N; j++) {
            a[j] = 1.0;
            b[j] = -2.0 + 2.0 * I / tilde_s - V[j];
            c[j] = 1.0;
        }

        calcula_alpha(Phi, a, b, c, gamma, alpha);

        int j;
        for (int n = 0; n < T; n++) {
            for (j = 0; j <= N; j++) {
                double x = j * h;
                double densidad = pow(cabs(Phi[j]), 2);
                fprintf(f_datos, "%.8f,%.8f,%.8f\n", x, densidad, V[j]);
            }
            fprintf(f_datos, "\n");

        
            calcular_coeficientes_tridiagonal(Phi, c, gamma, beta, chi, tilde_s, alpha);


            for (int j = 1; j < N; j++) {
                Phi_next[j] = chi[j] - Phi[j];
            }
            for (int j = 1; j < N; j++) {
                Phi[j] = Phi_next[j];
            }

            double norma = 0.0;
            for (int j = 1; j < N; j++) {
                norma += pow(cabs(Phi[j]), 2) * h;
            }
            fprintf(f_normas, "Paso %d: Norma = %.10f\n", n, norma);

    
            if (experimento == 0) {
                double pd = 0.0;
                for (int j = inicio_detector; j < inicio_detector + ancho_detector && j <= N; j++) {
                    pd += pow(cabs(Phi[j]), 2) * h;
                }
                PD[n] = pd;
                fprintf(f_pd, "%d %.20f\n", n, PD[n]);
            }
        }
        fclose(f_normas);
        fclose(f_datos);

         // Buscar el primer máximo local de PD(n)
         if(experimento==0) {
            n_D = maximo_global(PD, T);
          printf("Primer máximo local de PD(n) en n = %d\n", n_D);
           PD_nD = PD[n_D];
             printf("PD(n_D) = %.20f\n", PD_nD);
        }

        //T = n_D;
        if (experimento > 0) {
            double p = (double)rand() / ((double)RAND_MAX + 1.0);
            if (p < PD_nD) {
                m_T++;
            }
        }
    }
    fclose(f_pd);

    double K = (double)m_T / nexp;
    printf("Número de transmisiones m_T = %d\n", m_T);
    printf("Coeficiente de transmisión K = %.6f\n", K);

    return 0;
}