#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>


#define h 0.01 
#define PI 3.141592653589793
#define N 1000
#define I _Complex_I 
int T = 2000; 


// Función para inicializar el potencial V
void inicializar_potencial(double *V, double lambda, double k_tilde) {
    double V_j_h2 = lambda * k_tilde * k_tilde; // Altura del potencial
    int inicio = (2 * N) / 5;
    int fin = (3 * N) / 5;

    for (int j = 0; j <= N; j++) {
        if (j >= inicio && j <= fin) {
            V[j] = V_j_h2; 
        } else {
            V[j] = 0.0; 
        }
    }
}

// Función para inicializar la función de onda Phi
void inicializar_funcion_onda(double complex *Phi, double n_ciclos) {

    double k_tilde = (2.0 * PI * n_ciclos) / N; 
    double x0 = (N * h) / 4.0;                 
    double sigma = (N * h) / 16.0;         

    for (int j = 1; j < N; j++) {
        double x = j * h; 
        Phi[j] = cexp(I * k_tilde * j) * cexp(-pow((x - x0), 2) / (2 * sigma * sigma)); // Φ_{j,0} = e^{i\tilde{k}_0 x} e^{-(x - x_0)^2 / (2σ^2)}
    }

    Phi[0] = 0.0;
    Phi[N] = 0.0;

    // Normalización de Phi
    double norma = 0.0;
    for (int j = 0; j <= N; j++) {
        norma += pow(cabs(Phi[j]), 2) * h;
    }
    norma = sqrt(norma);

    for (int j = 1; j <N; j++) {
        Phi[j] /= norma; 
    }
}


// Cambiar la función para encontrar el máximo global
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

void Valor_esperado_j(double complex *Phi, double *promedio, double *var) {
    double suma = 0.0, suma2 = 0.0;

    for (int j = 0; j <= N; j++) {
        double p = pow(cabs(Phi[j]), 2) * h;
        suma += j * p; // Valor esperado de j
        suma2 += j * j * p;
    }
    *promedio = suma;
    *var = suma2 - suma * suma;
}

void Valor_esperado_cinetica(double complex *Phi, double *media, double *varianza) {
    double suma = 0.0, suma2 = 0.0;
    for (int j = 2; j < N-1; j++) {
        // Segunda dervida discreta
        double complex segunda_derivada = (Phi[j+1] - 2.0*Phi[j] + Phi[j-1]) / (h*h);
        double energia = -0.5 * creal(conj(Phi[j]) * segunda_derivada) * h;
        suma += energia;
        // Cuarta derivada central discreta 
        double complex cuarta_derivada = (Phi[j-2] - 4.0*Phi[j-1] + 6.0*Phi[j] - 4.0*Phi[j+1] + Phi[j+2]) / pow(h, 4);
        double energia2 = 0.25 * creal(conj(Phi[j]) * cuarta_derivada) * h;
        suma2 += energia2;
    }
    *media = suma;
    *varianza = suma2 - suma * suma;
}

int main() {
    srand((unsigned)time(NULL));

    int nexp= 0; // Número de experimentos
    int j;
    int m_T = 0;  // Número de transmisiones
    double PD[T]; 
    int n_D = 0;  // Paso del primer máximo global
    double PD_nD = 0.0;

    double complex Phi[N + 1], Phi_nueva[N + 1], chi[N + 1];
    double V[N + 1];
    complex double a[N + 1], b[N + 1], c[N + 1];
    complex double gamma[N -1];
    complex double alpha[N];
    complex double beta[N]; 

       // Parámetros iniciales
    int n_ciclos = N/16; // Número de ciclos
    double lambda = 0.3; // Altura del potencial
    double k_tilde = (2.0 * PI * n_ciclos) / N; // \tilde{k}_0

    double tilde_s = 1.0 / (4.0 * k_tilde * k_tilde);

    FILE *f_pd = fopen("probabilidades_pd.txt", "w");
    if (f_pd == NULL) {
    fprintf(stderr, "Error al abrir el fichero para guardar PD.\n");
    return 1;
    } 

    for (int experimento = 0; experimento <=nexp; experimento++) {
    a[0] = 0.0; 
    b[0] = 0.0; 
    c[0] = 0.0; 

    inicializar_funcion_onda(Phi, n_ciclos);
    inicializar_potencial(V, lambda, k_tilde);

     FILE *f_norma = fopen("norma.txt", "w");
     if (f_norma == NULL) {
         fprintf(stderr, "Error al abrir el fichero para guardar las normas.\n");
         return 1;
     }

    FILE *f_datos = fopen("schrodinger.dat", "w");
    if (f_datos == NULL) {
        fprintf(stderr, "Error al abrir el fichero para guardar los datos.\n");
        return 1;
    }

    FILE *f_esperadoj = fopen("esperado_j.dat", "w");
    if (f_datos == NULL) {
        fprintf(stderr, "Error al abrir el fichero para guardar los datos.\n");
        return 1;
    }
    FILE *f_esperadoE = fopen("esperado_E.dat", "w");
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
    
    for (int n = 0; n < T; n++) {

        for (j = 0; j <= N; j++) {
            double x = j * h;
            double densidad = pow(cabs(Phi[j]), 2);
            fprintf(f_datos, "%.8f,%.8f,%.8f\n", x, densidad, V[j]);
        }
        fprintf(f_datos, "\n"); 

        calcular_coeficientes_tridiagonal(Phi, c, gamma, beta, chi, tilde_s, alpha);

        for (int j = 1; j < N; j++) {
            Phi_nueva[j] = chi[j] - Phi[j];
        }

        // Actualizamos Phi
        for (int j = 1; j < N; j++) {
            Phi[j] = Phi_nueva[j];
        }

        // Calculamos la norma de la función de onda
        double norma = 0.0;
        for (int j = 1; j <N; j++) {
            norma += pow(cabs(Phi[j]), 2) * h; 
        }

        // Guardamos la norma en el fichero
        fprintf(f_norma, "Paso %d: Norma = %.10f\n", n, norma);
        

        if(experimento==0) {
            double pd = 0.0;
            double promedioj;
            double varj;
            double promedioE;
            double varE;
            for (int j = (4*N)/5; j <= N; j++) {
                pd += pow(cabs(Phi[j]), 2) * h;
            }
            PD[n] = pd;

             fprintf(f_pd, "%d %.20f\n", n, PD[n]);

           Valor_esperado_j(Phi, &promedioj, &varj);
           Valor_esperado_cinetica(Phi, &promedioE, &varE);

            fprintf(f_esperadoj, "%d %.10f %.10f\n", n, promedioj, sqrt(varj));
            fprintf(f_esperadoE, "%d %.10f %.10f\n", n, promedioE, sqrt(varE));
        }

    }
        fclose(f_norma);
        fclose(f_datos);
        fclose(f_esperadoj);
        fclose(f_esperadoE);

        // Buscar el primer máximo global de PD(n)
        if(experimento==0) {
            n_D = maximo_global(PD, T);
          printf("Primer máximo global de PD(n) en n = %d\n", n_D);
           PD_nD = PD[n_D];
             printf("PD(n_D) = %.20f\n", PD_nD);
        }
       
        T=n_D;

        int detectada_derecha;

       if(experimento>0) {

        double p = (double)rand() / ((double)RAND_MAX + 1.0);

             if (p < PD_nD) {
                 m_T++;
                 detectada_derecha=1;
            }else 
            {
                 detectada_derecha=0;
             };
     }
    

}
     fclose(f_pd);

    double K = (double)m_T / nexp;
    printf("Número de transmisiones m_T = %d\n", m_T);
    printf("Coeficiente de transmisión K = %.6f\n", K);


    return 0;
}




