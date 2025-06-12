#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#define N 100 // Tamaño de la red en la dirección x
#define PASOS 10000

double T = 1;

// Función para inicializar la red con magnetización inicial nula de forma más sencilla
void inicializarRed_aleatoria(int red[N][N]) {
    int espines_totales = (N - 2) * N; // Excluyendo las filas fijas
    int mitad = espines_totales / 2;

    // Crear un array temporal para almacenar los espines
    int espines[espines_totales];
    for (int i = 0; i < mitad; i++) {
        espines[i] = 1;  // Mitad positivos
    }
    for (int i = mitad; i < espines_totales; i++) {
        espines[i] = -1; // Mitad negativos
    }

    // Mezclar los espines aleatoriamente
    for (int i = espines_totales - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int temp = espines[i];
        espines[i] = espines[j];
        espines[j] = temp;
    }

    // Asignar los espines mezclados a la red
    int index = 0;
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            if (x == 0) {
                red[x][y] = -1; // Borde superior con espines negativos
            } else if (x == N - 1) {
                red[x][y] = 1; // Borde inferior con espines positivos
            } else {
                red[x][y] = espines[index++];
            }
        }
    }
}

void inicializarRed_mitad(int red[N][N]) {

    int espines_totales = (N-2) * N;
    int mitad = espines_totales / 2;
    int positivos = 0, negativos = 0;

    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            if (x == 0) {
                red[x][y] = -1; // Borde superior con espines negativos
            } else if (x == N - 1) {
                red[x][y] = 1; // Borde inferior con espines positivos
            } else {
                if (positivos < mitad) {
                    red[x][y] = -1;
                    positivos++;
                } else {
                    red[x][y] = 1;
                    positivos++;
                }
            }
        }
    }

}


// Imprime la red 
void escribirRed(FILE *file, int red[N][N]) {
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            fprintf(file, "%d", red[x][y]);
            if (y < N - 1) {
                fprintf(file, ", ");
            }
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
}

double calcularEnergia_Local(int red[N][N], int x1, int y1, int x2, int y2 ) {
    int E;
    E=0;

    // Vecinos del primer espín (x1, y1)
    int arriba1, abajo1, izquierda1, derecha1;


    // Vecino de arriba
    if (x1 == 1) {
        arriba1 = -1; // Bordes fijos en "x=0"
    } else {
        arriba1 = red[x1-1][y1];
    }

    // Vecino de abajo
    if (x1 == N - 2) {
        abajo1 = 1; // Bordes fijos en "x=N-1"
    } else {
        abajo1 = red[x1+1][y1];
    }

    // Vecino de la izquierda (condiciones periódicas en "y")
    if (y1 == 0) {
        izquierda1 = red[x1][N-1];
    } else {
        izquierda1 = red[x1][y1-1];
    }

    // Vecino de la derecha (condiciones periódicas en "y")
    if (y1 == N - 1) {
        derecha1 = red[x1][0];
    } else {
        derecha1 = red[x1][y1+1];
    }

    // Vecinos del segundo espín (x2, y2)
    int arriba2, abajo2, izquierda2, derecha2;

    // Vecino de arriba
    if (x2== 1) {
        arriba2 = -1; // Bordes fijos en "x=0"
    } else {
        arriba2 = red[x2-1][y2];
    }

    // Vecino de abajo
    if (x2 == N - 2) {
        abajo2 = 1; // Bordes fijos en "x=N-1"
    } else {
        abajo2 = red[x2+1][y2];
    }

    // Vecino de la izquierda (condiciones periódicas en "y")
    if (y2 == 0) {
        izquierda2 = red[x2][N-1];
    } else {
        izquierda2 = red[x2][y2-1];
    }

    // Vecino de la derecha (condiciones periódicas en "y")
    if (y2 == N - 1) {
        derecha2 = red[x2][0];
    } else {
        derecha2 = red[x2][y2+1];
    }

    //No es necesario calcular las interacciones entre los espines que se van a intercambiar, ya que se cancelan entre sí en la diferencia de energía
    // Anulamos las contribuciones de (x2, y2) en la energía de (x1, y1)
    //Para ello vemos a qué vecino de (x1, y1) corresponde (x2, y2) y lo anulamos

    //Si y1=N-1, y1+1=N y al aplicar el modulo N, se obtiene 0
    //2 es el vecino de la derecha de 1

    if ((x1 == x2 && (y1 + 1) % N == y2)) {
        derecha1 = 0;
        izquierda2 = 0;
    }

    //Si y1=0, y1-1+N=N-1 y al aplicar el modulo N, se obtiene N-1
     //2 es el vecino de la izquierda de 1

    if ((x1 == x2 && (y1 - 1 + N) % N == y2)) {
        izquierda1 = 0;
        derecha2 = 0;
    }
    //2 es el vecino de abajo de 1
    if ((x1 + 1 == x2 && y1 == y2)) {
        arriba2 = 0;
        abajo1 = 0;
    }
    //2 es el vecino de arriba de 1
    if ((x1 - 1 == x2 && y1 == y2)) {
        abajo2 = 0;
        arriba1 = 0;
    }

    // Energía del primer espín
    E = -red[x1][y1] * (arriba1 + abajo1 + izquierda1 + derecha1);
    // Energía del segundo espín
    E += -red[x2][y2] * (arriba2 + abajo2 + izquierda2 + derecha2);

  //Al final solo realizamos 6 evaluaciones antes del intercambio y 6 después del intercambio. 

    return E;
}

// Función para imprimir la red
void imprimirRed(int red[N][N]) {
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            printf("%2d ", red[x][y]);
        }
        printf("\n");
    }
    printf("\n");
}

// Función para calcular la magnetización promedio en la mitad superior del sistema
double calcularMagnetizacion_Superior(int red[N][N]) {
    int suma = 0;
    for (int x = 0; x < N / 2; x++) { // Solo filas de la mitad superior
        for (int y = 0; y < N; y++) {
            suma += red[x][y];
        }
    }
    return (double) fabs(suma) / ((N / 2) * N);
}


// Función para calcular la magnetización promedio en la mitad inferior del sistema
double calcularMagnetizacion_Inferior(int red[N][N]) {
    int suma = 0;
    for (int x = N / 2; x < N; x++) { // Solo filas de la mitad inferior
        for (int y = 0; y < N; y++) {
            suma += red[x][y];
        }
    }
    return (double) fabs(suma) / ((N / 2) * N);
}

// Función para calcular la energía de la configuración
double calcularEnergia_Total_Inicial(int red[N][N]) {
    double energia = 0.0;

    for (int x = 1; x < N-1; x++) {
        for (int y = 0; y < N; y++) {
            // Suma de los vecinos
            int sumaVecinos = red[x-1][y] + red[x+1][y] +
                              red[x][(y + 1) % N] + red[x][(y - 1 + N) % N];

            // Contribución de la celda (x, y) a la energía
            energia += red[x][y] * sumaVecinos;
        }
    }

    // Multiplicar por -1/2 según la fórmula
    return -0.5 * energia;
}


// Función para calcular y mostrar las densidades promedio de espines positivos y negativos en una columna
void calcularDensidadesFila(int red[N][N], int fila, int *densidadpositivo, int *densidadnegativo) {
    int positivos = 0;
    int negativos = 0;
    // Recorrer la columna especificada
    for (int y= 0; y < N; y++) {
        if (red[(int)fila][y] == 1) {
            positivos++;
        } else if (red[(int)fila][y] == -1) {
            negativos++;
        }
    }

    // Calcular las densidades
   *densidadpositivo = positivos;
   *densidadnegativo = negativos;
}

void calcular_CalorEspecifico(double sumaEnergia, double sumaEnergiaCuadrada, int conteo, double *cv) {

      double promedio_E = sumaEnergia/ (conteo);
      double promedio_E2 = sumaEnergiaCuadrada /(conteo);
  
      double varianza_E = promedio_E2 - (promedio_E * promedio_E);
      // Calcular el calor específico a volumen constante
      *cv = varianza_E / ((N-2)*N*T*T);
}
void calcular_Susceptibilidad(double sumaMag, double sumaMag2, int conteo, double T, double *chi) {
    double promedio_M = sumaMag / conteo;
    double promedio_M2 = sumaMag2 / conteo;
    double varianza_M = promedio_M2 - (promedio_M * promedio_M);
    *chi = varianza_M / (N * N * T);
}

int algoritmoMetropolis(int (*red)[N], 
    double *sumaEnergia, double *sumaEnergiacuadrada, 
    double *sumaMagnetizacionSuperior, double *sumaMagnetizacionInferior, 
    int *configuracionesCambiadas, double *sumadensidadpositivo, 
    double *sumadensidadnegativo, double *sumaMagnetizacion2) {//double magnetizaciones[], int *indiceMag) {

    double energias[PASOS];
    double energia_actual;

    int densidadpositivo = 0.0;
    int densidadnegativo = 0.0;
    int fila=N/4;

    double energiaConfiguracion = 0.0;

    // Variables para evaluar la convergencia
    double sumaEnergiaIntervalo = 0.0;
    double sumaEnergiaCuadradaIntervalo = 0.0;
    int pasosIntervalo = 100; // Número de pasos en el intervalo para evaluar la convergencia

    
    // Abrir el fichero para guardar las configuraciones
    FILE *file = fopen("configuraciones.txt", "w");
    if (file == NULL) {
        printf("Error al abrir el fichero.\n");
        return 1;
    }

    
    FILE *energiaFile = fopen("energia_pmontecarlo.txt", "w");
    if (energiaFile == NULL) {
        printf("Error al abrir el fichero de energía.\n");
        return 1;
    }

    //Calculamos la energía incial de la configuración 
    energia_actual = calcularEnergia_Total_Inicial(red);

        for (int j = 0; j < PASOS; j++){
            for(int i=0; i<N*N; i++){
                // Elegir un espín al azar
                int y1 = rand() % N;
                int x1 = rand() % (N-2) +1 ; // Elegir un espín en el rango [1, N-2] para evitar los bordes fijos
                
                // Calcular un vecino al azar (arriba, abajo, derecha o izquierda)
                int deltaX = 0, deltaY = 0, dirección=0; 

                if (x1 == 1) {
                    // Descartar el intercambio en el borde fijo x=1
                    int direccion = rand() % 3; // Generar un número aleatorio entre 0 y 1 para elegir la dirección
                    if (direccion == 0) {
                        deltaY = 1;  // Vecino a la derecha
                    } else if (direccion == 1) {
                        deltaY = -1; // Vecino a la izquierda
                    } else {
                        deltaX = 1;  // Vecino abajo
                    }

                }else if (x1 == N - 2) {
                // Descartar el intercambio en el borde fijo y=N-1
                    int direccion = rand() % 3; // Generar un número aleatorio entre 0 y 1 para elegir la dirección
                    if (direccion == 0) {
                        deltaY= -1; // Vecino a la izquierda
                    } else if (direccion == 1) {
                        deltaY = 1;  // Vecino a la derecha
                    } else {
                        deltaX = -1; // Vecino arriba
                    }
                } else {
                    int direccion = rand() % 4; // Generar un número aleatorio entre 0 y 3 para elegir la dirección
                    if (direccion == 0) {
                        deltaY= 1;  // Vecino a la derecha
                    } else if (direccion == 1) {
                        deltaY = -1; // Vecino a la izquierda
                    } else if (direccion == 2) {
                        deltaX = 1;  // Vecino abajo
                    } else if (direccion == 3) {
                        deltaX = -1; // Vecino arriba
                    }
                }

                // Calcular las coordenadas del vecino
                int x2 = x1 + deltaX;
                int y2 = y1 + deltaY;

                // Aplicar condiciones de contorno periódicas en "y"
                if (y2 < 0) {
                    y2 = N - 1; // Si se sale por la izquierda, vuelve al borde derecho
                } else if (y2 >= N) {
                    y2 = 0; // Si se sale por la derecha, vuelve al borde izquierdo
                }

                // Asegurar que los espines tengan distinto signo para intercambiarlos
                if (red[x1][y1] == red[x2][y2]) {
                    continue; // Si tienen el mismo signo, pasar a la siguiente iteración
                }

                // Calcular la energía antes del intercambio
                int energiaAntes = calcularEnergia_Local(red, x1, y1, x2, y2);

                // Intercambiar los espines
                int temp = red[x1][y1];
                red[x1][y1] = red[x2][y2];
                red[x2][y2] = temp;
            
            
                // Calcular la energía después del intercambio
                int energiaDespues = calcularEnergia_Local(red, x1, y1, x2, y2);

                // Calcular la diferencia de energía
                int difE = energiaDespues - energiaAntes;

                // Calcular la probabilidad de transición
                double probabilidad;

                if (difE > 0) {
                    probabilidad = exp(-difE / T);
                } else {
                    probabilidad = 1.0;
                }

                // Generar un número aleatorio entre 0 y 1
                double r = (double)rand() / RAND_MAX;

                if(r > probabilidad) { 
                    // Revertir el intercambio si no se acepta
                    temp = red[x1][y1];
                    red[x1][y1] = red[x2][y2];
                    red[x2][y2] = temp;
                }else {
                    // Si se acepta, actualizar la energía actual con la diferencia de energía
                    energia_actual += difE;
                   // break;
                     // Salir del bucle para pasar al siguiente paso Monte Carlo
                }
            }

            // Guardar la energía en el array
            energias[j] = energia_actual;

            fprintf(energiaFile, "%.6f %d\n", energia_actual, j + 1);

            escribirRed(file, red);

                //PARA EVALUAR LA CONVERGENCIA

                sumaEnergiaIntervalo += energia_actual;
                sumaEnergiaCuadradaIntervalo += energia_actual * energia_actual;

            // Evaluar la convergencia cada "pasosIntervalo" pasos
            if ((j + 1) % pasosIntervalo == 0) {
                double promedioEnergia = sumaEnergiaIntervalo / pasosIntervalo;
                double promedioEnergiaCuadrada = sumaEnergiaCuadradaIntervalo / pasosIntervalo;
                double varianzaEnergia = promedioEnergiaCuadrada - (promedioEnergia * promedioEnergia);
                double desviacionEstandarEnergia = sqrt(varianzaEnergia);

                // Reiniciar acumuladores para el siguiente intervalo
                sumaEnergiaIntervalo = 0.0;
                sumaEnergiaCuadradaIntervalo = 0.0;

                // Verificar si la varianza es suficientemente pequeña
                 if (desviacionEstandarEnergia < 30) {
                    printf("Convergencia alcanzada en el paso %d.\n", j + 1);
                   break; // Salir del bucle si se alcanza la convergencia
                    }
            }
        
            // PARA CALCULAR LA MAGNETIZACIÓN Y DENSIDAD CADA 100 PASOS MONTECARLO
            if ((j + 1) % (100) == 0) {

                double magnetizacionSuperior = calcularMagnetizacion_Superior(red);
                double magnetizacionInferior = calcularMagnetizacion_Inferior(red);

                //Sumamos para después hacer el promedio 
    
                *sumaMagnetizacionSuperior += magnetizacionSuperior;
                *sumaMagnetizacionInferior += magnetizacionInferior;
            
                calcularDensidadesFila(red, fila, &densidadpositivo, &densidadnegativo);

                *sumadensidadpositivo += densidadpositivo;
                *sumadensidadnegativo += densidadnegativo;
                
                energiaConfiguracion = energias[j];

                *sumaEnergia += energiaConfiguracion;

                //Necesaria para calcular el calor específico
                *sumaEnergiacuadrada += energiaConfiguracion * energiaConfiguracion;

                //Necesaria para calcular la susceptibilidad
                *sumaMagnetizacion2 += magnetizacionSuperior * magnetizacionSuperior;
                
                (*configuracionesCambiadas)++;

            }


        }

        fclose(file);
        fclose(energiaFile);
}

int main() {

   // Registrar el tiempo de inicio
   clock_t start = clock();


    int red[N][N];

    srand(time(NULL)); // Semilla para números aleatorios

    // Preguntar al usuario qué tipo de magnetización inicial desea
    int opcion;
    printf("Seleccione la condición inicial:\n");
    printf("1. Magnetización inicial nula (desordenada aleatoriamente)\n");
    printf("2. Mitad positivos, mitad negativos (aleatorio)\n");
    //printf("Opción: ");
    //scanf("%d", &opcion);
    opcion=1;

    // Inicializar la red según la elección del usuario
    if (opcion == 1) {
    inicializarRed_aleatoria(red);
    } else if (opcion == 2) {
    inicializarRed_mitad(red);
    } else {
    printf("Opción no válida. Saliendo del programa.\n");
    return 1;
    }

    // Abrir los ficheros para guardar las magnetizaciones
    FILE *magnetizacionSuperiorFile = fopen("promedio_magnetizacionsuperior.txt", "a");
    if (magnetizacionSuperiorFile == NULL) {
        printf("Error al abrir el fichero de magnetización superior.\n");
        return 1;
    }

    FILE *magnetizacionInferiorFile = fopen("promedio_magnetizacioninferior.txt", "a");
    if (magnetizacionInferiorFile == NULL) {
        printf("Error al abrir el fichero de magnetización inferior.\n");
        return 1;
    }

     FILE *filedensidadpositivo = fopen("promedio_densidadpositivo.txt", "a");
    if (filedensidadpositivo == NULL) {
        return 1;
    }
    FILE *filedensidadnegativo = fopen("promedio_densidadnegativo.txt", "a");
    if (filedensidadnegativo == NULL) {
        return 1;
    }
    FILE *filecv = fopen("filecv.txt", "a");
    if (filecv == NULL) {
        return 1;
    }

    FILE *filechi = fopen("susceptibilidad.txt", "a");
    if (filechi == NULL) {
        return 1;
    }

    // Imprimir la configuración inicial
    //printf("Configuración inicial de la red:\n");
    //imprimirRed(red);

      // Variables para acumular magnetización
    double sumaMagnetizacionSuperior = 0.0;
    double sumaMagnetizacionInferior = 0.0;

    int configuracionesCambiadas = 0;

    double sumaEnergia = 0.0;
    double sumaEnergiacuadrada = 0.0;

    double sumadensidadpositivo = 0.0;
    double sumadensidadnegativo = 0.0;
    
    double cv;
    double chi;

    double sumaMagnetizacion2 = 0.0;

    algoritmoMetropolis(red, &sumaEnergia, &sumaEnergiacuadrada, 
        &sumaMagnetizacionSuperior, &sumaMagnetizacionInferior, 
        &configuracionesCambiadas, &sumadensidadpositivo, 
        &sumadensidadnegativo, &sumaMagnetizacion2 ); 


        //Calculamos los promedios de magnetización y densidad
    
        double promedioMagnetizacionSuperior = sumaMagnetizacionSuperior / configuracionesCambiadas;
        double promedioMagnetizacionInferior = sumaMagnetizacionInferior / configuracionesCambiadas;

        double promedioDensidadPositiva = sumadensidadpositivo / (configuracionesCambiadas*N);
        double promedioDensidadNegativa = sumadensidadnegativo / (configuracionesCambiadas*N);


        fprintf(magnetizacionSuperiorFile, "%.6f %.1f\n", promedioMagnetizacionSuperior, T);
        fprintf(magnetizacionInferiorFile, "%.6f %.1f\n", promedioMagnetizacionInferior, T);
        

        fprintf(filedensidadnegativo, "%.6f %.1f\n", promedioDensidadNegativa, T);
        fprintf(filedensidadpositivo, "%.6f %.1f\n", promedioDensidadPositiva, T);


         // Calcular y mostrar la energía media por partícula
          double energiaMediaPorParticula = sumaEnergia / (configuracionesCambiadas * 2 * N * N);
            printf("Energía media por partícula: %.6f\n", energiaMediaPorParticula);
        
        // Calcular y mostrar el calor específico
        calcular_CalorEspecifico(sumaEnergia, sumaEnergiacuadrada, configuracionesCambiadas, &cv);
        fprintf(filecv, "%.6f %.1f\n", cv, T);

        //Calcular y mostrar la susceptibilidad
        calcular_Susceptibilidad(sumaMagnetizacionSuperior, sumaMagnetizacion2, configuracionesCambiadas, T, &chi);
        fprintf(filechi, "%.10f %.1f\n", chi, T);

        fclose(magnetizacionSuperiorFile);
        fclose(magnetizacionInferiorFile);
        fclose(filedensidadpositivo);
        fclose(filedensidadnegativo);
        fclose(filecv);
        fclose(filechi);

         // Registrar el tiempo de finalización
    clock_t end = clock();

    
    // Calcular el tiempo de CPU usado
    double cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecución: %.6f segundos\n", cpu_time_used);

     return 0;
}

