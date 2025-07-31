#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define N 25
#define tau 125
#define NRUN 10000   
//#define M_PI 3.14159265359
float **J;
//#define RELAX 5000 // Relaxation time steps

// Function to generate Gaussian-distributed random numbers using Box-Muller method
double gaussian_random() {
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return z0;
}

// Function to initialize the J_ij matrix with Gaussian values
void initialize_J() {
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            if (i == j) {
                J[i][j] = 0.0;  // No self-interaction
            } else {
                double Jij = gaussian_random() / sqrt(N);  // Variance 1/N
                J[i][j] = Jij;
                J[j][i] = Jij;  // Ensure symmetry: J_ij = J_ji
            }
        }
    }
}

void main() {
    clock_t start, end; 
    start = clock(); 
   
    int TMAX = tau,ii,jj,i,j;
    double mz[N],mx[N],md,md1,gamma,gamma0,ran, h,h0, T,T0, changex, changez, del, heff, sum, q2, engg, engg1,ent[TMAX],ent1[TMAX];
    srand(time(NULL));  // Seed for random numbers

    // Allocate memory for J matrix
    J = (float**)malloc(N * sizeof(float*));
    if (J == NULL) {
        printf("Memory allocation failed for J matrix.\n");
        return;
    }
    for (int i = 0; i < N; i++) {
        J[i] = (float*)malloc(N * sizeof(float));
        if (J[i] == NULL) {
            printf("Memory allocation failed for J[%d].\n", i);
            for (int k = 0; k <= i; k++) free(J[k]);
            free(J);
            return;
        }
    }
  
    // Initialize parameters
  //  h0 =0.;
    h=0;
    T0 = 0.;
    gamma0= 1.0;
    //T0 = 0.;  
    //del = 0.0001;

   // for(i=0;i<TMAX+1;i++) {ent[i]=0; ent1[i]=0;}

    for (int n = 0; n < NRUN; n++) {
        // Initialize spins
        for (int i = 0; i < N; i++) {
            mz[i] = (rand() / (double)RAND_MAX < 0.5) ? 1.0 : -1.0;
            mx[i] = (rand() / (double)RAND_MAX < 0.5) ? 0.0 :  0.0;
        }
   
        initialize_J();  // SK model initialization

        for (int t = 0; t < TMAX+1; t++) {
             T = T0 * (1.0 - (float)t / tau);
             //if(t>TMAX+1){gamma=(rand() / (double)RAND_MAX < 0.) ? 0.001 : 0.0;
             gamma=gamma0*(1.0 - (float)t/(tau));
             
             //else{gamma=gamma0*(1.0 - (float)t/tau);}
            // ran=rand() / (double)RAND_MAX;
            // h=-0.1+0.2*ran;
            // h=h0*(1.0 - (float)t/tau);
             
            for (int ii = 0; ii < N; ii++) {
                int i = rand() % N;  // Random index in range [0, N-1]
                sum = 0;
                q2 = 0;

                for (int j = 0; j < N; j++) {
                    if (i != j){ sum += J[i][j] * mz[j];}
                    q2 += mz[j] * mz[j];
                }
                q2 = q2 / N;

                heff = sqrt((sum + h - (1 - q2) * mz[i]) * (sum + h - (1 - q2) * mz[i]) + gamma * gamma);
              // heff = sqrt((sum + h - (1 - q2) * mz[i]) * (sum + h - (1 - q2) * mz[i]) + (gamma-(1-q2)*mx[i])*(gamma-(1-q2)*mx[i])); 
                changez = tanh(heff / T) * (sum + h - (1 - q2) * mz[i]) / heff;
                changex=tanh(heff/T)*(gamma)/heff;
               // changex=tanh(heff/T)*(gamma-(1-q2)*mx[i])/heff;
                
        
                if (changez * changez > 0) {
                  //  mx[i]=changex;
                    mz[i] = changez;
                }
                if (changex * changex > 0) {
                  //  mx[i]=changex;
                    mx[i] = changex;
                }
                
            }

            sum = 0; 
            engg = 0;engg1=0;
            for(ii=0;ii<N;ii++){
       sum+=mz[ii];
       q2+=mz[ii]*mz[ii];
      
        md=0; md1=0;
        for(jj=0;jj<N;jj++)
        {
          if(ii!=jj)
          {
            md+=J[ii][jj]*mz[jj];
            if(mz[jj]>0) {md1+=J[ii][jj];}
            else {md1+=-1*J[ii][jj];}
          }
        }
        md=md*mz[ii];

        if(mz[ii]<=0) {md1=md1*(-1);}  
        engg+=md;
        engg1+=md1;
      }
      
     // ent[t]+=-1*engg/(2*N);
     // ent1[t]+=-1*engg1/(2*N);
            if (t == TMAX) {
                printf("%d %f\n", t, -1*engg/(2 * N));
            }
            
        }
    } //nrun loop

  /* for(i=0;i<TMAX+1;i++){
       printf("%d %lf \n",i,ent[i]/NRUN);
   }*/

    // Free allocated memory
  /*  for (int i = 0; i < N; i++) {
        free(J[i]);
    } 
    free(J);*/

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC; // in seconds 
    printf("#Execution took %f seconds.\n", time_taken); 
    printf("#Config. av. %d\n",NRUN);
}

