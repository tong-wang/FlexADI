//============================================================================
// Name        : ADI_Hetero_I.cpp
// Author      : Tong WANG
// Email       : tong.wang@nus.edu.sg
// Version     : v2.0 (2015-12-16)
// Copyright   : MIT License
// Description : code for ADI Model 2 (Heterogeneous), case I (T <= L + 1)
//============================================================================
// Compile: [GCC] g++ -O3 -o ADI_Hetero_I.exe ADI_Hetero_I.cpp
//          [ICC] icpc -O3 -o ADI_Hetero_I.exe ADI_Hetero_I.cpp
//---------------------------------------------------------------------------

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> //required by setprecision()
#include <cmath>

using namespace std;

//---------------------------------------------------------------------------

#define N 12            // number of periods
#define T 2             // demand leadtime (fixed at 2)

#define DMAX 20         // upper bound of demand realization
#define YMAX 80
#define UMAX 1100
#define UMID 1000

//---------------------------------------------------------------------------
int L, I;
double	K, h, b;        // fixed cost, holding cost, backlog cost
double	lambda[T + 1];         // Poisson demand rate
double  lambda_I;

double  phi[T + 1][DMAX + 1];  // pre-calculated probability distribution for demand D 0, 1, 2
double  phi_I[12 * DMAX +1], phi_012[(T + 1) * DMAX + 1];  //D_I  is up to 3*L*DMAX
double	f[N + 1][UMAX + 1];  // for saving the DP value functions
int		Y[UMAX + 1];    // for saving Period 1 optimal (state-dependent) inventory positions

int		S_Opt, s_Opt;   // optimal (s, S) policy
ofstream file;          // output file

//============================================================================
// Probability Mass Function of Poisson Distribution
double Poisson_PMF(double lambda, int k)
{
    if ((k < 0) | (lambda < 0))
        return 0;
    else if (lambda == 0)
        return (k == 0) ? 1 : 0;
    else
    {
        double logP = - lambda + k * log(lambda);
        
        for (int i = k; i > 1; i--)
            logP -= log(i);
        
        return exp(logP);
    }
}

// pre-calculate probabilities
void phi_Init()
{
    for	(int t = 0; t <= T; t++)
        for	(int d = 0; d <= DMAX; d++)
            phi[t][d] = Poisson_PMF(lambda[t], d);

    for (int t = 0; t <= I * DMAX; t++)
        phi_I[t] = Poisson_PMF(lambda_I, t);
    
    for (int t = 0; t <= 3 * DMAX; t++)
        phi_012[t] = Poisson_PMF(lambda[0] + lambda[1] + lambda[2], t);
}


// single-period cost function
double l(int xx)
{
    return (xx > 0) ? (h * xx) : (-b * xx);
}


// DP value functions
double F(int ii, int uu);

double G(int ii, int yy)
{
    double g = 0;
    
    for (int d_I = 0; d_I <= I * DMAX; d_I++)
        g += l(yy - d_I) * phi_I[d_I];
    
    for (int d_012 = 0; d_012 <= (T + 1) * DMAX; d_012++)
        g += F(ii + 1, yy - d_012) * phi_012[d_012];
    
    return g;
}

double F(int ii, int uu)
{
    if (ii >= N - L + 1)
        return 0;
    // if the corresponding value function has been calculated before, load it directly from the memory
    else if ((uu + UMID >= 0) && (uu + UMID <= UMAX) && (f[ii][uu + UMID] != -7))
    {
        return f[ii][uu + UMID];
    }
    // if not, need to calculate by recursion
    else
    {
        double	min = 1.0e+30;
        int		y_opt = 0;
        
        for (int yy = uu; yy <= YMAX; yy++)
        {
            double temp = 0;
            
            if (yy > uu)
                temp += K;
            
            temp += G(ii, yy);
            
            if (temp <= min)
            {
                min = temp;
                y_opt = yy;
            }
            
            if (temp > min + K)
                break;
        }
        
        // save the newly calculated value function into the memeory
        f[ii][uu + UMID] = min;
        
        // record optimal policy if in period 1
        if (ii == 1)
            Y[uu + UMID] = y_opt;
        
        return min;
    }
}



int main()
{
    
    //Open output file
    file.open("Result_ADI_Hetero_I.txt", fstream::app|fstream::out);
    
    if (! file)
    {
        //if fail to open the file
        cerr << "can't open output file Result_ADI_Hetero_I.txt!" << endl;
        exit(EXIT_FAILURE);
    }
    
    file << setprecision(10);
    cout << setprecision(10);
    
    //--------------------------------------------------------------
    // output header to screen and file
    cout << "N\tL\tT\tlambda0\tlambda1\tlambda2\tK\th\tb\tCost\ts\tS" << endl;
    file << "N\tL\tT\tlambda0\tlambda1\tlambda2\tK\th\tb\tCost\ts\tS" << endl;
    
    
    // initialize parameters
    //K = 100;
    //h = 1.0;
    //b = 9.0;
    //lambda[0] = 4;
    //lambda[1] = 1;
    //lambda[2] = 1;
    
    for(L = 1; L <= 4; L++)
    for (lambda[0] = 0; lambda[0] <= 5; lambda[0]++)
    {
        lambda[1] = 1;
        lambda[2] = 5 - lambda[0];
        
        I = 3 * L;
        lambda_I = (L + 1) * lambda[0] + L * lambda[1] + (L - 1) * lambda[2];

        // initialize demand distributions
        phi_Init();
        
        for(K = 50; K <= 200; K *= 2)
        for(h = 1; h <= 5; h += 2)
        for(b = 9; b <= 29; b += 10)
        {
            
            cout << N << "\t" << L << "\t" << T << "\t" << lambda[0] << "\t" << lambda[1] << "\t" << lambda[2] << "\t" << K << "\t" << h << "\t" << b << "\t";
            file << N << "\t" << L << "\t" << T << "\t" << lambda[0] << "\t" << lambda[1] << "\t" << lambda[2] << "\t" << K << "\t" << h << "\t" << b << "\t";
            
            // initialize memory for saving the value function
            for (int i = 0; i <= N; i++)
                for (int j = 0; j <= UMAX; j++)
                        f[i][j] = -7;
            
            // cost via the DP
            //calculate optimal cost in period 1
            double C1 = 0;
            
            //uncontrollable
            double C_unCTRL[5] = {0, 0, 0, 0, 0};
            
            for (int d1 = 0; d1 <= DMAX; d1++)
                C_unCTRL[1] += l(- d1) * phi[0][d1];
            for (int d2 = 0; d2 <= 3 * DMAX; d2++)
                C_unCTRL[2] += l(- d2) * Poisson_PMF(2 * lambda[0] + lambda[1], d2);
            for (int d3 = 0; d3 <= 6 * DMAX; d3++)
                C_unCTRL[3] += l(- d3) * Poisson_PMF(3 * lambda[0] + 2 * lambda[1] + lambda[2], d3);
            for (int d4 = 0; d4 <= 9 * DMAX; d4++)
                C_unCTRL[4] += l(- d4) * Poisson_PMF(4 * lambda[0] + 3 * lambda[1] + 2 * lambda[2], d4);
            
            for (int l = 1; l <= L; l++)
                C1 += C_unCTRL[l];

            C1 += F(1, 0);
            
            cout << C1 << "\t";
            file << C1 << "\t";
            
            
            //search for optimal (s, S) in period 1
            for (int j = -40; j <= 40; j++)
                    F(1, j);
            
            for (int uu = UMAX; uu >= 0; uu--)
            {
                if (Y[uu] > uu - UMID)
                {
                    S_Opt = Y[uu];
                    s_Opt = uu - UMID;
                    
                    cout << s_Opt << "\t" << S_Opt << "\t";
                    file << s_Opt << "\t" << S_Opt << "\t";
                    
                    break;
                }
            }
            
            cout << endl;
            file << endl;
            
        }
    }

    file.close();
    return 0;
}