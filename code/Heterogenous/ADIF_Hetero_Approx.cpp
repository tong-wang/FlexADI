//============================================================================
// Name        : ADIF_Hetero_Approx.cpp
// Author      : Tong WANG
// Email       : tong.wang@nus.edu.sg
// Version     : v2.0 (2015-12-16)
// Copyright   : MIT License
// Description : code for ADI Model 2 (Heterogeneous) with flexible delivery, case II (T > L + 1)
//              --- Lower bound Approximation (assume delivery can be called back to satisfy urgent demand)
//============================================================================
// Compile: [GCC] g++ -O3 -o ADIF_Hetero_Approx.exe ADIF_Hetero_Approx.cpp
//          [ICC] icpc -O3 -o ADIF_Hetero_Approx.exe ADIF_Hetero_Approx.cpp
//---------------------------------------------------------------------------

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> //required by setprecision()
#include <cmath>

using namespace std;

//---------------------------------------------------------------------------

#define N 12            // number of periods
#define L 0             // supply leadtime (fixed at 0)
#define T 2             // demand leadtime (fixed at 2)

#define DMAX 20         // upper bound of demand realization
#define YMAX 80
#define UMAX 1100
#define UMID 1000
#define	VMAX 20

//---------------------------------------------------------------------------
double	K, h, b;        // fixed cost, holding cost, backlog cost
double	lambda[T + 1];         // Poisson demand rate

double  phi[T + 1][DMAX + 1];  // pre-calculated probability distribution for demand D 0, 1, 2
double	f[N + 1][UMAX + 1][VMAX + 1];  // for saving the DP value functions
int		Y[UMAX + 1][VMAX + 1];    // for saving Period 1 optimal (state-dependent) inventory positions

int		S_Opt[VMAX + 1], s_Opt[VMAX + 1];   // optimal state-dependent (s, S) policy
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
}


// single-period cost function
double l(int xx)
{
    return (xx > 0) ? (h * xx) : (-b * xx);
}

// on-hand inventory as a function of inventory position and advance demand
int X_AP_II(int yy, int vv, int d0, int d1, int d2)
{
    if (yy - d0 - d1 - d2 >= 0) return yy - d0 - d1 - d2;
    else if (yy + vv - d0 <= 0) return yy + vv - d0;
    else return 0;
}

// DP value functions
double F(int ii, int uu, int vv);

double G(int ii, int yy, int vv)
{
    double g = 0;
    
    for	(int d0 = 0; d0 <= DMAX; d0++)
        for (int d1 = 0; d1 <= DMAX; d1++)
            for (int d2 = 0; d2 <= DMAX; d2++)
                g += (l(X_AP_II(yy, vv, d0, d1, d2)) + F(ii + 1, yy - d0 - d1 - d2, d2)) * phi[0][d0] * phi[1][d1] * phi[2][d2];
    
    return g;
}

double F(int ii, int uu, int vv)
{
    if (ii >= N + 1)
        return 0;
    // if the corresponding value function has been calculated before, load it directly from the memory
    else if ((uu + UMID >= 0) && (uu + UMID <= UMAX) && (f[ii][uu + UMID][vv] != -7))
    {
        return f[ii][uu + UMID][vv];
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
            
            temp += G(ii, yy, vv);
            
            if (temp <= min)
            {
                min = temp;
                y_opt = yy;
            }
            
            if (temp > min + K)
                break;
        }
        
        // save the newly calculated value function into the memeory
        f[ii][uu + UMID][vv] = min;
        
        // record optimal policy if in period 1
        if (ii == 1)
            Y[uu + UMID][vv] = y_opt;

        return min;
    }
}



int main()
{
    
    //Open output file
    file.open("Result_ADIF_Hetero_Approx.txt", fstream::app|fstream::out);
    
    if (! file)
    {
        //if fail to open the file
        cerr << "can't open output file Result_ADIF_Hetero_Approx.txt!" << endl;
        exit(EXIT_FAILURE);
    }
    
    file << setprecision(10);
    cout << setprecision(10);
    
    //--------------------------------------------------------------
    // output header to screen and file
    cout << "N\tL\tT\tlambda0\tlambda1\tlambda2\tK\th\tb\tCost_AP";
    file << "N\tL\tT\tlambda0\tlambda1\tlambda2\tK\th\tb\tCost_AP";
    
    for (int vv = 0; vv <= VMAX; vv++)
    {
        cout << "\t" << "s" << vv << "_AP\t" << "S" << vv << "_AP";
        file << "\t" << "s" << vv << "_AP\t" << "S" << vv << "_AP";
    }
    cout << endl;
    file << endl;
    
    
    // initialize parameters
    //K = 100;
    //h = 1.0;
    //b = 9.0;
    //lambda[0] = 4;
    //lambda[1] = 1;
    //lambda[2] = 1;
    
    for (lambda[0] = 0; lambda[0] <= 5; lambda[0]++)
    {
        lambda[1] = 1;
        lambda[2] = 5 - lambda[0];
        
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
                    for (int k = 0; k <= VMAX; k++)
                        f[i][j][k] = -7;
            
            // cost via the DP
            double C1 = F(1, 0, 0);

            cout << C1 << "\t";
            file << C1 << "\t";

            
            //search for optimal (s, S) in period 1
            for (int j = -40; j <= 40; j++)
                for (int k = 0; k <= VMAX; k++)
                    F(1, j, k);

            for (int vv = 0; vv <= VMAX; vv++)
            {
                for (int uu = UMAX; uu >= 0; uu--)
                {
                    if (Y[uu][vv] > uu - UMID)
                    {
                        S_Opt[vv] = Y[uu][vv];
                        s_Opt[vv] = uu - UMID;
                        
                        cout << s_Opt[vv] << "\t" << S_Opt[vv] << "\t";
                        file << s_Opt[vv] << "\t" << S_Opt[vv] << "\t";
                        
                        break;
                    }
                }
                
            }

            cout << endl;
            file << endl;
        }
    }
    
    file.close();
    return 0;
}