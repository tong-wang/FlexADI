//============================================================================
// Name        : ADIF_Homo_II.cpp
// Author      : Tong WANG
// Email       : tong.wang@nus.edu.sg
// Version     : v2.0 (2015-12-16)
// Copyright   : MIT License
// Description : code for ADI Model 1 (Homogeneous) with flexible delivery, case II (T > L + 1), state-dependent
//============================================================================
// Compile: [GCC] g++ -O3 -o ADIF_Homo_II.exe ADIF_Homo_II.cpp
//          [ICC] icpc -O3 -o ADIF_Homo_II.exe ADIF_Homo_II.cpp
//---------------------------------------------------------------------------

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> //required by setprecision()
#include <cmath>

using namespace std;

//---------------------------------------------------------------------------

#define N 30            // number of periods
#define L 0             // supply leadtime (fixed at 0)
#define T 2             // demand leadtime (fixed at 2)

#define DMAX 30         // upper bound of demand realization
#define YMAX 80
#define UMAX 1100
#define UMID 1000
#define	VMAX 30

//---------------------------------------------------------------------------
double	K, h, b;        // fixed cost, holding cost, backlog cost
double	lambda;         // Poisson demand rate

double  phi[DMAX + 1];  // pre-calculated probability distribution for single-period demand D
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
    for	(int d = 0; d <= DMAX; d++)
        phi[d] = Poisson_PMF(lambda, d);
}


// single-period cost function
double l(int xx)
{
    return (xx > 0) ? (h * xx) : (-b * xx);
}

// on-hand inventory as a function of inventory position, state, and advance demand
int X(int yy, int vv, int dd)
{
    if (yy - dd > 0) return yy - dd;
    else if (yy + vv < 0) return yy + vv;
    else return 0;
}

// DP value functions
double F(int ii, int uu, int vv);

double G(int ii, int yy, int vv)
{
    double g = 0;
    
    for	(int d = 0; d <= DMAX; d++)
        g += (l(X(yy, vv, d)) + F(ii + 1, yy - d, d)) * phi[d];
    
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
    file.open("Result_ADIF_Homo_II.txt", fstream::app|fstream::out);
    
    if (! file)
    {
        //if fail to open the file
        cerr << "can't open output file Result_ADIF_Homo_II.txt!" << endl;
        exit(EXIT_FAILURE);
    }
    
    file << setprecision(10);
    cout << setprecision(10);
    
    //--------------------------------------------------------------
    // output header to screen and file
    cout << "N\tL\tT\tlambda\tK\th\tb\tCost";
    file << "N\tL\tT\tlambda\tK\th\tb\tCost";
    
    for (int vv = 0; vv <= VMAX; vv++)
    {
        cout << "\t" << "s" << vv << "\t" << "S" << vv;
        file << "\t" << "s" << vv << "\t" << "S" << vv;
    }
    cout << endl;
    file << endl;

    
    // initialize parameters
    lambda = 6;
    K = 100;
    h = 1.0;
    b = 9.0;
    
    //for(K = 50; K <= 200; K *= 2)
    //for(h = 1; h <= 5; h += 2)
    //for(b = 9; b <= 29; b += 10)
    {
        
        cout << N << "\t" << L << "\t" << T << "\t" << lambda << "\t" << K << "\t" << h << "\t" << b << "\t";
        file << N << "\t" << L << "\t" << T << "\t" << lambda << "\t" << K << "\t" << h << "\t" << b << "\t";
        

        // initialize demand distributions
        phi_Init();
        
        
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
    
    file.close();
    return 0;
}