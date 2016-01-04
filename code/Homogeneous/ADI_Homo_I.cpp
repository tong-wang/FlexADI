//============================================================================
// Name        : ADI_Homo_I.cpp
// Author      : Tong WANG
// Email       : tong.wang@nus.edu.sg
// Version     : v2.0 (2015-12-16)
// Copyright   : MIT License
// Description : code for ADI Model 1 (Homogeneous), case I (T <= L + 1)
//============================================================================
// Compile: [GCC] g++ -O3 -o ADI_Homo_I.exe ADI_Homo_I.cpp
//          [ICC] icpc -O3 -o ADI_Homo_I.exe ADI_Homo_I.cpp
//---------------------------------------------------------------------------

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> //required by setprecision()
#include <cmath>

using namespace std;

//---------------------------------------------------------------------------

#define N 30            // number of periods

#define DMAX 30         // upper bound of demand realization
#define YMAX 80
#define UMAX 1100
#define UMID 1000

//---------------------------------------------------------------------------
int     L, T;           // supply and demand leadtime
double	K, h, b;        // fixed cost, holding cost, backlog cost
double	lambda;         // Poisson demand rate

double  phi[DMAX + 1];  // pre-calculated probability distribution for single-period demand D
double	phi_I[5 * DMAX + 1];  // pre-calculated probability distribution for leadtime demand D1 (at least (L-T+1)*DMAX)
double	f[N + 1][UMAX + 1];  // for saving the DP value functions
int		Y[UMAX + 1];    // for saving Period 1 optimal (state-independent) inventory positions
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
		for	(int d_I = 0; d_I <= (L - T + 1) * DMAX; d_I++)
			phi_I[d_I] = Poisson_PMF((L - T + 1) * lambda, d_I);


		for	(int d = 0; d <= DMAX; d++)
			phi[d] = Poisson_PMF(lambda, d);
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
	
		for	(int d_I = 0; d_I <= (L - T + 1) * DMAX; d_I++)
			g += l(yy - d_I) * phi_I[d_I];
							
		
		if (ii < N - L)
            for	(int d0 = 0; d0 <= DMAX; d0++)
				g += F(ii + 1, yy - d0) * phi[d0];

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
    file.open("Result_ADI_Homo_I.txt", fstream::app|fstream::out);
    
    if (! file)
    {
        //if fail to open the file
        cerr << "can't open output file Result_ADI_Homo_I.txt!" << endl;
        exit(EXIT_FAILURE);
    }
    
    file << setprecision(10);
    cout << setprecision(10);

    //--------------------------------------------------------------
    // output header to screen and file
    cout << "N\tL\tT\tlambda\tK\th\tb\tCost\ts0\tS0" << endl;
    file << "N\tL\tT\tlambda\tK\th\tb\tCost\ts0\tS0" << endl;
    

    // initialize parameters
    lambda = 6;
    K = 100;
    h = 1.0;
    b = 9.0;
   
    
    for (L = 0; L <= 4; L++)
        for (T = 0; T <= min(L + 1, 2); T++)
        {
            cout << N << "\t" << L << "\t" << T << "\t" << lambda << "\t" << K << "\t" << h << "\t" << b << "\t";
            file << N << "\t" << L << "\t" << T << "\t" << lambda << "\t" << K << "\t" << h << "\t" << b << "\t";

            
            // initialize demand distributions
            phi_Init();


            // initialize memory for saving the value function
            for (int i = 0; i <= N; i++)
                for (int j = 0; j <= UMAX; j++)
                    f[i][j] = -7;
            
            // uncontrollable costs
            double C1 = 0;

            for (int j = T + 1; j <= L; j++)
            {
                double C_j = 0;

                for	(int d = 0; d <= (j - T) * DMAX; d++)
                    C_j += l(0 - d) * Poisson_PMF((j - T) * lambda, d);
                
                C1 += C_j;
            }

            // cost via the DP
            C1 += F(1, 0);
            
            cout << C1 << "\t";
            file << C1 << "\t";

            
            //search for optimal (s, S) in period 1 by bisectional search
            int s_l = -40, s_u = 40, s;

            for (; s_u - s_l > 2; )
            {
                s = (s_u + s_l) / 2;
                F(1, s);
                //cout << "\n" << s_l << "\t" << s_u << "\t" << s << "\t" << Y[s + UMID] << endl;

                if (Y[s + UMID] > s)
                    s_l = s;
                else
                    s_u = s;
            }

            for (s = s_u; s >= s_l; s--)
            {
                F(1, s);
                //cout << "\n" << s_l << "\t" << s_u << "\t" << s << "\t" << Y[s + UMID] << endl;

                if (Y[s + UMID] > s)
                {
                    s_Opt = s;
                    S_Opt = Y[s + UMID];
                    
                    cout << s_Opt << "\t" << S_Opt << endl;
                    file << s_Opt << "\t" << S_Opt << endl;

                    break;
                }
            }
        }


    file.close();
    return 0;
}