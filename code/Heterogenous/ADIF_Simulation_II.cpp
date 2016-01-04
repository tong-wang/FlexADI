//============================================================================
// Name        : ADIF_Simulation_II.cpp
// Author      : Tong WANG
// Email       : tong.wang@nus.edu.sg
// Version     : v2.0 (2015-12-16)
// Copyright   : MIT License
// Description : Simulation of the heuristics for the ADI model with flexible delivery, case II (T > L + 1)
//============================================================================
// Compile: [GCC] g++ -O3 -std=c++11 -o ADIF_Simulation_II.exe ADIF_Simulation_II.cpp
//          [ICC] icpc -O3 -std=c++11 -o ADIF_Simulation_II.exe ADIF_Simulation_II.cpp
//---------------------------------------------------------------------------

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> //required by setprecision()
#include <cmath>
#include <random>

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
double	lambda[T + 1];          // Poisson demand rate

double  phi[T + 1][DMAX + 1];   // pre-calculated probability distribution for demand D 0, 1, 2
double	f[N + 1][UMAX + 1][VMAX + 1];  // for saving the DP value functions
int		Y[N + 1][UMAX + 1][VMAX + 1];    // for saving Period 1 optimal (state-dependent) inventory positions

int		S_Opt[N+1][VMAX + 1], s_Opt[N+1][VMAX + 1];   // optimal state-dependent (s, S) policy in all the periods
int		Sigma, sigma[N+1];      //protection levels for PL(Sigma) and PL(sigma)

std::mt19937 random_engine(654321);     //std::mt19937 random engine;
double	C_AP, C_PL_0, C_PL_s, C_PL_S;			//simulated costs of different strategies

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

// approximation: on-hand inventory as a function of DP states and demand (case II)
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
        
        // record optimal policy
        Y[ii][uu + UMID][vv] = y_opt;
        
        return min;
    }
}


int Phi_Inverse(int t, double prob)
{
    int d;
    double Phi[DMAX+1];
    
    Phi[0] = phi[t][0];
    
    for (d = 0; d <= DMAX; d++)
    {
        if (Phi[d] >= prob)
            break;
        
        Phi[d + 1] = Phi[d] + phi[t][d + 1];
    }
    
    return d;
}

//Function H() that defines optimal protection-level
double H(int pl)
{
    double temp = h * pl;
    
    for (int dd = pl + 1; dd <= DMAX; dd++)
        temp += b * (dd - pl) * phi[0][dd];
    
    return temp;
}

//Protection Level: state transition
int X0_PL(int xz, int v0, int v1, int d0, int d1, int d2, int pl)
{
    if (xz - v0 - v1 - d0 - d1 - d2 >= pl) return xz - v0 - v1 - d0 - d1 - d2;
    else if (xz -v0 -v1 - d0 - d1 >= pl) return pl;
    else if (xz -v0 -v1 - d0 - d1 >= 0) return xz - v0 - v1 - d0 - d1;
    else if (xz -v0 -d0 >= 0) return 0;
    else return xz -v0 - d0;
}

int V0_PL(int xz, int v0, int v1, int d0, int d1, int d2, int pl)
{
    if (xz - v0 - v1 - d0 - d1 - d2 >= pl) return 0;
    else if (xz - v0 - v1 - d0 - d1 >= pl) return 0;
    else if (xz - v0 - v1 - d0 - d1 >= 0) return 0;
    else if (xz - v0 - d0 >= 0) return -(xz - v0 - v1 - d0 - d1);
    else return v1 + d1;
}

int V1_PL(int xz, int v0, int v1, int d0, int d1, int d2, int pl)
{
    if (xz - v0 - v1 - d0 - d1 - d2 >= pl) return 0;
    else if (xz - v0 - v1 - d0 - d1 >= pl) return -(xz - v0 - v1 - d0 - d1 - pl - d2);
    else if (xz - v0 - v1 - d0 - d1 >= 0) return d2;
    else if (xz - v0 - d0 >= 0) return d2;
    else return d2;
}

//one single run of simulation
void simu(void)
{
    int z[N+1];                 // order quantity in period n
    int u[N+2];                 // modified inventory position in period n
    int v_hat[N+2];             // state variable for the advance demand
    int d[N+1][T+1];            // demand vector in period n
    int x_0[N+2], v_0[N+2][T];  // inventory level and advance demand vector in period n under PL(0) heuristic
    int x_s[N+2], v_s[N+2][T];  // inventory level and advance demand vector in period n under PL(sigma) heuristic
    int x_S[N+2], v_S[N+2][T];  // inventory level and advance demand vector in period n under PL(Sigma) heuristic
    
    C_AP = 0;                   // simulated total cost under the Approximation
    C_PL_0 = 0;                 // simulated total cost under PL(0) heuristic
    C_PL_s = 0;                 // simulated total cost under PL(sigma) heuristic
    C_PL_S = 0;                 // simulated total cost under PL(Sigma) heuristic
    
    // initial state in period 1
    x_0[1] = 0;
    v_0[1][0] = 0;
    v_0[1][1] = 0;
    
    x_s[1] = 0;
    v_s[1][0] = 0;
    v_s[1][1] = 0;
    
    x_S[1] = 0;
    v_S[1][0] = 0;
    v_S[1][1] = 0;
    
    u[1] = 0; //x[1] - v[1][0] - v[1][1];
    v_hat[1] = 0; //v[1][1];
    
    for (int n = 1; n <= N; n++)
        z[n] = 0; //intialize arriving orders in period n
    
    
    // setup random demand distributions
    std::poisson_distribution<int> d0(lambda[0]);
    std::poisson_distribution<int> d1(lambda[1]);
    std::poisson_distribution<int> d2(lambda[2]);

    
    for (int n = 1; n <= N; n++)
    {
        //identify ordering quantity
        if (u[n] <= s_Opt[n][v_hat[n]])
        {
            z[n] = S_Opt[n][v_hat[n]] - u[n];
            C_AP += K;
            C_PL_0 += K;
            C_PL_s += K;
            C_PL_S += K;
        }
        
        
        //generate demand vector
        d[n][0] = d0(random_engine);
        d[n][1] = d1(random_engine);
        d[n][2] = d2(random_engine);
        
        
        //state transition under PL(0)
        x_0[n + 1] = X0_PL(x_0[n] + z[n], v_0[n][0], v_0[n][1], d[n][0], d[n][1], d[n][2], 0);
        v_0[n + 1][0] = V0_PL(x_0[n] + z[n], v_0[n][0], v_0[n][1], d[n][0], d[n][1], d[n][2], 0);
        v_0[n + 1][1] = V1_PL(x_0[n] + z[n], v_0[n][0], v_0[n][1], d[n][0], d[n][1], d[n][2], 0);
        //cost accounting under PL(0)
        C_PL_0 += l(x_0[n + 1]);
        
        //state transition under PL(sigma)
        x_s[n + 1] = X0_PL(x_s[n] + z[n], v_s[n][0], v_s[n][1], d[n][0], d[n][1], d[n][2], sigma[n]);
        v_s[n + 1][0] = V0_PL(x_s[n] + z[n], v_s[n][0], v_s[n][1], d[n][0], d[n][1], d[n][2], sigma[n]);
        v_s[n + 1][1] = V1_PL(x_s[n] + z[n], v_s[n][0], v_s[n][1], d[n][0], d[n][1], d[n][2], sigma[n]);
        //cost accounting under PL(s)
        C_PL_s += l(x_s[n + 1]);

        //state transition under PL(Sigma)
        x_S[n + 1] = X0_PL(x_S[n] + z[n], v_S[n][0], v_S[n][1], d[n][0], d[n][1], d[n][2], Sigma);
        v_S[n + 1][0] = V0_PL(x_S[n] + z[n], v_S[n][0], v_S[n][1], d[n][0], d[n][1], d[n][2], Sigma);
        v_S[n + 1][1] = V1_PL(x_S[n] + z[n], v_S[n][0], v_S[n][1], d[n][0], d[n][1], d[n][2], Sigma);
        //cost accounting under PL(S)
        C_PL_S += l(x_S[n + 1]);

        //simulate the Approximation
        C_AP += l(X_AP_II(u[n] + z[n], v_hat[n], d[n][0], d[n][1], d[n][2]));
        

        // update modified inventory position and state
        u[n + 1] = u[n] + z[n] - d[n][0] - d[n][1] - d[n][2];
        v_hat[n + 1] = d[n][2];
        
    }
    
    
}



int main()
{
    
    //Open output file
    file.open("Result_ADIF_Simulation_II.txt", fstream::app|fstream::out);
    
    if (! file)
    {
        //if fail to open the file
        cerr << "can't open output file Result_ADIF_Simulation_II.txt!" << endl;
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
    cout << "\tSigma\tsigma\tSimu_AP\tSimu_PL_0\tSimu_PL_s\tSimu_PL_S" << endl;
    file << "\tSigma\tsigma\tSimu_AP\tSimu_PL_0\tSimu_PL_s\tSimu_PL_S" << endl;
    
    
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
            {
                f[i][j][k] = -7;
                Y[i][j][k] = -7;
            }
            
            
            // cost via the DP
            double C1 = F(1, 0, 0);
            
            cout << C1 << "\t";
            file << C1 << "\t";
            
            
            //search for optimal (s, S) in period 1
            for (int j = -40; j <= 40; j++)
                for (int k = 0; k <= VMAX; k++)
                    F(1, j, k);
            
            //search for s(v),S(v) policies in all periods
            for (int i = 1; i <= N; i++)
            for (int vv = 0; vv <= VMAX; vv++)
            {
                for (int uu = UMAX; uu >= 0; uu--)
                {
                    if (Y[i][uu][vv] > uu - UMID)
                    {
                        S_Opt[i][vv] = Y[i][uu][vv];
                        s_Opt[i][vv] = uu - UMID;
                        
                        break;
                    }
                }
            }
            
            for (int vv = 0; vv <= VMAX; vv++)
            {
                cout << s_Opt[1][vv] << "\t" << S_Opt[1][vv] << "\t";
                file << s_Opt[1][vv] << "\t" << S_Opt[1][vv] << "\t";
            }

            
            
            //initialize protection level for PL(sigma)
            int ss = Phi_Inverse(0, 1 - h / b);
            if (H(ss) >= H(ss - 1))
                ss--;
            
            for (int i = 1; i < N; i++)
                sigma[i] = ss;
            sigma[N]=0;
            
            //initialize protection level for PL(Sigma)
            Sigma = Phi_Inverse(0, 0.999);

            cout << Sigma << "\t" << sigma[1] << "\t";
            file << Sigma << "\t" << sigma[1] << "\t";
            
            //====================================================================
            //start simulation
            //====================================================================
            //average over 100,000 runs
            int RUN = 100000;
            double sum_AP = 0, sum_0 = 0, sum_s = 0, sum_S = 0;
            for (int j = 1; j <= RUN; j++)
            {
                simu();
                sum_AP += C_AP;
                sum_0 += C_PL_0;
                sum_s += C_PL_s;
                sum_S += C_PL_S;
            }
            cout << sum_AP / RUN << "\t" << sum_0 / RUN << "\t" << sum_s / RUN << "\t" << sum_S / RUN << endl;
            file << sum_AP / RUN << "\t" << sum_0 / RUN << "\t" << sum_s / RUN << "\t" << sum_S / RUN << endl;
            
        }
    }
    
    file.close();
    return 0;
}