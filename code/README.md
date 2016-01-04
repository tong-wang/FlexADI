Code
-----

In /code/Homogeneous, there are

* <ADI_Homo_I.cpp>: the homogeneous model with Advance Demand Information, Case I: T <= L + 1;
* <ADI_Homo_II.cpp>: the homogeneous model with Advance Demand Information, Case II: T > L + 1;
* <ADIF_Homo_I.cpp>: the homogeneous model with Advance Demand Information and Flexible Delivery, Case I: T <= L + 1;
* <ADIF_Homo_II.cpp>: the homogeneous model with Advance Demand Information and Flexible Delivery, Case II: T > L + 1.
                    	
The code essentially solves the Dynamic Programs developed in the paper to obtain the optimal inventory control policies and costs.
                    
In /code/Heterogeneous, there are
                    	
* <ADI_Hetero_I.cpp>: the heterogeneous model with Advance Demand Information, Case I: T <= L + 1;
* <ADI_Hetero_II.cpp>: the heterogeneous model with Advance Demand Information, Case II: T > L + 1;
* <ADIF_Simulation_I.cpp>: the heterogeneous model with Advance Demand Information and Flexible Delivery, Case I: T <= L + 1;
* <ADIF_Simulation_II.cpp>: the heterogeneous model with Advance Demand Information and Flexible Delivery, Case II: T > L + 1.
                    	
The last two involve solving the approximated Dynamic Programs to obtain inventory control policies and then running Monte Carlo simulations to compare the various Protection-Level heuristics proposed in the paper.

In addition, <ADIF_Hetero_Approx.cpp> is the stand-alone code for only solving the approximated Dynamic Program (which is already embeded in ADIF_Simulation_II.cpp), and <ADIF_Hetero_PL_S.cpp> is the code for solving the DP under the Protection-Level($\Sigma$) heuristic (the result of which is not used because we already have the simulated result of the same heuristic).
