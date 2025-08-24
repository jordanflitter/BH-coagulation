# BH-coagulation
This is a Matlab code that solves the modified coagulation (Smoluchowski) equation in dense clusters that involve black holes (BHs) and neutron stars (NSs), given their initial mass function (IMF) and total number. The modified coagulation equation accounts for physical effects such as finite total mass, mass loss, ejections and delay time. The code enables tracking the cluster properties (number of black holes, their total mass, merger rate) as a function of time as well as the mass function. In addition, the code enables to combine the dynamic channel (of the cluster, where the mass function changes over time) with a static channel in which the mass function remains constant to compute the observed black hole mass function, as measured by any gravitational wave (GW) detector. Finally, the code contains a Fisher forecast section to predict the performance of GW-detectors in constraining the parameters of the model.

![BHMF](https://github.com/jordanflitter/BH-coagulation/blob/main/images/BHMF.png)
![N_vs_time](https://github.com/jordanflitter/BH-coagulation/blob/main/images/N_vs_time.png)
![M_vs_time](https://github.com/jordanflitter/BH-coagulation/blob/main/images/M_vs_time.png)
![P_vs_time](https://github.com/jordanflitter/BH-coagulation/blob/main/images/P_vs_time.png)

## Using the code
To run the code run `Main.m`. Read the documentation in that file to understand the meaning of the different parameters and adjust their values.

## Acknowledging
You are encouraged to use the code for your studies and modify it. If you use this code please cite:
* Jordan Flitter, Julian B. Munoz and Ely D. Kovetz, _"Outliers in the LIGO black hole mass function from coagulation in dense clusters"_, Mon. Not. Roy. Astron. Soc. 507 (2021) 1, 743-760 ([arXiv: 2008.10389](https://arxiv.org/pdf/2008.10389)).

## Contribution
While the code works perfectly fine, it can be significantly improved, starting from converting all the files into python.
Any contribution to the code will be greatly appreciated.
