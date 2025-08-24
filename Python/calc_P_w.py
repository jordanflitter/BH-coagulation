"""
Calculate the probability to detect merger event with mass m1,m2 and located at redshift z.

Usage:

calc_P_w(m1,m2,z,asdfile,approx=ls.IMRPhenomD)

Input--
m1,m2: binary component masses in solar mass
z: binary's redshift
asdfile: noise curve file [frequency, strain]

Optional inputs--
approx: frequency domain waveform. Default is IMRPhenomD.

Output--
Probability to detect a merger;

Author: 
Hsin-Yu Chen (hsin-yu.chen@ligo.org)

"""

import time
from numpy import *
import cosmolopy.distance as cd
import lal
import lalsimulation as ls
from scipy.interpolate import interp1d

global snr_th,cosmo
snr_th=6
cosmo = {'omega_M_0':0.308, 'omega_lambda_0':0.692, 'omega_k_0':0.0, 'h':0.678}

##generate the waveform    
def get_htildas(m1,m2,dist,
		   phase=0., 
		   df=1e-2,  
		   s1x=0.0, 
		   s1y=0.0, 
		   s1z=0.0, 
		   s2x=0.0, 
		   s2y=0.0, 
		   s2z=0.0, 
		   fmin=1.,
		   fmax=0.,
		   fref=1.,  
		   iota=0., 
		   lambda1=0., 
		   lambda2=0., 
		   waveFlags=None, 
		   nonGRparams=None, 
		   amp_order=0., 
		   phase_order=-1., 
		   approx=ls.IMRPhenomD,
		   ):
    hplus_tilda, hcross_tilda = ls.SimInspiralChooseFDWaveform(m1*lal.MSUN_SI, m2*lal.MSUN_SI, s1x, s1y, s1z, s2x, s2y, s2z, dist*(1E6 * ls.lal.PC_SI), iota, 0., 0., 0., 0., df, fmin, fmax, fref, nonGRparams, approx)
    freqs=linspace(hplus_tilda.f0,hplus_tilda.f0+(hplus_tilda.data.length-1.)*hplus_tilda.deltaF,hplus_tilda.data.length)
    return hplus_tilda.data.data,hcross_tilda.data.data,freqs

#take the h+ component to calculate SNR
def compute_horizonSNR(hplus_tilda,psd_interp,fsel,df):
	return sqrt(4.*df*sum(abs(hplus_tilda[fsel])**2./psd_interp))        
	
#calculate w (snr_th/snr)
def calc_w(m1,m2,z,asdfile,approx=ls.IMRPhenomD):

	fmin=1.
	fref=1.
	df=1e-2
	  
	input_freq,strain=loadtxt(asdfile,unpack=True,usecols=[0,1])
	interpolate_psd = interp1d(input_freq, (strain*1e-24)**2.)
	minimum_freq=maximum(min(input_freq),fmin)
	maximum_freq=minimum(max(input_freq),5000.)
	
	hplus_tilda,hcross_tilda,freqs = get_htildas((1.+z)*m1,(1.+z)*m2 ,cd.luminosity_distance(z,**cosmo) ,fmin=fmin,fref=fref,df=df,approx=approx)
	fsel=logical_and(freqs>minimum_freq,freqs<maximum_freq)
	psd_interp = interpolate_psd(freqs[fsel])
	snr=compute_horizonSNR(hplus_tilda,psd_interp,fsel,df)
	w=snr_th/snr
	return w

start = time.time()
import scipy.io as sio
Mmin = 1.0; Mmax = 300.0; dM = 1.0
zmin = 0.005; zmax = 3.1; Nz = 1
M = arange(Mmin,Mmax+dM,dM)
z = logspace(log10(zmin),log10(zmax),50); z=concatenate([[0],z])
N=size(M)
Nz=size(z)
W=zeros((N,N,Nz))
for k in arange(1,Nz):
    for i in arange(N):
        for j in arange(i,N):
            W[i][j][k]=calc_w(M[i],M[j],z[k],"aLIGO_O5_Noise_Curve.dat")
            W[j][i][k]=W[i][j][k]
            if j==i:
                print("Now calculating for k="+str(k)+",i="+str(i+1)+",j="+str(j+1)+" (out of Nz="+str(Nz-1)+",N="+str(N)+")\n")
                percents=100.*(k-1+(j+1.+i*(2.*N-i-1.)/2.)/(N*(N+1.)/2.))/Nz
                end=time.time(); Elapsed_Time=end-start; Remaining_Time=Elapsed_Time*(100/percents-1)*0.45
                Elapsed_Hours=int(floor(Elapsed_Time/3600)); Elapsed_Minutes=int(floor((Elapsed_Time-3600*Elapsed_Hours)/60)); Elapsed_Seconds=int(floor(Elapsed_Time-3600*Elapsed_Hours-60*Elapsed_Minutes))
                Remaining_Hours=int(floor(Remaining_Time/3600)); Remaining_Minutes=int(floor((Remaining_Time-3600*Remaining_Hours)/60)); Remaining_Seconds=int(floor(Remaining_Time-3600*Remaining_Hours-60*Remaining_Minutes))        
                Elapsed_Time_str=bool(Elapsed_Hours)*(str(Elapsed_Hours)+" hours, ")+bool(Elapsed_Minutes)*(str(Elapsed_Minutes)+" minutes, ")+bool(Elapsed_Seconds)*(str(Elapsed_Seconds)+" seconds")
                Remaining_Time_str=bool(Remaining_Hours)*(str(Remaining_Hours)+" hours, ")+bool(Remaining_Minutes)*(str(Remaining_Minutes)+" minutes, ")+bool(Remaining_Seconds)*(str(Remaining_Seconds)+" seconds")
                print(str(int(round(percents)))+"% of calculations have been completed.")
                print("Elapsed time: "+Elapsed_Time_str+".")
                print(bool(Remaining_Time)*("Estimation of remaining time: "+Remaining_Time_str+"."))
w_sample,P_sample=genfromtxt("Pw_single.txt",unpack=True)
P=interp1d(w_sample, P_sample,bounds_error=False,fill_value=0.0)
sio.savemat('P_values_aLIGO_O5.mat',{'M_vec':M,'z_vec':z,'P_values':P(W)})