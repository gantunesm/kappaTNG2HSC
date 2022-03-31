import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import os
from multiprocessing import Pool
from astropy import units as u
from astropy.convolution import Gaussian2DKernel, convolve, convolve_fft
import pymaster as nmt 
from hsc_utils import map_utils as mu
from hsc_utils import stats_code as stats 
from hsc_utils.ST_Jan22 import *
import time
 

#Input parameters/information:
sigma_e = 0.28 ## per component/ sqrt(2) 
zs = ['z1','z2','z3','z4' ,'singlez']

side = 1024  # number of grids
theta = 5.0  # opening angle in deg
map_side_deg = theta*u.degree
pixel_angular_side = map_side_deg /side
theta_smooth = [1,2,3,5,10, 15]  #arcmin 


##################For Stats computation#####################################
nbins = 20 #number of bins to compute peaks, mfs, pdf and minima
l0_bins = np.linspace(2,2900,20)# bins to compute the power spectrum- Going up to 1900 to avoid spurious signal at higher ells
binner = nmt.NmtBinFlat(l0_bins[:-1], l0_bins[1:]) #Namaster power spectrum binner 
# The effective sampling rate for these bandpowers can be obtained calling:
ells_uncoupled = binner.get_effective_ells()

#Important to compute the power spectrum:
lx_rad = theta* np.pi/180
ly_rad = theta* np.pi/180
masks_s = np.ones((side,side))
#compute workspace for power spectrum computation- this only depends on the mask!:
f0_ = nmt.NmtFieldFlat(lx_rad,ly_rad,masks_s,[masks_s]) 
wsp = nmt.NmtWorkspaceFlat()
wsp.compute_coupling_matrix(f0_,f0_,binner)
#Importat to compute the Scattering-Transform:
# define ST calculator
M = side
N = side
J = 6
L = 4

# filters_set = FiltersSet(M, N, J, L)
# filters_set.generate_morlet(if_save=True, save_dir='/global/homes/g/gmarques/hsc_ng/', precision='single')


filters_set = np.load('/home/gmarques/ng_hsc/filters_set_mycode_M'+str(M)+'N'+str(N)+'J'+str(J)+'L'+str(L)+'_single.npy',allow_pickle=True)[0]['filters_set']
ST_calculator = ST_mycode_new(filters_set, J, L )
path_stats = '/lustre/work/gmarques/kappa_TNG/stats_TNG_Dark/'


# filters_set = np.load('/global/homes/g/gmarques/hsc_ng/filters_set_mycode_M'+str(M)+'N'+str(N)+'J'+str(J)+'L'+str(L)+'_single.npy',allow_pickle=True)[0]['filters_set']
# ST_calculator = ST_mycode_new(filters_set, J, L )
# path_stats = '/global/cscratch1/sd/gmarques/tng/stats_TNGDark/'

###################################################################
def get_sigmae(zz):
    #Read neff/arcmin2
    # neff = np.load('/global/homes/g/gmarques/hsc_ng/for_baryons/neff/neff_Heymans2012_'+zz+'.npy')
    neff = np.load('/home/gmarques/ng_hsc/for_baryons/neff/neff_Heymans2012_'+zz+'.npy' )
    sigma_pix_arr =  (sigma_e / (pixel_angular_side * sqrt(neff/ u.arcmin**2))).decompose().value 
    return sigma_pix_arr

sigma_e = [get_sigmae(z) for z in zs]



def get_kernel(sig_x, sig_y, Gaussian = False):
    """
    Get kernel for a image with sig_x and sig_y (assuming assymetric). If the kernel is simply Gaussian, set 
    Gaussian = True. Otherwise, it will be the same as Oguri+18
    """
    sigma_x, sigma_y  = sig_x, sig_y
    #print(sigma_x)
    x, y = np.meshgrid(np.linspace(-15,15, int(20*sigma_x)), np.linspace(-15,15, int(20*sigma_y)))
    d = np.sqrt(x*x+y*y)

    exp_part = x**2/(sigma_x**2)+ y**2/(sigma_y**2)
    kernel = 1./(np.pi*sigma_x*sigma_y) * np.exp(-exp_part)
    kernel /= kernel.sum()

    if Gaussian ==True: 
        kernel = Gaussian2DKernel(sig_x, sig_y).array
        kernel/= np.sum(kernel)

    return kernel


def smooth_map(map_in, theta_g):
    smooth_in_pixel = (theta_g/60*u.degree * [np.shape(map_in)[0],np.shape(map_in)[1]] /[theta,theta]).value #.decompose().
    # print(smooth_in_pixel, np.shape(smooth_in_pixel))
    kernel =  mu.get_kernel(smooth_in_pixel[0],smooth_in_pixel[1], Gaussian = False)
    # print(kernel)
    smooth_image = convolve_fft(map_in,kernel)
    return smooth_image


def map_stats(out_Emode, out_Bmode, masks_s, nbins, sigmas, theta_smooth, lx_rad, ly_rad, wsp_all,binner,J, L, ST_calculator):
    ''' Warning: Function defined with global variables- this is not the best, but can be easier to parallelize with multiprocessing, for example.
    From KE and KB already computed, compute the  statistics for all smoothing scales
    '''        
    #######
    #Compute the statistics:
    
    #V0,v1,V2, pdf, peak-counts, minima  :  
    stats_threshold = [stats.get_stats(out_Emode[theta_g],masks_s ,np.linspace(-4,4,nbins)*sigmas[theta_g]) for theta_g in range(len(theta_smooth)) ]
    
    V0 = np.array(stats_threshold)[:,0,:]
    V1 = np.array(stats_threshold)[:,1,:]
    V2 = np.array(stats_threshold)[:,2,:]
    pdff = np.array(stats_threshold)[:,3,:]
    peak = np.array(stats_threshold)[:,4,:]
    minima = np.array(stats_threshold)[:,5,:]


    #Power-spectrum:
    cls_all = [stats.compute_ps_spin0_all_terms(map01=out_Emode[theta_g],map02=out_Bmode[theta_g],mask1=masks_s,lx_rad=lx_rad,ly_rad= ly_rad,b=binner,wsp=wsp_all) for theta_g in range(len(theta_smooth))]
    clee = np.array(cls_all)[:,0]
    clbb = np.array(cls_all)[:,1]
    cleb = np.array(cls_all)[:,2]
 
    #moments
    moments_all  = [stats.get_stats_moments(data =out_Emode[theta_g], mask= masks_s) for theta_g in range(len(theta_smooth))]
    sigma0= np.array(moments_all)[:,0]
    sigma1= np.array(moments_all)[:,1]
    S0= np.array(moments_all)[:,2]
    S1= np.array(moments_all)[:,3]
    S2= np.array(moments_all)[:,4]
    kur0= np.array(moments_all)[:,5]
    kur1= np.array(moments_all)[:,6]
    kur2= np.array(moments_all)[:,7]
    kur3= np.array(moments_all)[:,8]



    # calculate ST coef.
    S_out = np.empty((0,1+J+J*J*L))
    temp, _,_,_ = ST_calculator.forward(np.array(out_Emode).astype(np.float32) , J, L)
    S_out = np.concatenate((S_out, temp.cpu().numpy()), axis=0)

 
    return [V0,V1,V2,pdff, peak,minima, clee,clbb,cleb, moments_all, S_out]
         




 

 

def add_noise_stats(counter):
    print(counter ,'counter')
    for zz in range(len(zs)):

        # print(zs[zz], counter,'---zz')
        sigma_ei = sigma_e[zz]
        # print(sigma_ei)
        #Add white noise:
        # filename = '/global/cscratch1/sd/gmarques/tng/kappaTNG-Dark_HSCzs/'+zs[zz]+'/kappa%i.npy' % (counter)

        filename= '/lustre/work/gmarques/kappa_TNG/kappa_TNG_HSC/kappaTNG-Dark_HSCzs/'+zs[zz]+'/kappa%i.npy' % (counter)
        kappa_map = np.load(filename)
        noise_map = np.random.normal(loc=0.0, scale= sigma_ei, size=(side, side))

        kappa_total = kappa_map + noise_map
        #Smooth map:
        k_smoothed = [smooth_map(kappa_total, thetai) for thetai in theta_smooth]

        #Read std to compute stats:

        std_bar = np.load('/home/gmarques/ng_hsc/mean_std_fid_mocks/mean_std_wide12hallsmoothscales_'+zs[zz]+'.npy')

        # std_bar = np.load('/global/homes/g/gmarques/hsc_ng/for_baryons/std/mean_std_wide12hallsmoothscales_'+zs[zz]+'.npy')


        #Compute statistics:
        V0,V1,V2, pdff, peak,minima, clee,clbb,cleb, moments_all, S_out = map_stats(k_smoothed, np.zeros_like(k_smoothed), masks_s, nbins, std_bar, theta_smooth, lx_rad, ly_rad, wsp, binner,J, L, ST_calculator)


        np.save(path_stats+'power_spectrum/clee_'+zs[zz]+'_nsim_'+str(counter)+'.npy', clee)
        np.save(path_stats+'pdf/pdf_'+zs[zz]+'_nsim_'+str(counter)+'.npy',pdff)
        np.save(path_stats+'peaks/peaks_'+zs[zz]+'_nsim_'+str(counter)+'.npy',  peak)
        np.save(path_stats+'minima/minima_'+zs[zz]+'_nsim_'+str(counter)+'.npy', minima)
        np.save(path_stats+'moments/moments_'+zs[zz]+'_nsim_'+str(counter)+'.npy',  moments_all)
        np.save(path_stats+'MFs/V0_'+zs[zz]+'_nsim_'+str(counter)+'.npy', V0)
        np.save(path_stats+'MFs/V1_'+zs[zz]+'_nsim_'+str(counter)+'.npy', V1)
        np.save(path_stats+'MFs/V2_'+zs[zz]+'_nsim_'+str(counter)+'.npy', V2)
        np.save(path_stats+'ST/ST_'+zs[zz]+'_nsim_'+str(counter)+'.npy', S_out)

    return
 
        # istds_noisy = [k_smoothed[i].std() for i in range(5)]
        # #### plot the no ise maps
        # plt.figure(figsize=(16,2))
        # for i in range(5):
        #     plt.subplot(1,5,i+1)
        #     istd=istds_noisy[i]
        #     plt.imshow(k_smoothed[i].data, vmin=-3*istd, vmax=3*istd, origin='lower', extent=[0,10,0,10], cmap="Blues")
        #     # plt.title('Noisy sim z=%.1f'%(zz))
        #     plt.colorbar()
        #     plt.tight_layout()
        #     plt.savefig('/global/homes/g/gmarques/hsc_ng/for_baryons/kappa_smoothed_z2.pdf')

        #Compute statistics:

 
 
 

# pool = Pool() 
# out = pool.map(add_noise_stats, np.arange(1, 500))
# pool.close()


def main():
    start= time.time()
    args = np.arange(1, 1999)
    nthreads = 52
    with Pool(nthreads) as pool:
        pool.map(add_noise_stats,  args)
        pool.close()
        pool.join()

    runtime = time.time() -start 
    print(f'Time: {runtime}s')

if __name__== '__main__':
    main() 


# def main():
#     start= time.time()
#     args = np.arange(1,1000)
#     nthreads = 52
#     with Pool(nthreads) as pool:
#         pool.map(add_noise_stats,  args)
#         pool.close()
#         pool.join()

#     runtime = time.time() -start 
#     print(f'Time: {runtime}s')

# if __name__== '__main__':
#     main()