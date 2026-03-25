import multiprocessing
import sys
import h5py
import numpy as np
#from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors

#### Input Parameters ####
ColourMapL = 4
ColourMapR = -2
ColourMapC = 1
QL = 4
QR = 1
QC = 0
#0 = Density; 1 = Pressure; 2 = UR; 3 = UP; 4 = UU; 5 = Gamma; 6 = GammaBeta; 7 = Luminosity; 8 = Gammah

def printname(name):
    print(name)

print(sys.argv)
print(int(sys.argv[1]))
num_pross = int(sys.argv[1])
checkpoint_initial = 0
checkpoint_final = 300

def plot_2D(chkpts_min,chkpts_max):
    check_points = np.arange(chkpts_min,chkpts_max,1)
    total_chkpts = len(check_points)
    
    klNe_vmax = 0.3
    
    for fn in range(total_chkpts):
        check_point_num = check_points[fn]
        if(check_point_num<10):
            f = h5py.File("Checkpoints/checkpoint_000%d.h5" %check_point_num,'r')
        elif(check_point_num<100):
            f = h5py.File("Checkpoints/checkpoint_00%d.h5" %check_point_num,'r')
        elif(check_point_num<1000):
            f = h5py.File("Checkpoints/checkpoint_0%d.h5" %check_point_num,'r')
        else:
            f = h5py.File("Checkpoints/checkpoint_%d.h5" %check_point_num,'r')
            
        print("reading checkpoint: ", f)
        
        #f.visit(printname)
        
        #Reading data
        Data_cells = f["Data/Cells"]
        Nr = f["Grid/Nr"]   #Size of each radial track
        T = f["Grid/T"]     #Time of the data
        p_kph = f["Grid/p_kph"] #1D phi data
        t_jph = f["Grid/t_jph"] #1D theta data
        Index = f["Grid/Index"] #Index of each region
        
        radial_track_length = Nr[:,]
        theta = t_jph[:,]
        phi = p_kph[:,]
        index = Index[:,]
        data = Data_cells[:,:]
        time = T[:,]
        
        rmin = data[0,-1]
        rmax = data[index[1,0]-1,-1]
        radial_res_max = max(radial_track_length)
        radial_res_max = radial_res_max[0]
        RADIUS = np.geomspace(rmin, rmax, radial_res_max)
        THETA = np.zeros(len(theta)-1)
        
        RHO = np.zeros([len(RADIUS),len(THETA)])
        PPP = np.zeros([len(RADIUS),len(THETA)])
        UR = np.zeros([len(RADIUS),len(THETA)])
        UP = np.zeros([len(RADIUS),len(THETA)])
        UU = np.zeros([len(RADIUS),len(THETA)])
        GAM = np.zeros([len(RADIUS),len(THETA)])
        GAM_Beta = np.zeros([len(RADIUS),len(THETA)])
        LUM = np.zeros([len(RADIUS),len(THETA)])
        GAM_h = np.zeros([len(RADIUS),len(THETA)])
        
        QUANT = [RHO,PPP,UR,UP,UU,GAM,GAM_Beta,LUM,GAM_h]
        colourmap = ['jet','tab20c','plasma_r','inferno','magma','seismic','PiYG','hot','jet_r']
        
        #Arrays to store data
        for n in range(len(radial_track_length)):
            THETA[n] = 0.5*(theta[n]+theta[n+1])
            x = np.zeros(radial_track_length[n,0])  #x-axis for that radial region
            param = np.zeros([radial_track_length[n,0],6])  #fluid parameters
            
            for i in range(radial_track_length[n,0]):
                j = i+index[n,0]
                if i==0:
                    x[i] = data[j,-1]/2.
                    param[i,:] = data[j,:]
                else:
                    x[i] = data[j-1,-1]+(data[j,-1]-data[j-1,-1])/2.       #r plus half/ upper edge of radial region
                    param[i,:] = data[j,:]        #data at r
                    
            rho = param[:,0]
            ppp = param[:,1]
            ur = param[:,2]
            up = param[:,3]
            rho = param[:,0]
            gam = np.zeros(len(x))
            v = np.zeros(len(x))
            lum = np.zeros(len(x))
            
            uu = np.sqrt(ur*ur+up*up)   #total bulk velocity
            
            for i in range(len(x)):
                #if(uu[i] >= 1.0):
                    gam[i] = np.sqrt(1.0+uu[i]*uu[i])
                    v[i] = np.sqrt( 1. - (1./(gam[i]*gam[i])))
            gam_beta = np.where(uu > np.ones(len(uu)),uu,0.0)
            
            h = 1.+4.*np.divide(ppp,rho)
            gamh = gam*h
            
            rel_indx = np.where(uu>1.0)
            lum[rel_indx] = 4*np.pi*np.square(x[rel_indx])*rho[rel_indx]*h[rel_indx]*gam[rel_indx]*gam[rel_indx]*v[rel_indx]
            
            #Interpolate
            RHO[:,n] = np.interp(RADIUS,x,rho)
            PPP[:,n] = np.interp(RADIUS,x,ppp)
            UR[:,n] = np.interp(RADIUS,x,ur)
            UP[:,n] = np.interp(RADIUS,x,up)
            UU[:,n] = np.interp(RADIUS,x,uu)
            GAM[:,n] = np.interp(RADIUS,x,gam)
            GAM_Beta[:,n] = np.interp(RADIUS,x,gam_beta)
            LUM[:,n] = np.interp(RADIUS,x,lum)
            GAM_h[:,n] = np.interp(RADIUS,x,gamh)
    
        #Plot
        #**** Plot Beautification with rcParams ****

        plt.rcParams["font.weight"] = "bold"
        #plt.rcParams["axes.labelweight"] = "bold"
        plt.rcParams["text.color"] = "white"
        plt.rcParams["axes.labelcolor"] = "white"
        plt.rcParams["axes.labelsize"] = 16
        plt.rcParams["axes.linewidth"] = 0
        plt.rcParams["xtick.color"] = "white"
        plt.rcParams["xtick.major.size"] = 8
        plt.rcParams["xtick.major.width"] = 3
        plt.rcParams["xtick.minor.size"] = 4
        plt.rcParams["xtick.minor.width"] = 2
        plt.rcParams["xtick.labelsize"] = 12
        plt.rcParams["xtick.direction"] = "out"
        plt.rcParams["ytick.color"] = "white"
        #plt.rcParams["ytick.labelcolor"] = "white"
        plt.rcParams["ytick.major.size"] = 8
        plt.rcParams["ytick.major.width"] = 3
        plt.rcParams["ytick.labelsize"] = 12
        plt.rcParams["ytick.direction"] = "out"
        plt.rcParams["legend.facecolor"] = "black"
        plt.rcParams["legend.fontsize"] = 16

        fig = plt.figure(figsize=[24,20],facecolor='k')
        fig.subplots_adjust(bottom=-1.0)
        plt.style.use('dark_background')
        ax = fig.add_axes([0.001,0.20,0.98,0.95],polar=True)
        ax_1d = fig.add_axes([0.35,0.065,0.4,0.27])
        
        
        #What to Plot
        yL = QUANT[QL]
        yR = QUANT[QR]
        yC = QUANT[QC]
        
        yL_1D = yL[:,0]
        yR_1D = yR[:,0]
        yC_1D = yC[:,0]
        
        RMAX_1 = 1.2*RADIUS[np.where(QUANT[4][:,0] < np.ones(len(QUANT[4][:,0])))[0][-1]]
        RMAX_2 = klNe_vmax*time[0]
        RMAX = max(RMAX_1,RMAX_2)
        
        #Making meshgrid
        TT, RR = np.meshgrid(THETA,RADIUS)
        
        R_photo = (klNe_vmax*time[0])*np.ones(len(THETA))
        ax.plot(THETA,R_photo,'white',linestyle='--',linewidth=4)
        ax.plot((2*np.pi-THETA),R_photo,'white',linestyle='--',linewidth=4)
        
        ax_1d.vlines(klNe_vmax*time[0],min(min(yL_1D),min(yR_1D)),1e2*max(max(yL_1D),max(yR_1D)),'white','--',linewidth=4)
        
        #ax.grid(True)
        ax.grid(False)
        
        vmin = yL.min()
        vmax = yL.max()
        if(vmin<=0):
            vmin = 1e-10
        if(vmax<=0):
            vmax = vmin
        #vmin = 1e-10*vmax
            
        LP = ax.pcolormesh(THETA,RADIUS,yL,norm=colors.LogNorm(vmin, vmax),cmap=colourmap[ColourMapL])
        fig.colorbar(LP,ax=ax,location='left',shrink=0.5,pad=6e-2)
        
        vmin = yR.min()
        vmax = yR.max()
        if(vmin<=0):
            vmin = 1e-10
        if(vmax<=0):
            vmax = vmin
        if(QR==4):
            vmin = 1.0

        RP = ax.pcolormesh((2*np.pi-THETA),RADIUS,yR,norm=colors.LogNorm(1.0, vmax),cmap=colourmap[ColourMapR])
        fig.colorbar(RP,ax=ax,location='right',shrink=0.5,pad=8e-2)
        
        #--- Plotting 1D along axis ---
        ax_1d.plot(RADIUS,yR_1D,linewidth=4,color='r',label="Right")
        ax_1d.plot(RADIUS,yL_1D,linewidth=4,color='b',label="Left")
        ax_1d.plot(RADIUS,np.divide(yC_1D,max(yC_1D)),linewidth=4,color='g',label=r"$(\rho/\rho_{\rm{max}})$")
        ax_1d.legend(frameon=False)
        
        ax_1d.set_xscale('log')
        ax_1d.set_yscale('log')
        
        #yC_filtered = np.where(yL > 1e-2*np.max(yL), yL, 1e-10)
        #levels = np.logspace(np.log10(1e-2*yC_filtered.max()), np.log10(yC_filtered.max()), 2)
        #print(levels)
        if(np.max(yL)>0.0):
            levels = np.array([1e-2*np.max(yL),1e-1*np.max(yL),np.max(yL)])
            ax.contour(THETA,RADIUS,yL,levels=levels,linewidths=1.5,cmap=colourmap[ColourMapC])
        if(np.max(yR)>1.0):
            c1 = min(1.0,1e-1*np.max(yR))
            c2 = max(1.0,1e-1*np.max(yR))
            levels = np.array([c1,c2,np.max(yR)])
            #print(vmax,np.max(yR),yC_filtered.max(),1e-2*np.log10(yC_filtered.max()),np.log10(yC_filtered.max()),levels)
            ax.contour((2*np.pi-THETA),RADIUS,yR,levels=levels,linewidths=1.5,cmap=colourmap[ColourMapC])
        #cs = ax.contour(THETA,RADIUS,yC_filtered, locator=ticker.LogLocator(), cmap='viridis')
        #cbar = fig.colorbar(cs)
        #ax.contour(THETA,RADIUS,yC)
        
        #ax.set_title("Time = %.2f"%time[0],fontweight='bold',fontsize=18,color='white')
        fig.text(0.45,0.95,"Time = %.2e"%time[0],fontweight='bold',fontsize=18,color='white')
        fig.text(0.80,0.30,r"$T_{obs}$ = %.2e"%(time[0]*(1-klNe_vmax)),fontweight='bold',fontsize=14,color='white')
        #ax.set_theta_direction(-1)
        ax.set_theta_offset(np.pi/2.0)
        #ax.set_thetalim((3./2)*np.pi,np.pi/2.0)
        ax.set_thetamin(-90)
        ax.set_thetamax(90)
        ax.set_rmax(1.1*RMAX)
        ax.set_rmax(1.1*RMAX)
        ax_1d.set_xlim(min(RADIUS),1.1*RMAX)
        MAX_val = [max(yL_1D),max(yR_1D),max(yC_1D)]
        if(min(MAX_val)>0):
            ax_1d.set_ylim(1e-4*min(MAX_val),10.0*max(MAX_val))
        else:
            ax_1d.set_ylim(1e-4*max(MAX_val),10.0*max(MAX_val))
        ax_1d.set_xlabel("Distance (c.u.)")
        ax_1d.set_ylabel("Value (c.u.)",color='white')
        ax.grid(True)
        fig.legend(["Luminosity"],loc=[0.05,0.92],frameon=False)
        fig.legend([r"$\Gamma\beta_{rel}$"],loc=[0.82,0.92],frameon=False)
        ax.set_facecolor("black")
        
        fig.set_size_inches(13,9)
        plt.savefig('VIS_2D/density_4v_2D_%04d.png' %check_point_num, dpi=240)
        #plt.show()
        plt.close()
        
        f.close()

if __name__ == '__main__':
    pool = multiprocessing.Pool(num_pross)
    chkpts_range = np.linspace(checkpoint_initial,checkpoint_final,num_pross,dtype=int)
    result_async = [pool.apply_async(plot_2D,args=(chkpts_range[i],chkpts_range[i+1])) for i in range(len(chkpts_range)-1)]
    pool.close()
    results = [r.get() for r in result_async]
    pool.join()
