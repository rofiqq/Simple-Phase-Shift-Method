import numpy as np
import math
import matplotlib.pyplot as plt


class PhaseShift(object):
    def __init__(self, filename, freq, vr, xo, dx, dt):
        """
        Atributes :
            filename = txt file, matrix 
                vertikal-axes = time
                horizontal-axes = distance
            freq = (fmin, fmax)
            vr = (vrmin, vrmax, vrinc)
            xo = distance source - 1st geophone (m)
            dx = distance between geophones (m)
            dt = sampling periode (s)
            fmin = minimal frequency  (Hz)
            fmax = maximal frequency (Hz)
            vrmin = min Vphase (m/s)
            vrmax = max Vphase (m/s)
            vrinc = increment for Vphase (m/s)
        """
        self.filename = filename
        self.freq = freq
        self.vr = vr
        self.xo = xo
        self.dx = dx
        self.dt = dt
        
    def Matrix(self):
        inp = open(self.filename,'r')
        matrix2 = []
        inp = inp.readlines()
        for i in range(len(inp)):
            row=[]
            for j in range(len(inp[i].split())):
                row.append(float((inp[i].split())[j]))
            matrix2.append(row)
        return np.array(matrix2)
    
    def run(self):
        
        matrix = self.Matrix()
        
        fmin, fmax = self.freq[0], self.freq[1]
        vrmin, vrmax, vrinc = self.vr[0], self.vr[1],  self.vr[2]
        xo, dx = self.xo, self.dx
        dt=self.dt;
        #Time series
        t=np.array(range(len(matrix)))*dt
        
        #Loads delta t 
        [ m, n] = matrix.shape
        
        #time series for next power of two
        np2 = (m - 1).bit_length()
        u = np.zeros((2**np2,n))
        u[0:m,0:n]=matrix
                
        #fft
        U = np.fft.fft(u, axis=0)
        
        #Use half of data (Nyquist frequency) 
        U = U[0:2**(np2-1)+1,:]
        
        #Nyquist frequency
        fnyq = 1/(2*dt)
        fvec = fnyq*(np.arange(0,(2**(np2-1))+1,1))/(2**(np2-1))
        
        
        fmini=min(np.argwhere(abs(fvec-fmin)<1))
        fmaxi=min(np.argwhere(abs(fvec-fmax)<1))
        Vtrial = np.arange(vrmin, vrmax+vrinc, vrinc)
        fi = np.arange(fmini, fmaxi+1)
        
        V = np.zeros((len(Vtrial),len(fi)),dtype=complex)
        xnf = 0
        for i in fi :
            xnf += 1
            keq = (2*math.pi*fvec[i])/Vtrial
            for j in range(len(U[1])):
                e = np.exp(1j*keq*abs(j)*dx+xo)
                V[:,xnf-1]=V[:,xnf-1]+(U[i,j]/abs(U[i,j]))*e
        
        normalized=1
        if normalized==1:
            [tmpm, tmpn]=V.shape
            for im in range(tmpn):
                V[:,im]=V[:,im]/max(abs(V[:,im]))
        ffi = fvec[fi]
        
        fgrid, vgrid = np.meshgrid(ffi, Vtrial)
        V=abs(V)
        
        # =============================================================================
        # Plot result
        # =============================================================================
        fig=plt.figure(figsize=(6,8),tight_layout=True)
        # Plot Trace
        ax = fig.add_subplot(211)
        xtick = []
        for i in range(len(matrix[0])):
            ax.plot(matrix[:,i]+((matrix.max()-matrix.min())/2)*i,t,'k')
            xtick.append(matrix[len(matrix)-1,i]+((matrix.max()-matrix.min())/2)*i)
        xticks, xticklabels = [], []
        for i in range(len(xtick)):
            if i%2==0 :
                xticklabels.append(str(i+1))
                xticks.append(xtick[i])
        plt.gca().invert_yaxis()
        plt.title('Trace')
        ax.set_xlabel('Geophone')
        ax.set_ylabel('Time (s)')
        ax.set_xticks(np.array(xticks))
        ax.set_xticklabels(xticklabels)
        # Plot Dispersion Curve
        ax2 = fig.add_subplot(212)
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel('Phase Velocity (m/s)')
        plt.contourf(fgrid, vgrid, V,50, cmap='jet')
        plt.title('Dispersion Curve')
        fig.savefig('PhaseShift.png',bbox_inches="tight",dpi=fig.dpi)
        return
    

def main():
    PS = PhaseShift('example.txt',(0,50),(10,1000,5),12,2,0.001)
    PS.run()

if __name__ == "__main__":
  main()
