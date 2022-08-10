# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 20:31:45 2021
@author: erich

FDTD transmison line simulator:
    This software is capable of simulating the time_dependent transmission line 
    currents and voltages on a uniform lossy transmision line by using second 
    order accurate FDTD approximation to the transmision line equations.
     
Written by: Erich Wanzek
Written 2/17/2021
Intro to Computational Electromangtics
University of Colrado Denver
"""
##############################################################################
##############################################################################
'''Import Required Libraries'''
import numpy as np
import math as m
import matplotlib.pyplot as plt

plt.rcParams['figure.dpi']=500
plt.style.use('default')

#############################################################################
import matplotlib.animation as animation
##############################################################################
def Trapezoidal_pulse_inline(t):
    Vs=0
    t=dt*(t-0.5)
    
    if(t <= 0):
        Vs=0;
    elif( t< risetime):
        Vs = t * 1/(risetime)
    elif( t < risetime + duration):
        Vs = 1
    elif( t < risetime + duration + falltime):
        Vs =(risetime +duration +falltime -t)/falltime
    return Vs_amp*Vs         
##############################################################################
def Gaussian_pulse_inline(t):
    Vs=0
    t=dt*(t+0.5)
    Vs=np.sin(t*10**10) #*np.exp(-(((t)-pulse_delay)**2)/((pulse_hw)**2))
    return Vs_amp*Vs
##############################################################################
def Trapezoidal_pulse(timestart,risetime,duration,falltime,nt,dt,Vs_amplitude):
    trap_pulse=np.zeros((1,nt))
    ts= m.floor(timestart/dt);
    rt= m.floor(risetime/dt);
    d= m.floor(duration/dt);
    ft= m.floor(falltime/dt);
    
    if ((ts+rt+d+ft) > nt):
        print('Notice: total pulse time longer than simulaiton time')
    for t in range(0,ts):
        trap_pulse[0][t]=0;
        
    for t in range(ts,(ts+rt)):
        trap_pulse[0][t]=((1)/((ts+rt)-(ts+1)))*(t-ts);
        
    for t in range((ts+rt),(ts+rt+d)):
        trap_pulse[0][t]=1;
        
    for t in range((ts+rt+d),(ts+rt+d+ft)):
        trap_pulse[0][t]=1-(((1)/((ts+rt+d+ft)-(ts+rt+d+1)))*(t-(ts+rt+d)));
        
    for t in range((ts+rt+d+ft),nt-1):
        trap_pulse[0][t]=0;
         
    return Vs_amplitude*trap_pulse
##############################################################################
def Gaussian_pulse():
    gauss_pulse=np.zeros((1,nt))
    to=pulse_delay
    tw=pulse_hw
 

    #if ((to) > nt):
    #    print('Notice: total pulse time longer than simulaiton time')
    for t in range(0,nt):
        
        gauss_pulse[0][t]=np.exp(-(((t*dt)-to)**2)/((tw)**2))
    
    return Vs_amp*gauss_pulse           
##############################################################################
def TimeAdvance(signature):
    '''This funciton time-advance the line voltages and currents with in 
    line update funcitons.'''
    for t in range(0,nt):
        V_line_update();    #update line voltages:interior node voltages'
        #V_source_update(signature[0][t]);  #update source node voltage'
        V_source_update2(t)
        V_load_update();     #update load node voltage '
        
        I_line_update();      #update the interior line currents'
        
        save_output((V[0][0]),(V[0][nx-1]),t); #'Save and output probe voltages and currents'
        
    return None
##############################################################################
def V_line_update():
    """THis fucniton isa finciton the computes the time update of the line
        voltage in the interior non of the homogenous transmission line"""
    for k in range(1,nx-1):  
        V[0][k] = cv1*V[0][k] - cv2*(I[0][k]-I[0][k-1]) 
        
    return None
##############################################################################
def I_line_update():
    """THis fucniton isa finciton the computes the time update of the line
        current in the interior of the homogenous transmission line"""
    for k in range(0,nx-1):
        I[0][k] = ci1*I[0][k] - ci2*(V[0][k+1]-V[0][k])
   
    return None         
##############################################################################
def V_source_update(Vs_sig):
    """This function is a function that computes the time update of the source
    voltage at the source end of the transmission line, V(1)"""
    if(Rg > 0):
        V[0][0] = c1g * (c2g * V[0][0] - I[0][0] + Vs_sig/Rg)
    else:
        V[0][0] = Vs_sig  
    return None         
##############################################################################
def V_source_update2(t):
    """This function is a function that computes the time update of the source
    voltage at the source end of the transmission line, V(1)"""

    if (source_sig == 'trapezoidal'):
        Vs_sig=Trapezoidal_pulse_inline(t)
        if(Rg > 0):
            V[0][0] = c1g * (c2g * V[0][0] - I[0][0] + Vs_sig/Rg)
        else:
            V[0][0] = Vs_sig 
        
    if (source_sig == 'gaussian'):
        Vs_sig=Gaussian_pulse_inline(t)
        if(Rg > 0):
            V[0][0] = c1g * (c2g * V[0][0] - I[0][0] + Vs_sig/Rg)
        else:
            V[0][0] = Vs_sig 
     
    return None   
##############################################################################
def V_load_update():
    """This function is a function that computes the time update of the source
    voltage at the source end of the transmission line, V(1)"""
    if(RL== float('inf')):
        V[0][nx-1] = V[0][nx-2];
        I[0][nx-2]=0
    if(RL==0):
        V[0][nx-1]=0.0;
       
       
    elif(load_type==0.0):   # R Load
        V[0][nx-1] = c1l * (c2l * V[0][nx-1] + I[0][nx-2]);  
    
    elif(load_type==1.0):   #Parallel R-C Load
        V[0][nx-1] = (cPRC1 * V[0][nx-1]) +  (cPRC2 *I[0][nx-2]);
        
    elif(load_type==2.0):   #Series R-L load
        Y[0,0]=I[0][nx-2]
        Y[1,0]=0.0
      
        X[:,:]=np.matmul(Pinv,(np.matmul(Q_matrix,X[:,:]) + Y))
      
        V[0][nx-1]=X[0,0]
       
    return None         
##############################################################################
def write_array_source_voltage(Vs,t):
    V_source[0][t]=Vs
    return None
##############################################################################
def write_array_load_voltage(Vl,t):
    V_load[0][t]=Vl
    return None
##############################################################################
def write_txline_V_I(t):
    V_line_matrix[t,:]=V[0]
    I_line_matrix[t,:]=I[0]
    return None
##############################################################################
def save_output(Vs,Vl,t):
    if  plot_load:
        write_array_load_voltage(Vl,t)
    if  plot_source:
        write_array_source_voltage(Vs,t)
    if  make_vid or plot_line_V:
        write_txline_V_I(t)
    return None    
##############################################################################
def plot_Txline(): 
    print(len(V[0]))
    fig, ax =plt.subplots(2)
    ax[0].plot(X1[0],V[0],label='Line Voltage',color='b')
    ax[0].grid(linestyle='--')
    ax[0].legend()
    ax[0].set_ylim(-Vs_amp,+Vs_amp)
    ax[0].set_xlim(-0.01,d)
    ax[0].set_xlabel('Distance m')
    ax[0].set_ylabel('Voltage Volts')
    ax[0].set_title('Line Voltage' + ' Zo=' + str(Zo))
    
    ax[1].plot(X2[0],I[0],label='Line Current',color='r')
    ax[1].grid(linestyle='--')
    ax[1].legend()
    ax[1].set_ylim(-0.01*Vs_amp,+0.01*Vs_amp)
    ax[1].set_xlim(-0.01,d)
    ax[1].set_xlabel('Distance m')
    ax[1].set_ylabel('Current Amps')
    ax[1].set_title('Line Current')
    
    fig.set_figheight(10)
    fig.set_figwidth(12)
    #fig.suptitle('Uniform TransmissionLine V and I')
    plt.show()

    return fig

##############################################################################
def plot_line_voltage(): 
    t=int(np.floor(snap_shot_time/dt))
    time_here=((t*dt)/10**-9)
    if t ==nt:
        t=t-1
            
    if(load_type==0.0):   # R Load
        title='R Load:' + ' RL= ' +str(RL) + ' \u03A9 '  
    
    elif(load_type==1.0):   #Parallel R-C Load
        title=('Parallel RC Load:' + ' CL= ' +str(Cl/(10**-12)) + ' pF'  + ' RL= ' +str(RL) + ' \u03A9')
        
    elif(load_type==2.0):   #Series R-L load
        title=('Series RL Load:' + ' LL= ' +str(Ll/(10**-9)) + ' nH'  + ' RL= ' +str(RL) + ' \u03A9')
    
    fig, ax =plt.subplots()
    ax.plot(X1[0],V_line_matrix[t,:],label='Line Voltage' +' %.2f' %time_here + ' ns',color='k',linewidth=2)
    ax.grid(linestyle='--')
    ax.legend()
    ax.set_ylim(-Vs_amp,+Vs_amp)
    ax.set_xlim(0.0,d)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('Line Voltage (V)')
    ax.set_title(title + ' Rg=' + str(Rg)+ ' \u03A9 ')
    fig.suptitle('Uniform Transmission Line' + ' Zo=' + str(Zo) + ' \u03A9')
    plt.show()
    
    return None 
##############################################################################
def plot_source_signature(tpulse): 
    
   
    fig, ax =plt.subplots()
    ax.plot((time[0])*10**9,tpulse[0],label='Line V',color='b')
    ax.grid(linestyle='--')
    ax.legend()
    ax.set_ylim(-Vs_amp*2,+Vs_amp*2)
    ax.set_xlim(0,simtime)
    ax.set_xlabel('Time')
    ax.set_ylabel('Voltage Volts')
    ax.set_title('Source Voltage Pulse Signature')
    plt.show()
    
    return None
##############################################################################
def plot_source_voltage(): 
    if(load_type==0.0):   # R Load
        title='R Load:' + ' RL= ' +str(RL) + ' \u03A9 '  
    
    elif(load_type==1.0):   #Parallel R-C Load
        title=('Parallel RC Load:' + ' CL= ' +str(Cl/(10**-12)) + ' pF'  + ' RL= ' +str(RL) + ' \u03A9')
        
    elif(load_type==2.0):   #Series R-L load
        title=('Series RL Load:' + ' LL= ' +str(Ll/(10**-9)) + ' nH'  + ' RL= ' +str(RL) + ' \u03A9')
    
    y_max=1.1*np.amax(V_source[0])
    fig, ax =plt.subplots()
    ax.plot((time[0])*10**9, V_source[0],label='Line V at Source',color='k')
    ax.grid(linestyle='--')
    ax.legend()
    ax.set_ylim(-y_max,+y_max)
    ax.set_xlim(0,simtime*(10**9))
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Line Voltage (V)')
    ax.set_title(title + ' Rg=' + str(Rg)+ ' \u03A9 '  + ' Zo=' + str(Zo) + ' \u03A9')
    fig.suptitle('Line Voltage at Source vs Time')
    
    return None
##############################################################################
def plot_load_voltage(): 
    if(load_type==0.0):   # R Load
        title='R Load:' + ' RL= ' +str(RL) + ' \u03A9 '  
    
    elif(load_type==1.0):   #Parallel R-C Load
        title=('Parallel RC Load:' + ' CL= ' +str(Cl/(10**-12)) + ' pF'  + ' RL= ' +str(RL) + ' \u03A9')
        
    elif(load_type==2.0):   #Series R-L load
        title=('Series RL Load:' + ' LL= ' +str(Ll/(10**-9)) + ' nH'  + ' RL= ' +str(RL) + ' \u03A9')
    y_max=1.1*np.amax(V_load[0])   
    fig, ax =plt.subplots()
    ax.plot((time[0])*(10**9),V_load[0],label='Line V at Load',color='k')
    ax.grid(linestyle='--')
    ax.legend()
    ax.set_ylim(-y_max, y_max)
    ax.set_xlim(0,simtime*(10**9))
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Line Voltage (V)')
    ax.set_title(title + ' Rg=' + str(Rg)+ ' \u03A9 '  + ' Zo=' + str(Zo) + ' \u03A9')
    fig.suptitle('Line Voltage at Load vs Time')
    
    return None
##############################################################################
def V_source_signature_setup(signature,nt,dt,Vs):
    
    if signature == 'trap':
        Vs_pulse_signature=Trapezoidal_pulse()
    elif signature == 'gauss':
        Vs_pulse_signature=Gaussian_pulse()
    elif signature != ('gauss' or 'trap') :
        print('NOTICE: Must select trap or gauss for signature! Now defaulting to gauss signature')
        Vs_pulse_signature=Gaussian_pulse(1*10**-11,1*10**-15,nt,dt,Vs)
    
    return Vs_pulse_signature  
##############################################################################
def readin_parameter_data_file():
    """The function readin_parameter_data_file() is a function that reads in 
    the inital simualtion variables from the text file called
    Tx_Line_simulation_parameters.txt. The software user must edit the value 
    of the variables in this txt file and save the file before running the 
    python simulation script. The function readin_parameter_data_file() is 
    called once at the beginning of the main() function"""
    
    variables={}
    with open('Tx_Line_simulation_parameters.txt') as f:
         for line in f:
             name, value = line.strip().split('=')
           
             if name == 'source_signature':
                 variables[name]= value
             elif name == 'make_video':
                 variables[name]= value
             elif name == 'graph_source_voltage':
                 variables[name]= value 
             elif name == 'graph_load_voltage':
                 variables[name]= value
             elif name == 'graph_line_voltage':
                 variables[name]= value     
             else:    
                 variables[name]= float(value)
             
    global L; L = variables['L']
    global C; C = variables['C']
    global R; R = variables['R']
    global G; G = variables['G']
    global d; d = variables['d']
    global ncells; ncells = int(variables['ncells'])
    global cfln; cfln = variables['cfln']
    global simtime; simtime = variables['simtime']
    global Rg; Rg = variables['Rg']
    global RL; RL = variables['RL']
    global Vs_amp; Vs_amp = variables['Vs_amplitude']
    global pulse_delay; pulse_delay = variables['pulse_delay']
    global pulse_hw; pulse_hw = variables['pulse_half_width']
    global risetime; risetime = variables['risetime']
    global duration; duration = variables['duration']
    global falltime; falltime = variables['falltime']
    global load_type; load_type = variables['load_type']
    global Cl; Cl = variables['C_load']
    global Ll; Ll = variables['L_load']
    global source_sig; source_sig= variables['source_signature']
    global make_vid; make_vid= variables['make_video']
    global plot_source; plot_source= variables['graph_source_voltage']
    global plot_load; plot_load= variables['graph_load_voltage']
    global plot_line_V; plot_line_V= variables['graph_line_voltage']
    global snap_shot_time; snap_shot_time=variables['snap_shot_time']

    return None
##############################################################################
def initialize_parameters():
    """initialize_parameters() precalculates the transimision line parameters
    from the inputted data.This funciton is performed once when 
    initialize_parameters() is called in setup_FDTD_sim()"""
    
    '''precompute transmission line parameters'''
    '''compute transmission line paramters'''
    global Zo; Zo=np.sqrt(L/C);
    global vo; vo=1.0/(np.sqrt(L*C));
    
    'compute grid discretization parameters'
    global dx; dx= d/ncells;
    global dt; dt= cfln * (dx / vo);
    global nx; nx= ncells + 1;
    global nt; nt= m.floor(simtime/dt);
    
    return None
##############################################################################
def initialize_arrays():
    """initialize_arrays() initializes and sets up the data arrays for space,
    time, current and voltage variables. This funciton is performed once when 
    initialize_arrays() is called in setup_FDTD_sim()"""
    
    '''Initialize distance and time array'''
    global X1; X1=np.zeros((1,ncells+1))
    global X2; X2=np.zeros((1,ncells))
    global time; time=np.zeros((1,nt))
   
    '''Populate distance and time arrays with space/time values'''
    for i in range(0,ncells+1):
        X1[0][i]=(dx*i)
        
    for i in range(0,ncells):    
        X2[0][i]=(dx*(i-0.5))

    for t in range(0,nt):    
        time[0][t]=(dt*(t))
    
    'Initialize the line voltages and current arrays'
    global V; V=np.zeros((1,nx));   global I; I=np.zeros((1,nx-1)); 
    
    "Initialize Source Voltage and Load voltage vs time storage arrays"
    global V_source; V_source=np.zeros((1,nt)); 
    global V_load; V_load=np.zeros((1,nt));
    
    "Initialize Source Voltage and interior transmision line storage array Matrix"
    global V_line_matrix; V_line_matrix=np.zeros((nt,nx)); 
    global I_line_matrix; I_line_matrix=np.zeros((nt,nx-1));
    
    return None
##############################################################################
def calculate_prefactors():
    """calculate_prefactors() function precalcualtes the FDTD update equation
    prefactors. This funciton is performed once when calcualte_prefactors()
    is called in setup_FDTD_sim()"""
    
    'Precalcualte multiplier coefficients'
    'Precalcualte V/I update equation coefficients'
    global cv1; cv1=((C/dt)-(G/2.0))/((C/dt)+(G/2.0));
    global cv2; cv2=((1.0/((C/dt)+(G/2.0))))*(1/dx);

    global ci1; ci1=((L/dt)-(R/2.0))/((L/dt)+(R/2.0));
    global ci2; ci2=((1.0/((L/dt)+(R/2.0))))*(1/dx);

    'Precalculate Source generator update equation coefficients'
    if (Rg > 0.0):
        global b1g; b1g= C * dx * 0.5/dt;  
        global b2g; b2g= 0.5 / Rg;
        global c1g; c1g= 1.0 / (b1g + b2g);
        global c2g; c2g= b1g - b2g;
    

    'Precalculate Source generator update equation coefficients'
    if (RL > 0.0):
        global b1l; b1l= C * dx * 0.5/dt;  
        global b2l; b2l= 0.5 / RL;
        global c1l; c1l= 1.0 / (b1l + b2l);
        global c2l; c2l= b1l - b2l;  
    
    return None
##############################################################################
def load_setup():
    
    'Precalcualte multiplier coefficients for Parallel RC load update'
    if (RL > 0.0):
        bPRC1=(((C*dx)+(2*Cl))/(2*dt))
        bPRC2=((0.5*G*dx)+(1/RL))
        global cPRC1; cPRC1=((bPRC1-(0.5*bPRC2))/(bPRC1+(0.5*bPRC2)))
        global cPRC2; cPRC2=(1/(bPRC1+(0.5*bPRC2)))
    
    
    'Precalcualte multiplier coefficients for Parallel RC load update'
    P11=(((C*dx)/(2*dt))+(0.25*(G*dx)))
    P12=0.5
    P21=0.5
    P22=-1*((Ll/dt)+(0.5*RL))
    
    Q11=(((C*dx)/(2*dt))-(0.25*(G*dx)))
    Q12=-0.5
    Q21=-0.5
    Q22=-1*((Ll/dt)-(0.5*RL))
    
    global P_matrix; P_matrix=np.zeros((2,2))
    global Q_matrix; Q_matrix=np.zeros((2,2))
  
    P_matrix[0,0]=P11
    P_matrix[0,1]=P12
    P_matrix[1,0]=P21
    P_matrix[1,1]=P22
    
    Q_matrix[0,0]=Q11
    Q_matrix[0,1]=Q12
    Q_matrix[1,0]=Q21
    Q_matrix[1,1]=Q22
    
    
    global Pinv
    Pinv=np.linalg.inv(P_matrix)
    global X; X=np.zeros((2,1))
    global Y; Y=np.zeros((2,1))

    
    return None
##############################################################################
def setup_FDTD_sim():
    """setup_FDTD_sim function is the main setup fucntion that runs once in 
    the main funciton, when main() is called. This function runs the setup 
    subroutines to initialize the data arrays, calcualte simulation 
    parameters, and calcualte FDTD update equation prefactors"""
    
    initialize_parameters()
    initialize_arrays()
    calculate_prefactors()
    load_setup()              # Rl, Rl+Cl, RL+LL

    return None
##############################################################################
def post_processing():
    """The function post_processing() is a function that calls on subfunctions
    to produce various plots and/or videos of the date outputed by the FDTD
    simulation performed by Time_advance()fucniton. The function 
    post_processing() is called at the end of the main() function when all 
    FDTD simulation functions have been performed"""
    
    #plot_source_signature(tpulse)

    if plot_source == 'yes':
        plot_source_voltage()
    if plot_load == 'yes':
        plot_load_voltage()
    if plot_line_V == 'yes':
        plot_line_voltage()
    if make_vid == 'yes':
        video_writer()
    
    
    return None
##############################################################################
def video_writer():
    """video_writer() is a function that writes a video of the current and 
    voltage on the transmision line. Its reads the sotred values from the 
    simualtion from the V_line_matrix and I_line_matrix data arrays to 
    produce a matplotlib frame for each time step. The video_writer() function
    is called from within the function post_processing() if the videowrite 
    is true."""
    
    if(load_type==0.0):   # R Load
        title='R Load:' + ' RL= ' +str(RL) + ' \u03A9 '  
    
    elif(load_type==1.0):   #Parallel R-C Load
        title=('Parallel RC Load:' + ' CL= ' +str(Cl/(10**-12)) + ' pF'  + ' RL= ' +str(RL) + ' \u03A9')
        
    elif(load_type==2.0):   #Series R-L load
        title=('Series RL Load:' + ' LL= ' +str(Ll/(10**-9)) + ' nH'  + ' RL= ' +str(RL) + ' \u03A9')
    
    fig, ax =plt.subplots()
    fig.suptitle('Uniform Transmission Line Voltage' + ' Zo=' + str(Zo)+ ' \u03A9')
    ax.grid(linestyle='--')

    ax.set_ylim(-Vs_amp,+Vs_amp)
    ax.set_xlim(0,d)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('Line Voltage (V)')
    ax.set_title(title + ' Rg=' + str(Rg)+ ' \u03A9 ')
   
    # ims is a list of lists, each row is a list of artists to draw in the
    # current frame; here we are just animating one artist, the image, in
    # each frame
    ims = []
   
    for i in range(0,nt):  
        time=(i*dt)/(10**-9)
        im1=ax.text(0.06*d,0.77* Vs_amp, '%.2f' %time + ' ns')
        im2,=ax.plot(X1[0],V_line_matrix[i,:],label='Line Voltage',color='c',linewidth=2)
        
        print('\r Video Writing ' + str(((i+1)/nt)*100) + ' complete \r ',end='')
        ims.append([im1,im2])
       
    ani = animation.ArtistAnimation(fig, ims, interval=100, blit=True,repeat_delay=1000)
    writer = animation.FFMpegWriter(fps=20, metadata=dict(artist='Me'), bitrate=18000)
    
    filename="FDTD_Simulation_Video.mp4"
    
    ani.save(filename, writer=writer)
    
    plt.show()
    
    return None
##############################################################################
##############################################################################   
def main():
    
    '''Readin data for input PRogram Control'''
    readin_parameter_data_file()
    
    '''Setup the simulation: Initiliaze Data Arrays and Constants'''
    setup_FDTD_sim()
    
    #TimeAdvance()
    tpulse=V_source_signature_setup('gauss',nt,dt,Vs_amp)
    TimeAdvance(tpulse)


    '''Perform post_proccessing of data to produce plots and video'''
    post_processing()
 

    return None
##############################################################################
##############################################################################

main()

##############################################################################
##############################################################################
##############################################################################


















