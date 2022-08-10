
FDTD transmison line simulator:
    This software is capable of simulating the time_dependent transmission line 
    currents and voltages on a uniform lossy transmision line by using second 
    order accurate FDTD approximation to the transmision line equations.

Language: Python 

scriptname: Txsim.py

    
Written by: Erich Wanzek
Written 2/17/2021
Intro to Computational Electromangtics
University of Colrado Denver


To run properly, the libraries:

numpy
matplotlib

must be installed.

The media package ffmpeg is included to perfrom video writing. 


This code was developed in and performed within the spyder IDE for python. 
The spyder IDE has the requried libraries already available with its installation.


To input data into the software, there is a txt file that contains all the parameters and input controlls.
This txt file is called Tx_Line_simulation_parameters.txt
This file must be in the same path folder as the python script Txsim.py so it can be read in by the Txsim.py script.


Here is the layout of the input txt file, with example data:

L=250e-9                          %Unit per length inductance
C=100e-12                         %unit per length capacitance
R=0                               %unit per length resitance 
G=0                               %unit per length conductance
d=0.5                             %Transmission line length
ncells=50                         %number of spatial cells
cfln=1                            %cfln number
simtime=5e-9                      %simulation duration
Rg=50                             %source resitance
RL=50                             %load resitance
source_signature=trapezoidal      %source signature type, either gausssian or trapezoidal
Vs_amplitude=2                    %voltage source amplitude
pulse_delay=1000e-12              % pulse delay time for gaussian pulse
pulse_half_width=50e-11           %pulse half width for gaussian pulse
risetime=200e-12                  %rise time for trapezoidal pulse
duration=500e-12                  %duration for trapezoidal pulse
falltime=200e-12                  %fall time for trapezoidal pulse
load_type=0                       %Load type: 0 for purely R, 1 for parallel RC, 2 for series RL load
C_load=5e-12                      %Load capacitance
L_load=10e-9                      %load inductance
make_video=yes                    %option to write video: yes or no, Video is saved in the same folder as script as FDTD_simulation_video.mp4
graph_source_voltage=no           %option to plot the source voltage vs time: yes or no
graph_load_voltage=no		  %option to plot the load voltage vs time: yes or no
graph_line_voltage=yes            %option to graph line voltage at snap_shot_time: yes or no
snap_shot_time=1.250e-9		  %snap_shot time to plot the line voltage at














 