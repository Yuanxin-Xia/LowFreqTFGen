# Transfer function generator

Install Transfer function generator in DTU high performance computer system

Implementation of COMSOL MATLAB LiveLink

This version is developed for random sources and receivers; contact Author if you need to specify source and receiver.

contact: xiayuanxin98@gmail.com | s213256@dtu.dk

OUTPUT: TF.txt

OUTPUT format <br>
 ```
Lx1 Ly1 Lz1 ... Transfer functions in defined frequency range/ interval 
...
LxN LyN LzN ... Transfer functions in defined frequency range/ interval 

 ```
## INSTALLATION

### To access DTUs HPC system, log in using the command (see section about HPC system as well) <br>
`> ssh username@login1.gbar.dtu.dk`
### Clone the code from git in the home user directory ~/ <br>
 `> git clone https://github.com/Spoke-arrows/LowFrecTFGen`
 
## RUN the simulation
### Make sure enter the folder   <br>
`bsub < RUN.lsf`
 
 
 # SIMULATION SETTING
### Frequency range   <br>
 ```
UFreq = 354; % Upper frequency
LFreq = 1;  % Lower frequency
IFreq = 1;  % Frequency interval
 ```
 
### Room size
In the MATLAB file MAIN.m, you can change room size  <br>

```
Room = [

Lx1 Ly1 Lz1

Lx2 Ly2 Lz2

...

LxN LyN LzN

]
```
### Impedance of different walls
Copy your frequency dependent impedance to the ImpedanceFile folder and also edit impedance in the MAIN.m
<br>
```
S1 = 'A';
S2 = 'R';
S3 = 'R';
S4 = 'R';
S5 = 'R';
S6 = 'R';
```
![image](https://user-images.githubusercontent.com/42115062/185906660-15faf9ab-1471-4346-80ee-5991430e6736.png)

There are three kinds of impedance you can add in one room, if you want more impedance for one room, contact Author

### SOURCE AND RECEIVER position

![image](https://user-images.githubusercontent.com/42115062/185909373-4a152752-85cf-439a-9159-88419b34d8f6.png)

The Point_generation.m will generate source and receiver points, in order to make the source and receiver enenly distributed in the room
It is defined in the Point_generation.m that how many block/grid in X,Y,Z direction
<br>
```
%% Source grid
nx = 3;
ny = 1;
nz = 2;
```

```
%% Receiver grid
nx = 6;
ny = 6;
nz = 6;
```
In this case, the total source amount is 6, and total receiver amount is 216.
The script will loop each room, and each source in the room and each receiver.

## Useful command

### Run time and effciency
bstat –C

### Job information in detail
bjobs –l

### Check the HPC available CPU (You can also change CPU in RUN.lsf)

nodestat -F hpc
