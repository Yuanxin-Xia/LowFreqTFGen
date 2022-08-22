# TFs

Install Transfer function generator in DTU high performance computer system

Implementation of COMSOL MATLAB LiveLink

contact: xiayuanxin98@gmail.com | s213256@dtu.dk

## INSTALLATION

### To access DTUs HPC system, log in using the command (see section about HPC system as well) <br>
`> ssh username@login1.gbar.dtu.dk`
### Clone the code from git in the home user directory ~/ <br>
 `> git clone https://github.com/Spoke-arrows/TFs`
 
## RUN the simulation
### Make sure enter the folder   <br>
`bsub < RUN.lsf`
 
## Room size
In the MATLAB file MAIN.m, you can change room size  <br>

```
Room = [

Lx1 Ly1 Lz1

Lx2 Ly2 Lz2

...

LxN LyN LzN

]
```
## Impedance of different walls
Copy your frequency independent impedance to the ImpedanceFile folder and also edit impedance in the MAIN.m
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

## Useful command

### Run time and effciency
bstat –C

### Job information in detail
bjobs –l

### Check the HPC available CPU (You can also change CPU in RUN.lsf)

nodestat -F hpc
