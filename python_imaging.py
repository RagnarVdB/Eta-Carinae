######################################################
#
# FUNCTIONS TO RUN BEFORE WE START
#
#######################################################

# Function for reading the parameter file to check that you have the right "convert_type" set (vtuCCmpi) in your run
# and to work out the size of the blocks in the experiment
#
# you may need to apt install some of the modules, e.g. sudo apt install numpy
#
# INPUTS:
#
# parfile:    the name of the parameter file you want to read
#
# RETURNS:
#
# block_nx1, block_nx1:  standard experiment block dimensions

def readparfile(parfile):
    import os
    # need to run a convert with type \'vtuCCmpi\' first
    # read par file and get block information
    f = open(parfile, "r")
    read = f.read()
    for line in read.splitlines():
        if ('convert_type=' in line):
            data_type=line.split('=',1)[1]
            # in vtuCCmpi, the cell center values of variables are outputed, but
            # the provided coordinates are the locations of cell corners. It means
            # that the number of coordinates is greater than that of cells. The data
            # is outputed block by block, firstly cell values and then cell corner 
            # coordinates.
            if (data_type!='\'vtuCCmpi\''):
                print('ERROR! The convert_type is not vtuCCmpi, found: '+data_type)
                print(os.getcwd())
                print(parfile)
                sys.exit()
            
        if ('block_nx1=' in line):
            block_nx1=int(line.split('=',1)[1])
        
        if ('block_nx2=' in line):
            block_nx2=int(line.split('=',1)[1])
    
    f.close()
    return block_nx1, block_nx2


# Function to read the .vtu files
#
# INPUTS:
#
# Filebase:   the base name of the experiments
# Snapno:     the file number you want to read
# block_nx1, block_nx2: the block numbers define the sizes of the block you want to read (can be done with readparfile)
# namevars:   list of the variables you want read
# varunits:   units values of the variables you want to read
# lengthunit: length units, 2D array as this is set up for 2D experiments
# folder:     an optional subfolder
# fillwidth:  digits used in the save file numbering system, defauly value of 4
#
# RETURNS
#
# var:        [numgrid,numofvariables,block_nx1, block_nx2]: variable values
# x_corner, y_corner: cell coordinates of the corners of each block 
# x_center, y_center: cell coordinates of the centres of each block
# numgrid     total number of blocks in the snapshot
def readvtu(filebase, snapno, block_nx1, block_nx2, namevars, varunits, lengthunit, folder='', fillwidth=4):
    import numpy as np
    import math
    import os
    if folder!='':
        currentdir=os.getcwd()
        os.chdir(folder)
    
    numvar=len(namevars)     # number of variables
    x_corner=np.zeros((1,block_nx1+1))      # cell corner coordinates will be stored here block by block
    y_corner=np.zeros((1,block_nx2+1))
    x_block=np.zeros(block_nx1+1)     # array for storing cell corner coordinates in a block
    y_block=np.zeros(block_nx2+1)
    var=np.zeros((1,numvar,block_nx1,block_nx2))      # variables will be stored in this array block by block
    var_block=np.zeros((numvar,block_nx1,block_nx2))    # array for storing variables in a block
    numPoint=(block_nx1+1)*(block_nx2+1)    # number of points (cell corners) in each block 
    len_line=200      # maximum number of cell values in each line, see convert.t of amrvac
    numCell=block_nx1*block_nx2     # number of cells in each block
    numline=math.ceil(1.0*numCell/len_line)     # number of lines needed to out a variable in one block
    temp=np.zeros(numCell)
    fill='0'
    snaptempstr=f'{snapno:{fill}{fillwidth}}'
    vtufile=filebase+snaptempstr+'.vtu'
    # read vtu file
    print('reading {} ...'.format(vtufile))
    f=open(vtufile, 'r')
    
    # time
    for i in range(0,5):
        f.readline()
    
    sss=f.readline()
    time=float(sss)     # read time
    
    numgrid=0
    sss=f.readline()
    while (sss!=''):
        # read block variable data
        for ivar in range(0,numvar):
            # e.g. when 'Name=rho' appear in a line, following one/two lines
            # will be the data of rho (in a block)
            namestr='Name="{}"'.format(namevars[ivar])  
            if (sss.find(namestr)!=-1):
                icell=0
                temp[:]=0
                for iline in range(0,numline):
                    sss=f.readline()
                    lens=len(sss.split())
                    for ir in range(0,lens):
                        temp[icell]=sss.split()[ir]   # read data of a variable in the block
                        icell+=1
              
                icell=0
                for iy in range(0,block_nx2):
                    for ix in range(0,block_nx1):
                        var_block[ivar,ix,iy]=temp[icell]   # put variables in the same block together
                        icell+=1
    
        # read cell corner coordinate
        # when 'NumberOfComponents="3"' appear in a line, following 
        # (block_nx1+1)*(block_nx2+1) lines will provide the cell
        # corner coordinates
        namestr='NumberOfComponents="3"'
        if (sss.find(namestr)!=-1):
            for iy in range(0,block_nx2+1):
                for ix in range(0,block_nx1+1):
                    sss=f.readline()
                    if (iy==0):
                        x_block[ix]=sss.split()[0]  # get x of the block
                    
                    if (ix==0):
                        y_block[iy]=sss.split()[1]  # get y of the blocck
            
            if (numgrid==0):
                var[0,:,:,:]=var_block
                x_corner[0,:]=x_block
                y_corner[0,:]=y_block
            else:
                var=np.insert(var,numgrid,var_block,axis=0)   # add block data to the end of var
                x_corner=np.insert(x_corner,numgrid,x_block,axis=0)
                y_corner=np.insert(y_corner,numgrid,y_block,axis=0)
            
            numgrid+=1
        
        sss=f.readline()
    
    f.close()
    x_corner*=lengthunit[0]
    y_corner*=lengthunit[1]
    for loop2 in range(0,numvar):
        var[:,loop2,:,:]*=varunits[loop2]
    
    # calculate cell center coordinates
    x_center=np.zeros((numgrid,block_nx1))
    y_center=np.zeros((numgrid,block_nx2))
    for igrid in range(0,numgrid):
        for ix in range(0,block_nx1):
            x_center[igrid,ix]=0.5*(x_corner[igrid,ix]+x_corner[igrid,ix+1])
        
        for iy in range(0,block_nx2):
            y_center[igrid,iy]=0.5*(y_corner[igrid,iy]+y_corner[igrid,iy+1])
    
    if folder!='':
        os.chdir(currentdir)
    
    return var, x_corner, y_corner, x_center, y_center, numgrid



def getunits(filename, quiet=True):
    import os
    import numpy as np
    import matplotlib as plt
    FileRead = open(filename, 'r')
    Lines=FileRead.readlines()
    FileRead.close()
    linelooper=0
    strtemp=''
    while strtemp != "length":
        strtemp=Lines[linelooper]
        strtemp=strtemp[1:7]
        linelooper+=1
    
    unitstart=linelooper-1
    while strtemp != 'Reading':
        strtemp=Lines[linelooper]
        strtemp=strtemp[1:8]
        linelooper+=1
    
    unitend=linelooper-1
    unitnames=[]
    unitvals=np.zeros([unitend-unitstart])
    for looper in range(unitstart, unitend):
        strtemp=Lines[looper]
        strtemp=strtemp.replace("\n", "")
        strtemp=strtemp.replace(",", "")
        strtemp=strtemp.replace("  ", " ")
        strtemp=strtemp.replace("  ", " ")
        strtemp=strtemp.replace("  ", " ")
        strtemp=strtemp.replace("  ", " ")
        #print(strtemp)
        if strtemp[-1]==' ':
            strtemp=strtemp[:-1]
        if strtemp[0]==' ':
            strtemp=strtemp[1:]
        
        #print(strtemp)
        getline=np.array(strtemp.split(' '))
        if quiet==False:
            print(getline)
        
        unitvals[looper-unitstart]=float(getline[1])
        unitnames.append(getline[0])
   
    return unitnames,unitvals

def getunit(var, unitfile):
    looper=0
    unitnames, unitvals=getunits(unitfile)
    while unitnames[looper] != var:
        looper+=1
    
    return unitvals[looper]


###################################################################################
#
# Loading and plotting example
#
###################################################################################
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap
import numpy as np
from subprocess import run

# Note Te is not always saved, so you may want to edit it out for your example
#namevars=['rho','p','Te','v1','v2']   # variables we want to read
#rho_=0                                # specify index for variables
#p_=1
#Te_=2
#v1_=3
#v2_=4
namevars=['rho','p','v1','v2']   # variables we want to read
rho_=0                                # specify index for variables
p_=1
v1_=2
v2_=3
varunits=[1,1,1,1]                    # units of variables
lengthunit=[1,1]                      # length units in 2D experiment



def create_video(snapNumbers, parfilename, filebase, foldername, framerate, video_file):

    for snapno in snapNumbers:

        #
        # File info
        #
        # These are the files I checked the procedures with, they should work for your if you used the right 
        # convert_type in your run
        #
        # tips: 
        # 1) in your run folder create an "output" directory
        # 2) in the amrvac.par file add output to the base_filename, e.g. base_filename='output/wind_2d_eta_'
        # 3) I ran the simulation with 
        #     mpirun -np 16 ./amrvac -i amrvac_carinae.par | tee output/terminaloutput.txt
        # -i specifies the name of the parameter file to usee
        # | tee output/your_terminal_output_filename.txt stores the terminal output so that you can collect the unit length info

        # If you store the terminal output and your run prints the unit values then you can load it with this routine.
        # If not, just comment it out and use the unit numbers above!
        # This will convert the experiment into SI or CGS units depending on what was in the output file parameters
        unitfilename=foldername+"terminaloutput.txt"
        # unitnames, unitvals = getunits(unitfilename)
        # unitlength=getunit("length",unitfilename)
        # unitrho=getunit("rho",unitfilename)
        # unitp=getunit("p",unitfilename)
        # unitv=getunit("v",unitfilename)
        # varunits=[unitrho,unitp,unitv,unitv]   # units of variables
        # lengthunit=[unitlength,unitlength]           # length units in 2D experiment

        varunits=[1,1,1,1]
        lengthunit=[1,1]

        block_nx1, block_nx2=readparfile(parfilename)

        var, x_corner, y_corner, x_centre, y_centre, numgrid = readvtu(filebase, snapno, block_nx1, block_nx2, namevars, varunits, lengthunit, folder=foldername, fillwidth=4)

        #check some values to see if it has worked
        np.min(var[:,v2_,:,:])
        np.max(var[:,v2_,:,:])
        np.min(var[:,p_,:,:])
        np.max(var[:,p_,:,:])

        #plotting exaple
        plotlimsx=(np.min(x_corner),np.max(x_corner))
        plotlimsy=(np.min(y_corner),np.max(y_corner))
        textUL=(plotlimsx[0]+0.1*(plotlimsx[1]-plotlimsx[0]), plotlimsy[1]-0.1*(plotlimsy[1]-plotlimsy[0]))
        textUR=(plotlimsx[0]+0.6*(plotlimsx[1]-plotlimsx[0]), plotlimsy[1]-0.1*(plotlimsy[1]-plotlimsy[0]))

        fig_size=(8.8,8.8)
        ## if you want to plot multiple subplots in the same figure change the 1,1 inputs to 
        ## e.g. 3,2 for 3 cols and 2 rows of panels
        fig, axs = plt.subplots(1,1, constrained_layout=True, figsize=fig_size)
        ## specify the axes of the subplot if using subplots... don't if you aren't!
        #ax=axs[1,1]
        ax=axs
        varno=rho_
        myvar="rho"
        minplot=np.min(var[:,varno,:,:])
        maxplot=np.max(var[:,varno,:,:])
        varlabel='N$_H$ (cm$^{-3}$)'
        findit=np.where( np.array(namevars) == myvar )
        findit=findit[0]
        findit=findit[0]
        # plot data block by block. 
        for igrid in range(0,numgrid):
            pcm=ax.pcolormesh(x_corner[igrid,:],y_corner[igrid,:], var[igrid,findit,:,:].transpose(), norm=colors.LogNorm(vmin=minplot,vmax=maxplot), antialiased=False)

        ax.set_xlim(plotlimsx)
        ax.set_ylim(plotlimsy)
        cb=fig.colorbar(pcm, ax=ax)
        cb.set_label(varlabel)


        fillwidth=4
        fill='0'
        snaptempstr=f'{snapno:{fill}{fillwidth}}'

        # Save the figure or show it to screen, I just saved it.
        # If you save a load of them you can make a video with "ffmpeg" I can show you how.
        plt.savefig(foldername+filebase+snaptempstr+myvar+".png",dpi=300)
        #plt.show()
        plt.close()

    # Animate
    run(f"cd {foldername} && ffmpeg -framerate {framerate} -pattern_type glob -i '*.png'   -c:v libx264 -pix_fmt yuv420p {video_file}", shell=True)

if __name__ == "__main__":
    snapNumbers = [0, 1, 2]
    parfilename = "amrvac_carinae_double.par"
    filebase = "wind_2d_eta_"
    foldername = "output_double"
    framerate = 1
    video_file = "out.mp4"
    create_video(snapNumbers, parfilename, filebase, foldername, framerate, video_file)