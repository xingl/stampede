#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from subprocess import call
import os
from interp import *
from finite_differences import *
import sys

kx_center_scan = False
setup_global = True

x0_scan = False
setup_lilo = False

#####Modify
case_name = '02b'
case = '1120907032.01012'
efit_file_name = 'g_901_901_1415_imode'
x0_values=[0.965,0.975,0.985,0.99]
ky_scan_string = '0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.7, 1.0'
template_prob_dir_loc  = 'prob_template_local'
template_prob_dir_glob = 'prob_template_global'
submit_runs = True
x0_glob = '0.98'
lx_a = '0.024'
nx0_glob = 96
batch_script_job_prefix = '#SBATCH -J '
submit_command = 'sbatch'
ExB_glob = 0.0   #'scan' for both 0.0 and -1111.0
### For kx_center scan
num_kxcenter = 5
gene_dirname = 'gene-dev_xing/'
######Modify

#####Setup
#####Setup
#####Setup
x0_scan_string = str(x0_values[0])
for i in range(1,len(x0_values)):
    x0_scan_string += ', '+str(x0_values[i])
print "x0_scan_string:",x0_scan_string
#basedir='/scratch/02289/xl3676/cmodlin/'
basedir='/scratch/02289/xl3676/cmod_global'
diagdir = basedir+case_name
homedir = '/home1/02289/xl3676/'
genedir = homedir + gene_dirname
probdirloc = genedir + 'prob_loc_' + case_name
probdirglob = genedir + 'prob_glob_' + case_name
probdirlilo = genedir + 'prob_lilo_' + case_name
probdirkxc = []
for i in range(len(x0_values)):
    probdirkxc.append(genedir+'prob_kxc_'+case_name+'_'+str(x0_values[i]))

eqs_dir = genedir+'cmod/'
print "Checking existence of efit file."
if os.path.isfile(eqs_dir+efit_file_name):
    print "Efit file exists:",efit_file_name
else:
    sys.exit("Efit file does not exist.  Select different efit file.") 
if os.path.isfile(eqs_dir+'profiles_wb'+case+'.iterdb'):
    print "Iterdb file exists:",'profiles_wb'+case+'.iterdb'
else:
    sys.exit("Iterdb file does not exist.  Select different iterdb file.") 

#####Setup
#####Setup
#####Setup

def calc_shat(q,rhot):
    rhot0 = np.arange(10000.0)/9999.0
    q0 = interp(rhot,q,rhot0)
    qprime = fd_d1_o4(q0,rhot0)
    shat = rhot0/q0*qprime
    return rhot0,q0,shat 

#####Execute
#####Execute
#####Execute
if not os.path.exists(diagdir):
    print diagdir," does not exist."  
    os.makedirs(diagdir)
    print "Now it does."

if x0_scan:
    if os.path.exists(probdirloc):
        call(['rm','-r',probdirloc])
    call(['cp','-r',template_prob_dir_loc,probdirloc])
    call(['cp','-r',homedir+'/scripts/setup_linear_runs.py',probdirloc])

if setup_global:
    if os.path.exists(probdirglob):
        call(['rm','-r',probdirglob])
    call(['cp','-r',template_prob_dir_glob,probdirglob])
    call(['cp','-r',homedir+'/scripts/setup_linear_runs.py',probdirglob])

if setup_lilo:
    if os.path.exists(probdirlilo):
        call(['rm','-r',probdirlilo])
    call(['cp','-r',template_prob_dir_glob,probdirlilo])
    call(['cp','-r',homedir+'/scripts/setup_linear_runs.py',probdirlilo])

if kx_center_scan:
    for dir in probdirkxc:
        if os.path.exists(dir):
            call(['rm','-r',dir])
        call(['cp','-r',template_prob_dir_loc,dir])
        call(['cp','-r',homedir+'/scripts/setup_linear_runs.py',dir])

#### LOCAL ####
if x0_scan:
    os.chdir(probdirloc)
    call(['ln','-s',eqs_dir+efit_file_name])
    call(['ln','-s',eqs_dir+'profiles_wb'+case+'.iterdb'])
    #call(['ln','-s',eqs_dir+'gene_profiles_e'])
    #call(['ln','-s',eqs_dir+'gene_profiles_i'])
    call(['ln','-s',eqs_dir+'q_profile.dat'])

    #change submit.cmd file
    f=open('submit.cmd','r')
    cmd_data=f.read()
    f.close()
    cmd_data_split = cmd_data.split('\n')
    cmd_data_split[1] = batch_script_job_prefix+case_name+'         # Job Name'
    cmd_file_out='\n'.join(cmd_data_split)
    f=open('submit.cmd','w')
    f.write(cmd_file_out)
    f.close()

    #change diagdir
    f=open('parameters','r')
    parfile=f.read()
    f.close()
    parfile_split = parfile.split('\n')
    for i in range(len(parfile_split)):
        if 'diagdir' in parfile_split[i]:
            parfile_split[i] = 'diagdir = \''+diagdir+'\''
        if 'x0' in parfile_split[i] and 'nx0' not in parfile_split[i]:
            parfile_split[i] = 'x0 = 0.95  !scanlist: '+x0_scan_string
        if 'kymin' in parfile_split[i]: 
            parfile_split[i] = 'kymin = 0.05  !scanlist: '+ky_scan_string
        if 'geomfile' in parfile_split[i]: 
            parfile_split[i] = 'geomfile = ' + '\''+efit_file_name + '\''
        if 'iterdb_file' in parfile_split[i]: 
            parfile_split[i] = 'iterdb_file = ' + '\''+'profiles_wb'+ case+'.iterdb' + '\''
    parfile_out='\n'.join(parfile_split)
    f=open('parameters','w')
    f.write(parfile_out)
    f.close()

    if submit_runs:
       call([submit_command,'submit.cmd',])

##### GLOBAL #####
if setup_global:
    os.chdir(probdirglob)
    call(['ln','-s',eqs_dir+efit_file_name])
    call(['ln','-s',eqs_dir+'profiles_wb'+case+'.iterdb'])
    #call(['ln','-s','gene_profiles_e'])
    #call(['ln','-s','gene_profiles_i'])


    #change submit.cmd file
    f=open('submit.cmd','r')
    cmd_data=f.read()
    f.close()
    cmd_data_split = cmd_data.split('\n')
    cmd_data_split[1] = batch_script_job_prefix+case_name+'_glob     # Job Name'
    cmd_file_out='\n'.join(cmd_data_split)
    f=open('submit.cmd','w')
    f.write(cmd_file_out)
    f.close()

    #change diagdir
    f=open('parameters','r')
    parfile=f.read()
    f.close()
    parfile_split = parfile.split('\n')
    for i in range(len(parfile_split)):
        if 'diagdir' in parfile_split[i]:
            parfile_split[i] = 'diagdir = \''+diagdir+'\''
        if 'x0' in parfile_split[i] and 'nx0' not in parfile_split[i]:
            parfile_split[i] = 'x0 = '+x0_glob
        if 'nx0' in parfile_split[i]:
            parfile_split[i] = 'nx0 = '+str(nx0_glob)
        if 'lx_a' in parfile_split[i]:
            parfile_split[i] = 'lx_a = '+lx_a
        if 'kymin' in parfile_split[i]: 
            parfile_split[i] = 'kymin = 0.05  !scanlist: '+ky_scan_string
        if 'geomfile' in parfile_split[i]: 
            parfile_split[i] = 'geomfile = ' + '\''+efit_file_name + '\''
        if 'iterdb_file' in parfile_split[i]: 
            parfile_split[i] = 'iterdb_file = ' + '\''+ 'profiles_wb'+case+'.iterdb' + '\''
        if 'ExBrate' in parfile_split[i] and ExB_glob: 
            if ExB_glob == 'scan':
                parfile_split[i] = 'ExBrate = -1111.0 !scanlist: 0.0, -1111.0' 
            else:
                parfile_split[i] = 'ExBrate = -1111.0' 
    parfile_out='\n'.join(parfile_split)
    f=open('parameters','w')
    f.write(parfile_out)
    f.close()

    if submit_runs:
       call([submit_command,'submit.cmd'])

##### LILO #####
if setup_lilo:
    os.chdir(probdirlilo)
    call(['ln','-s',eqs_dir+efit_file_name])
    call(['ln','-s',eqs_dir+'profiles_wb'+case+'.iterdb'])
    call(['ln','-s','gene_profiles_e'])
    call(['ln','-s','gene_profiles_i'])

    #change submit.cmd file
    f=open('submit.cmd','r')
    cmd_data=f.read()
    f.close()
    cmd_data_split = cmd_data.split('\n')
    cmd_data_split[1] = batch_script_job_prefix+case_name+'_lilo     # Job Name'
    cmd_file_out='\n'.join(cmd_data_split)
    f=open('submit.cmd','w')
    f.write(cmd_file_out)
    f.close()

    #change diagdir
    f=open('parameters','r')
    parfile=f.read()
    f.close()
    parfile_split = parfile.split('\n')
    for i in range(len(parfile_split)):
        if 'diagdir' in parfile_split[i]:
            parfile_split[i] = 'diagdir = \''+diagdir+'\''
        if 'x0' in parfile_split[i] and 'nx0' not in parfile_split[i]:
            parfile_split[i] = 'x0 = '+x0_glob
        if 'lx_a' in parfile_split[i]:
            parfile_split[i] = 'lx_a = '+lx_a
        if 'kymin' in parfile_split[i]: 
            parfile_split[i] = 'kymin = 0.05  !scanlist: '+ky_scan_string
        if 'geomfile' in parfile_split[i]: 
            parfile_split[i] = 'geomfile = ' + '\''+efit_file_name + '\''
        if 'iterdb_file' in parfile_split[i]: 
            parfile_split[i] = 'iterdb_file = ' + '\''+ 'profiles_wb'+case+'.iterdb' + '\''
        if 'ExBrate' in parfile_split[i] and ExB_glob: 
            if ExB_glob == 'scan':
                parfile_split[i] = 'ExBrate = -1111.0 !scanlist: 0.0, -1111.0' 
            else:
                parfile_split[i] = 'ExBrate = -1111.0' 
        if 'lilo' in parfile_split[i]: 
            parfile_split[i] = 'lilo = T' 
        if 'rad_bc_type' in parfile_split[i]: 
            parfile_split[i] = 'rad_bc_type = -1' 
        if 'drive_buffer' in parfile_split[i]: 
            parfile_split[i] = 'drive_buffer = F' 
        if 'mag_prof' in parfile_split[i]: 
            parfile_split[i] = 'mag_prof = F' 
    parfile_out='\n'.join(parfile_split)
    f=open('parameters','w')
    f.write(parfile_out)
    f.close()

    if submit_runs:
       call([submit_command,'submit.cmd'])

##### kx_center_scan #####
if kx_center_scan:
    call(['ln','-s',homedir+'scripts/scan_info.py',diagdir])
    call(['ln','-s',homedir+'scripts/plot_scan_info.py',diagdir])
    #call(['ln','-s',homedir+'scripts/plot_contour_x0_ky.py',diagdir])
    for j in range(len(x0_values)):
        os.chdir(probdirkxc[j])
        call(['ln','-s',eqs_dir+efit_file_name])
        call(['ln','-s',eqs_dir+'profiles_wb'+case+'.iterdb'])
        #call(['ln','-s',eqs_dir+'gene_profiles_e'])
        #call(['ln','-s',eqs_dir+'gene_profiles_i'])
        call(['ln','-s',eqs_dir+'q_profile.dat'])

        #Calculating shat
        rbs = np.genfromtxt('q_profile.dat')
        q = rbs[:,1]
        rhot = rbs[:,0]
        rhot0,q0,shat = calc_shat(q,rhot)
        xind = np.argmin(abs(rhot0-x0_values[j]))
        shat_out = abs(shat[xind])
        print "Assuming shat = ",shat_out
        if j==0:
            print "Saving shat.dat"
            np.savetxt('shat.dat',np.column_stack((rhot0,q0,shat)))

        #change submit.cmd file
        f=open('submit.cmd','r')
        cmd_data=f.read()
        f.close()
        cmd_data_split = cmd_data.split('\n')
        cmd_data_split[1] = batch_script_job_prefix+case_name+'kxc_x0'+str(x0_values[j])+'         # Job Name'
        cmd_file_out='\n'.join(cmd_data_split)
        f=open('submit.cmd','w')
        f.write(cmd_file_out)
        f.close()

        kx_center_scan_string = '  !scanlist: 0.0 '
        for i in range(num_kxcenter-1):
            kx_center_scan_string += ', '+str(0.8*(i+1)/float(num_kxcenter)*2*np.pi)+'*'+str(shat_out)+'*kymin(1)'

        print "kx_center_scan_string",kx_center_scan_string

        #change diagdir
        f=open('parameters','r')
        parfile=f.read()
        f.close()
        parfile_split = parfile.split('\n')
        for i in range(len(parfile_split)):
            if 'diagdir' in parfile_split[i]:
                parfile_split[i] = 'diagdir = \''+diagdir+'\''
            if 'x0' in parfile_split[i] and 'nx0' not in parfile_split[i]:
                parfile_split[i] = 'x0 = '+str(x0_values[j])
            if 'kymin' in parfile_split[i]: 
                parfile_split[i] = 'kymin = 0.05  !scanlist: '+ky_scan_string
            if 'geomfile' in parfile_split[i]: 
                parfile_split[i] = 'geomfile = ' + '\''+efit_file_name + '\''
            if 'iterdb_file' in parfile_split[i]: 
                parfile_split[i] = 'iterdb_file = ' + '\''+ 'profiles_wb'+case+'.iterdb' + '\''
            if 'kx_center' in parfile_split[i]: 
                parfile_split[i] = 'kx_center = 0.0  ' + kx_center_scan_string
        parfile_out='\n'.join(parfile_split)
        f=open('parameters','w')
        f.write(parfile_out)
        f.close()

        if submit_runs:
           call([submit_command,'submit.cmd',])

