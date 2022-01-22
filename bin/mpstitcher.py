#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This program uses a master slave approach to consume a queue 
of elaborations using terastitcher
Copyright (c) 2014:
Massimiliano Guarrasi (1), Giulio Iannello (2), Alessandro Bria (2)
(1): CINECA
(2): University Campus Bio-Medico of Rome
The program was made in the framework of the HUMAN BRAIN PROJECT.
All rights reserved.

EXAMPLE of usage (X is the major version, Y is the minor version, Z is the patch):
mpirun -np XX python parastitcherX.Y.Z.py -2 --projin=import_xml_file -projout=displacement_xml_file OTHER_OPTIONS

where:
- OTHER_OPTIONS is the list of are any other required option (possibly empty)
See terastitcher documentation for more details

*******************************
*        Change Log           *
*******************************

v2.0.1 2017-10-07
- revisted platform dependent instructions
"""

import os
import sys
import shutil
import time
import datetime
import operator
import subprocess

from math import *
from glob import glob
from pathlib import Path
import __future__
import multiprocessing
from collections import deque
from subprocess import *

"""
Set prefix='./' if your terastitcher executable is present in the running director.
If  the terastitcher executable is included in the $PATH dir select prefix=''
"""
if (Path(__file__).parent/"terastitcher").exists():
    prefix = mdprefix = str(Path(__file__).parent) + os.sep
elif "TERASTITCHER" in os.environ:
    prefix = mdprefix = os.path.dirname(os.environ["TERASTITCHER"]) + os.sep
else:
    prefix = os.path.dirname(subprocess.check_call(["which", "terastitcher"]))
#prefix = "C:\\Progra~1\\TeraSt~1.61\\bin\\"

###############################################
# functions
###############################################

def partition ( m, n, N ) : 
    """
    return the number of partitions along V and H, respectively that are optimal to
    partition a block of size m_V x n_H in at least P sub-blocks

    m: block size along V
    n: block size along H
    N:   number of required partitions

    return:
    p_m, p_n: the  number of partitions along V and H, respectively

    PRE:
    """
    m = float(m)
    n = float(n)
    N = float(N)
    p_m_min = sqrt(N*m/n)
    p_n_min = sqrt(N*n/m)
    c_min = m * p_n_min + n * p_m_min

    # test point p
    p_m = ceil(p_m_min)
    p_n = floor(p_n_min)
    if p_n * p_m < N : # p does no satisfy the conditions
        # test point p'
        p_m = floor(p_m_min)
        p_n = ceil(p_n_min)
        if p_n * p_m < N : # p' does no satisfy the conditions
            # set (p_m, p_n) = p''
            p_m = ceil(p_m_min)
    if 1 <= p_m and p_m <= m/2.0 and 1 <= p_n and p_n <= n/2.0 : # (p_m, p_n) is a candidate solution
        # save p_m0 and set the current cost
        p_m0 = p_m
        c = m * p_n + n * p_m
    else : # the problem does not have an acceptable solution
        return (-1,-1)

    # optimal cost is in interval [c_min, c]

    # look for a better solution by increasing p_m
    p_m_cur = p_m + 1 # test a new candidate p_m
    p_n_cur = (c - n * p_m_cur) / m
    while p_m_cur <= m/2.0 and (p_m_cur * p_n_cur) >= N : # there may be a better solution
        if 1 <= floor(p_n_cur) and (p_m_cur * floor(p_n_cur)) >= N : # this is such a better solution
            # update solution
            p_m = p_m_cur
            p_n = floor(p_n_cur)
            c   = m * p_n + n * p_m
        # look for a better solution further increasing p_m
        p_m_cur = p_m_cur + 1
        p_n_cur = (c - n * p_m_cur) / m

    # look for a better solution by decreasing p_m
    p_m_cur = p_m0 - 1 # start from the initial candidate p_m
    p_n_cur = (c - n * p_m_cur) / m
    while floor(p_n_cur) <= p_m_cur and (p_m_cur * p_n_cur) >= N : # there may be a better solution
        if floor(p_n_cur) <= n/2.0 and (p_m_cur * floor(p_n_cur)) >= N : # this is such a better solution
            # update solution
            p_m = p_m_cur
            p_n = floor(p_n_cur)
            c   = m * p_n + n * p_m
        # look for a better solution further decreasing p_m
        p_m_cur = p_m_cur - 1
        p_n_cur = (c - n * p_m_cur) / m

    return (int(p_m), int(p_n))


##################################################################################


def read_input(inputf,nline=0):
    """
    Reads the file included in inputf at least up to line number nline (if declared).

    """
    i = 0
    data = [] 
    f = open(inputf, 'r')
    for line in f:
        line = line.strip()
        l = line.split(' ', 1)
        data.append(l)
        if (nline != 0) and (i > nline):
            break
        i += 1
    f.close()
    return data


def extract_np(inputf):
    """
    extract the number of slices along z from the input xml file.

    """
    data = read_input(inputf, 8)
    for tmp_string in data :
        #tmp_string = data[8][1]
        if tmp_string[1].find('stack_slices="') != -1 :
            start = tmp_string[1].find('stack_slices="')+14
            break
    end = tmp_string[1].find('" />')
    sel_string = tmp_string[1][start:end]
    nslices = int(sel_string)
    
    start = tmp_string[1].find('stack_rows="')+12
    end = start + tmp_string[1][start:].find('"')
    sel_string = tmp_string[1][start:end]
    nrows = int(sel_string)
    
    start = tmp_string[1].find('stack_columns="')+15
    end = start + tmp_string[1][start:].find('"')
    sel_string = tmp_string[1][start:end]
    ncols = int(sel_string)
    
    return (nrows, ncols, nslices)


def find_last_slash(string):
    """
    Search for / in a string. If one or more / was found, divide the string in a list of two string:
    the first containf all the character at left of the last / (included),
    and the second contains the remanent part of the text.
    If no / was found, the first element of the list will be set to ''
    """
    len_string = len(string)
    check = 0
    index = []
    i = 0
    for chara in string:
        if chara == '/' or chara == '\\':
            index.append(i)
            check =1
        i += 1
    if check == 1:
        last_slash = max(index)
        output_string = [string[0:last_slash+1],string[last_slash+1:]]
    else:
        output_string = ['',string]

    return output_string


def extract_params():
    """
    extract parameter from line of commands.

    """
    params = (sys.argv)
    return params

def check_flag(params, string, delete):
    """
    check if a parameter was beeen declared in the line of commands and return the associated value.
    If delete is true the related string will be deleted
    If it is not present, return 0

    """
    i = 0
    value = None
    size = len(string)
    for line in params:
        tmp = line.find(string)

        if tmp != -1:
            start = tmp + size
            sel_string = line[start:]
            if delete :
                params.pop(i)
            value = sel_string
        i += 1
    return value

def add_chars(params):
    string =['volin_plugin=','imin_plugin']
    i= 0
    for line in params:
        for local_string in string:
            tmp = line.find(local_string)
            if tmp != -1:
                size = len(local_string)
                sel_string = line.split('=')
                input_string = sel_string[1]
                mod_string = '"'+input_string+'"'
                tot_mod_string = sel_string[0]+'='+mod_string
                params[i] = tot_mod_string
        i += 1
    return params



def score_function(params):
    """
    Assigns a score value with the formula 
    score = N_of_slices
    params is a dictionary {input_name : [Np,Ny]}
    returns a dictionary of the form  {input_name : score}
    """
    tmp_scores = {}
    for name, Nz in params.iteritems():
        score = Nz
        tmp_scores[name] =score

    scores = {}
    # apply formula 
    den = max(tmp_scores.values()) # calculates the denominator of the formula
    for name, score in tmp_scores.iteritems():
        score = score*100.0/den
        scores[name] = score
    return scores


def sort_elaborations(scores):
    """
    scores is a dictionary of the form  {input_name : score}
    returns a list of input_name sorted by score
    """
    sorted_list = sorted(scores.iteritems(), key=operator.itemgetter(1), reverse=True)
    return [name[0] for name in sorted_list] 


def sort_work(params,priority):
    """
    params is a dictionary of the form  {input_name : value}
    priority is the list of input_name ordered by score

    returns a dictionary as params but ordered by score
    """	
    sorted_dict = {}
    for index in priority:
        sorted_dict.update({index:params[index]})
    return sorted_dict


def pop_left(dictionary):
    """
    Cuts the first element of dictionary and returns its first element (key:value)
    """	
    if len(dictionary) > 0:
        keys = list(dictionary.keys())
        first_el ={ keys[0] : dictionary[keys[0]]}
        dictionary.pop(keys[0])
    else:
        first_el = None
    return first_el



def worker(input_file, output_file):
    """
    perform elaboration for each element of the queue.

    """

    t1 = time.time()
    execution_string = prefix+"terastitcher " + input_file + " > "+  output_file
    print(execution_string)
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Sleep da cancellare!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #time.sleep((5.+100./(input_file.keys()[0]+1.))/100.)
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    os.system(execution_string)
    t2 = time.time()
    print (" ---> processor has calculated for " , t2-t1)
    return input_file


def num_processors():
    return int(os.environ.get("N_PROCESSORS", multiprocessing.cpu_count()))

def master(queue):
    """
    dispatch the work among processors

    queue is a list of job input

    """
    # constants to use as tags in communications
    WORKTAG = 1   
    DIETAG = 2

    nprocs = num_processors()

    pool = multiprocessing.Pool(nprocs)
    closed = False
    try:
        futures = []
        keys = sorted(queue)
        for key in keys:
            input_file = queue[key]
            print("Dispatching %s" % input_file)
            output_file = "output_"+str(key)+".out"
            futures.append(pool.apply_async(worker, (input_file,output_file)))
        pool.close()
        closed = True
        for future in futures:
            future.get()
    except:
        if not closed:
            pool.terminate()
        raise
    finally:
        pool.join()
        
def do_additional_partition ( nprocs, nrows, ncols, n_ss ) :
    """
    All parameters should be float
    """

    print('do_additional_partition ( nprocs =', nprocs, ', nrows =', nrows, ', ncols =', ncols, ', nlayers =', n_ss, ')')

    if n_ss >= 2*nprocs : # dataset has not to be further partitioned
        return (1,1)

    elif floor(nrows/2)*floor(ncols/2)*n_ss < 2*nprocs :
        print("WARNING: not enough paritions for", nprocs, "processors")
        print("Dataset partitioned in", int(floor(nrows/2)), "x", int(floor(ncols/2)), "x", n_ss, "=", int(floor(nrows/2))*int(floor(ncols/2))*n_ss, "partitions")
        return (floor(nrows/2), floor(ncols/2))

    else :                 # dataset has to be further partitioned
        m = max(nrows,ncols)
        n = min(nrows,ncols)

        (p_m, p_n) = partition ( m, n, ceil(2*nprocs/n_ss) ) # n_ss partitions are already available along Z

        if nrows < ncols : # the greater dimensione is the second: p_m and p_n have to exchanged
            temp = p_m
            p_m  = p_n
            p_n  = temp

        return (p_m, p_n)

# 		#initialize global variables (tables)
# 		initTables(m,n)
# 		
# 		c = costo(m,n,n_ss,2*nprocs) # each tile is already partitioned in n_ss layers
# 						
# 		return partition(nrows,ncols,n_ss,2*nprocs) # each tile is already partitioned in n_ss layers




###############################################################################
#                                                                             #
# --------------------------------- main ------------------------------------ #       
#                                                                             #
###############################################################################

if __name__ == '__main__':

    nprocs = num_processors()
    # timing
    t1 = time.time()
    print ('*'*80)
    print (str(datetime.datetime.utcnow()), " -- Calculation started on ", num_processors(), " cores.")
    print ('*'*80)

    exit_msg = 0

    # get the input list
    params = extract_params()
    print(params)
    #params = str("terastitcher -2 --projin=/Users/iannello/Home/Windows/myTeraStitcher/TestData/unstitched_RGB_2D/xml_import.xml --projout=/Users/iannello/Home/Windows/myTeraStitcher/TestData/unstitched_RGB_2D/xml_compdispl.xml --subvoldim=200 --imin_channel=G").split()
    params.pop(0)
    execution_flag = params.pop(0)

    if execution_flag == '-2':
        print ("Alignment will be performed")
    #elif execution_flag == '-6':
    #	print ("Merging will be performed.")
    else:
        print (execution_flag,": Option not allowed! No operations will be performed.")
        exit_msg = -1
        sys.exit()

    # select the minimum number of slices
    if execution_flag == '-2':
        tmp = check_flag(params, 'subvoldim=',True)
        if tmp == None:
            nsubstring_min = 200
            print ("Number of substring was not declared. It will be set to",  nsubstring_min, "by default.")
        else:
            tmp = int(tmp)
            nsubstring_min = tmp

    # Declare input file
    tmp = check_flag(params, 'projin=',False)
    if tmp == None:
        if execution_flag == '-2':
            input_name = 'xml_import.xml'
        print ("Name of the input file was not declared. It will be set to",  input_name, "by default.")
    else:
        input_name = tmp

    # Find standard output file
    if execution_flag == '-2':
        tmp = check_flag(params, 'projout=',True)
        if tmp == None:
            output_name = 'xml_compdispl.xml'
            print ("name of the output file was not declared. It will be set to",  output_name, "by default.")
        else:
            output_name = tmp
        len_out = len(output_name)
        output_name = output_name[0:len_out - 4]

    params = add_chars(params)

    # Reads the size of the tile matrix and the number slices
    (nrows, ncols, nslices) = extract_np(input_name)

    # Calculate the size of slices per task
    n_ss = int(ceil(float(nslices) / float(nsubstring_min))) # Number of substacks
    last_size = int(floor(float(nslices) / float(n_ss)))
    first_size = last_size + 1
    n_of_first_stacks = nslices % n_ss
    n_of_last_stacks = n_ss - n_of_first_stacks

    # look for further partitions
# 		nprocs = 6
# 		nrows = 11
# 		ncols = 12
    (p_r,p_c) = do_additional_partition(float(nprocs),float(nrows),float(ncols),float(n_ss)) # master should not be counted
    # compute start, end of rows partitions
    s_r = int(floor(nrows / p_r))
    r_r = nrows % p_r
    r_start = [0]
    r_end   = []
    for i in range(1,int(p_r - r_r)) :
        r_end.append(i*s_r)
        r_start.append(i*s_r)
    offset = (p_r - r_r) * s_r
    for i in range(int(r_r)) :
        r_end.append(offset + i*s_r)
        r_start.append(offset + i*s_r)
    r_end.append(nrows - 1)
    print(r_start, r_end)
    # compute start, end of columns partitions
    s_c = int(floor(ncols / p_c))
    r_c = ncols % p_c
    c_start = [0]
    c_end   = []
    for i in range(1,int(p_c - r_c)) :
        c_end.append(i*s_c)
        c_start.append(i*s_c)
    offset = (p_c - r_c) * s_c
    for i in range(int(r_c)) :
        c_end.append(offset + i*s_c)
        c_start.append(offset + i*s_c)
    c_end.append(ncols - 1)
    print(c_start, c_end)
    # generate temporary directory names
    if len(r_start)>1 or len(c_start)>1 :
        gr_dir_names = []
        for i in range(len(r_start)) :
            for j in range(len(c_start)) :
                gr_dir_names.append('gr_R['+str(r_start[i])+','+str(r_end[i])+']_C['+str(c_start[j])+','+str(c_end[j])+']/')

    # Create some dictionaries containing:
    # - the number of strips per substak (new_params),
    # - the index of the starting strip for each stack (start_dict)
    # - the index of the ending strip for each stack (start_dict)
    # - the entire command line for each stacks
    cmd_string = {}
    params_str = ' '.join(params)
    #if execution_flag == '-2': # it has already been checked

    tmp_out_name = find_last_slash(output_name)

    # Create temporary directory for temporary xml files
    tmp_xml_dir = tmp_out_name[0]+'tmp/'
    #execution_string = 'mkdir -p ' + tmp_xml_dir
    #print (execution_string)
    #os.system(execution_string)
    try :
        shutil.rmtree(tmp_xml_dir)
    except OSError :
        pass
    print('removed directory tree ' + tmp_xml_dir + ' if any')
    os.mkdir(tmp_xml_dir)
    print('created directory ' + tmp_xml_dir)
    #execution_string = 'rm -f '+tmp_xml_dir+'*'
    #print (execution_string)
    #os.system(execution_string)

    if len(r_start) == 1 and len(c_start) == 1 : # no groups

        # generate commands if there are no groups
        start_dict = {}
        end_dict = {}
        start_tmp = 0
        end_tmp = 0
        new_params = {}
        new_output_name = tmp_xml_dir+tmp_out_name[1]
        for i in range(n_ss):
            start_dict.update({i:end_tmp})
            if i < n_of_first_stacks:
                new_params.update({i:first_size})
                end_tmp += first_size
            else:
                new_params.update({i:last_size})
                end_tmp += last_size
            end_dict.update({i:end_tmp-1})
            tmp_string = execution_flag+' '+params_str+' --projout='+new_output_name+'-'+str(start_dict[i]).zfill(6)+'-'+str(end_dict[i]).zfill(6)+'.xml'+' --subvoldim='+str(new_params[i])+' --D0='+str(start_dict[i])+' --D1='+str(end_dict[i])+' --noprogressbar'
            cmd_string.update({i:tmp_string})

    else : # there are groups

        # create additional temporary directories
        for r in range(len(r_start)) :
            for c in range(len(c_start)) :
                # Create temporary directory for temporary xml files
                gr_xml_dir = tmp_out_name[0]+gr_dir_names[r*len(c_start)+c]
                #execution_string = 'mkdir -p ' + gr_xml_dir
                #print (execution_string)
                #os.system(execution_string)
                try :
                    shutil.rmtree(gr_xml_dir)
                except OSError :
                    pass
                print('removed directory tree ' + gr_xml_dir + ' if any')
                os.mkdir(gr_xml_dir)
                print('created directory ' + gr_xml_dir)
                #execution_string = 'rm -f '+gr_xml_dir+'*'
                #print (execution_string)
                #os.system(execution_string)

                # generate commands for this group
                start_dict = {}
                end_dict = {}
                start_tmp = 0
                end_tmp = 0
                new_params = {}
                new_output_name = gr_xml_dir+tmp_out_name[1]
                for i in range(n_ss):
                    start_dict.update({i:end_tmp})
                    if i < n_of_first_stacks:
                        new_params.update({i:first_size})
                        end_tmp += first_size
                    else:
                        new_params.update({i:last_size})
                        end_tmp += last_size
                    end_dict.update({i:end_tmp-1})
                    tmp_string = execution_flag+' '+params_str+' --projout='+new_output_name+'-'+str(start_dict[i]).zfill(6)+'-'+str(end_dict[i]).zfill(6)+'.xml'+' --subvoldim='+str(new_params[i]).zfill(6)+' --D0='+str(start_dict[i]).zfill(6)+' --D1='+str(end_dict[i]).zfill(6)
                    tmp_string = tmp_string+' --R0='+str(r_start[r])+' --R1='+str(r_end[r])+' --C0='+str(c_start[c])+' --C1='+str(c_end[c])+' --noprogressbar --parallel'
                    if c < (len(c_start) - 1) : # it not the last group of columns: disable displacement computation of last column
                        tmp_string = tmp_string+' --disable_last_col'
                    if r < (len(r_start) - 1) : # it not the last group of columns: disable displacement computation of last column
                        tmp_string = tmp_string+' --disable_last_row'
                    cmd_string.update({(r*len(c_start)+c)*n_ss+i:tmp_string})

    ## start scoring
    #scores = score_function(new_params)

    # Sort tasks by score
    #elaborations = sort_elaborations(scores)
    #work_list = sort_work(cmd_string,elaborations)
    work_list = cmd_string

    # Call the routine to distribute the work between cpus
    master(work_list)

    if len(r_start) == 1 and len(c_start) == 1 : # no groups

        #Collect all the xml files corresponding to subvolumes in a unique xml file (only for step 2)
        slash_pos = len(tmp_xml_dir) - 1
        execution_string = mdprefix+'mergedisplacements -d="'+tmp_xml_dir[0:slash_pos]+'" -o="'+tmp_out_name[0]+tmp_out_name[1]+'.xml"'
        print (execution_string)
        os.system(execution_string)

    else : # there are groups

        #Collect xml files of each group in a single xml file
        for dname in gr_dir_names :
            slash_pos = dname.find('/')
            suffix = dname[dname.find('_'):slash_pos]
            execution_string = mdprefix+'mergedisplacements -d="'+tmp_out_name[0]+dname[0:slash_pos]+'" -o="'+tmp_out_name[0]+'tmp/'+tmp_out_name[1]+suffix+'.xml"'
            print (execution_string)
            os.system(execution_string)

        #Collect all the xml files corresponding to groups in a unique xml file (only for step 2)
        slash_pos = len(tmp_xml_dir) - 1
        execution_string = mdprefix+'mergedisplacements --mgroups -d="'+tmp_xml_dir[0:slash_pos]+'" -o="'+tmp_out_name[0]+tmp_out_name[1]+'.xml"'
        print (execution_string)
        os.system(execution_string)

    #clean up all temporary directories and files
    if len(r_start)>1 or len(c_start)>1 :
        for r in range(len(r_start)) :
            for c in range(len(c_start)) :
                gr_xml_dir = tmp_out_name[0]+gr_dir_names[r*len(c_start)+c]
                try:
                    shutil.rmtree(gr_xml_dir)
                except OSError :
                    pass
                print( 'deleted directory ' + gr_xml_dir + ' and all files in it')
    try:
        shutil.rmtree(tmp_xml_dir)
    except OSError :
        pass
    print( 'deleted directory ' + tmp_xml_dir + ' and all files in it')
    t2 = time.time()
    print ('*'*80)
    print (str(datetime.datetime.utcnow()), "-- Calculation ended after ", (t2-t1), " seconds")
    print ('*'*80)
