# $Id$

import numpy as N
from pencil.files.dim import read_dim
from os import path

def read_power(file, datadir='data/'):
    """ 
    29-apr-2009/dintrans: coded
    t,dat=read_power(name_power_file)
    Read a power spectra file like 'data/poweru.dat'
    """ 
    filename = path.join(datadir, file)
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
#
#  find the number of blocks (t,power) that should be read
#
    dim=read_dim(datadir=datadir)
    nblock=int(len(lines)/int(N.ceil(dim.nxgrid/2/8.)+1))
#
    with open(filename, 'r') as infile:
        t=N.zeros(1, dtype='Float32')
        data=N.zeros(1, dtype='Float32')
        for i in range(nblock):
            st=infile.readline()
            t=N.append(t, float(st))
            for ii in range(int(N.ceil(dim.nxgrid/2/8.))):
                st=infile.readline()
                data=N.append(data, N.asarray(st.split()).astype('f'))

    t=t[1:] ; data=data[1:]
    nt=len(t) ; nk=int(len(data)/nt)
    data=data.reshape(nt, nk)
    return t, data

