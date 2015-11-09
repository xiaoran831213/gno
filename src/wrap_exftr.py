import pdb
import os
import os.path as pt
from glob import glob as gg
import hlp
from itertools import izip

def write_script(dst, n = 50000):
    """
    """
    dst = pt.expandvars(dst)
    hlp.mk_dir(dst)

    ## surface region to be sampled per task run.
    cmd = './exftr.R --rng {f0},{f1} --wgs wgs --ped all.ped &>{i:03d}.log\n'
    s = 200 # step size
    for fo, i in hlp.hpcc_iter(
            xrange(0, n, s), dst, npb=1, mpn=8, tpp=2.0,
            mds=['R/3.2.0'],
            cpy=['exftr.R'],
            lnk=['.', '../wgs', '../dat/all.ped'],
            debug=False):

        ## save the working material specification for one nodes line
        fo.write(cmd.format(f0=i+1, f1=i+s, i=i/s))

def test():
    write_script('../dsg')
    
if __name__ == "__main__":
    pass
