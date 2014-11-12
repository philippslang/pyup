# screen and file output model
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
import numpy as np
import math
from itertools import cycle
import sys

def prompt_for_exit():
    print ''
    raw_input('done, hit key to exit...')
    sys.exit()

def screen_out(result_data):
    for i in range(len(result_data)):
        x0, y0, x1, y1, kface_w_mean, kface_w_hmean, kxx, kxy, kyy, kmin, kmax, eigw, eigv = result_data[i]
        if i == 0:
            cwidth = 9
            print ''
            print str('window').rjust(cwidth), str('x0').rjust(cwidth), str('y0').rjust(cwidth), str('x1').rjust(cwidth), str('y1').rjust(cwidth), str('kmean').rjust(cwidth), str('khmean').rjust(cwidth), str('kxx').rjust(cwidth), str('kxy').rjust(cwidth), str('kyy').rjust(cwidth), str('kmin').rjust(cwidth), str('kmax').rjust(cwidth)
        print '{0:9d} {1:9.2f} {2:9.2f} {3:9.2f} {4:9.2f} {5:9.2e} {6:9.2e} {7:9.2e} {8:9.2e} {9:9.2e} {10:9.2e} {11:9.2e}'.format(i, x0, y0, x1, y1, kface_w_mean, kface_w_hmean, kxx, kxy, kyy, kmin, kmax )


def file_out( result_data, fname='keff_windows.txt' ):
    result_lines = ['window\tx0\ty0\tx1\ty1\tkmean\tkhmean\tkxx\tkxy\tkyy\tkmin\tkmax\n']
    for i in range(len(result_data)):
        x0, y0, x1, y1, kface_w_mean, kface_w_hmean, kxx, kxy, kyy, kmin, kmax, eigw, eigv = result_data[i]
        result_lines += [ str(str(i) + '\t' + str(x0) + '\t' + str(y0) + '\t' + str(x1) + '\t' + str(y1) + '\t' + str(kface_w_mean) + '\t' + str(kface_w_hmean) + '\t'+ str(kxx) + '\t' + str(kxy) + '\t' + str(kyy) + '\t' + str(kmin) + '\t' + str(kmax) +'\n') ]
    with open( fname, 'w' ) as f:
        f.writelines(result_lines)


def ellipse_plot(result_data,xmax,ymax,xmin=0.,ymin=0.):
    plt.clf()
    a = plt.subplot(111, aspect='equal')

    h_ell_max = (xmax-xmin)/7
    kmax_all = 0.
    for i in range(len(result_data)):
        x0, y0, x1, y1, kface_w_mean, kface_w_hmean, kxx, kxy, kyy, kmin, kmax, eigw, eigv = result_data[i]
        kmax_all = max( kmax_all, kmax )

    color_sim = cycle(['red','blue','green','magenta'])
    text_offset = (x1-x0)/15
    for i in range(len(result_data)):
        color = next(color_sim)
        x0, y0, x1, y1, kface_w_mean, kface_w_hmean, kxx, kxy, kyy, kmin, kmax, eigw, eigv = result_data[i]
        kratio = kmin/kmax
        h_ell = h_ell_max * kmax/kmax_all
        w_ell = h_ell*kratio
        wmax = np.argmax(eigw)
        vmax = np.array(eigv[:,wmax])
        angledeg = math.acos(np.dot( vmax, [1.,0.]))*180/np.pi
        xcenter, ycenter = x0+(x1-x0)/2, y0+(y1-y0)/2
        a.scatter(xcenter,ycenter,marker='.',s=30,color='black')
        a.text( x0+text_offset, y0+text_offset, str(i), color=color )
        a.add_artist( Ellipse((xcenter, ycenter), h_ell, w_ell, angledeg, fill=False, color=color ) )    
        a.add_artist( Rectangle((x0,y0),(x1-x0),(y1-y0),fill=False,linestyle='dashed',color=color ) )

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.savefig('keff_windows.png',dpi=300)
    plt.close()
