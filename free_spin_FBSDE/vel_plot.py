import sys
import os

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def main(argv):
    # filepath = sys.argv[1]
    # if not os.path.isfile(filepath):
    #   print("File path {} does not exist. Exiting...".format(filepath))
    #   sys.exit()
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig = plt.gcf() 
    fig.set_size_inches(20,10)

    folder = "pltdata/"
    if not argv:
        file_filter = "_f.txt"
    elif argv[0] == "hyd":
        file_filter = "hydrogen_f.txt"
    elif argv[0] == "har":
        file_filter = "harm_osci_f.txt"
    elif argv[0] == "dou":
        file_filter = "double_well_f.txt"
    elif argv[0] == "dou_b":
        file_filter = "double_well_fb.txt"
    elif argv[0] == "har_b":
        file_filter = "harm_osci_fb.txt"
    elif argv[0] == "hyd_b":
        file_filter = "hydrogen_fb.txt"
    elif argv[0] == "tel":
        file_filter = "teller_f.txt"
    elif argv[0] == "tel_b":
        file_filter = "teller_fb.txt"
    elif argv[0] == "spin_b":
        file_filter = "spin_fb.txt"
    else:
        print("wrong option")
        sys.exit()
    if len(argv) > 1:
        if argv[1] == "num":
            numerov = True
    else:
        numerov = False
    a = []
    for filename in os.listdir(folder):
        #if filename.endswith(file_filter):
        if filename.endswith(file_filter):# or filename.endswith("well_f.txt"):
            data = pd.read_csv(folder+filename, sep=" ", header=None)
            data.columns = ["x", "p", "u", "q", "index", "newl"]
            #ax1.plot(data.x,data.u*np.cos(0.5*data.x)/np.sin(0.5*data.x), label=filename.replace("_"+file_filter,""))
            
            if file_filter == "spin_fb.txt":
                ax1.plot(data.x,data.u*np.cos(data.x/2)/np.sin(data.x/2), "o", label= filename.replace("_"+file_filter,"").replace("m_",r'$\mu=$').replace("_n_",r' | $\nu=$'))
                ax2.plot(data.x,data.q, label=filename.replace("_"+file_filter,"").replace("m_",r'$\mu=$').replace("_n_",r' | $\nu=$'))
                ax1.set_xlabel(r'$\vartheta$')
                ax1.set_title(r'$\frac{u(\vartheta)}{\tan(\vartheta/2)}$')
                ax1.set_ylabel(r'$\frac{u(\vartheta)}{\tan(\vartheta/2)}$')
                ax2.set_xlabel(r'$\vartheta$')
                ax2.set_ylabel(r' $\partial_\vartheta u(\vartheta)$')
            else:
                ax1.plot(data.x,data.u, label=filename.replace("_"+file_filter,""))
                ax2.plot(data.x,data.q, label=filename.replace("_"+file_filter,""))
                ax1.set_xlabel('x')
                ax1.set_title('$u(x)$')
                ax2.set_xlabel('x')
                ax2.set_title('grad $u(x)$')
        else:
            continue

    if file_filter == "spin_fb.txt":
        plt.legend(title="", fancybox=True)
    else:
        plt.legend(title="# of iterations", fancybox=True)

    plt.show()


if __name__== "__main__":
    # first argument: which problem should be displayed
    main(sys.argv[1:])
