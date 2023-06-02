import sys
import os

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plt.rcParams.update({'font.size': 26})

def main(argv):
    # filepath = sys.argv[1]
    # if not os.path.isfile(filepath):
    #   print("File path {} does not exist. Exiting...".format(filepath))
    #   sys.exit()
    plt.style.use(['seaborn-deep'])
    fig, ax1 = plt.subplots(1, 1)
    fig = plt.gcf() 
    fig.set_size_inches(13,10)

    folder = "pltdata/"
    
    
    
    
    if not argv:
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
        if filename.endswith(file_filter):
            data = pd.read_csv(folder+filename, sep=" ", header=None)
            data.columns = ["x", "p", "u", "q", "index", "newl"]
            #ax1.plot(data.x,data.u*np.cos(0.5*data.x)/np.sin(0.5*data.x), label=filename.replace("_"+file_filter,""))
            ax1.plot(data.x,data.u, "o", label= filename.replace("_"+file_filter,"").replace("m_",r'$\mu=$').replace("_n_",r', $\nu=$'), alpha = 0.5)
            ax1.set_xlabel(r'$\vartheta$')
            ax1.set_title(r'osmotic velocities for different $\mu, \nu$')
            ax1.set_ylabel(r'${\tilde{\Omega}^u_\vartheta(\vartheta)}$')
            ax1.legend(frameon=True, loc='lower center', ncol=1)

        else:
            continue

    plt.gca().set_prop_cycle(None)

    file_filter = "spin_fb_exact.txt"
    for filename in os.listdir(folder):
        if filename.endswith(file_filter):
            data = pd.read_csv(folder+filename, sep=" ", header=None)
            data.columns = ["x", "p", "u", "q", "index", "newl"]
            #ax1.plot(data.x,data.u*np.cos(0.5*data.x)/np.sin(0.5*data.x), label=filename.replace("_"+file_filter,""))
            ax1.plot(data.x,data.u, "--", linewidth=2, label= filename.replace("_"+file_filter,"").replace("m_",r'$\mu=$').replace("_n_",r', $\nu=$') + " exact")
            ax1.set_xlabel(r'$\vartheta$')
            ax1.set_title(r'osmotic velocities for different $\mu, \nu$')
            ax1.set_ylabel(r'$\tilde{\Omega}^u_\vartheta(\vartheta)$')
            ax1.legend(frameon=True, loc='lower center', ncol=1)

        else:
            continue
    
    #plt.legend(title="", fancybox=True)
    plt.show()


if __name__== "__main__":
    # first argument: which problem should be displayed
    main(sys.argv[1:])
 
