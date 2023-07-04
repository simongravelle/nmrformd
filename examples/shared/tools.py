from matplotlib import pyplot as plt
# from matplotlib.ticker import AutoMinorLocator
# from mpl_toolkits.axes_grid1 import make_axes_locatable

def figure_improvement(ax1, mygray, font, fontsize, xlabel=None, ylabel=None, xlim=None, ylim=None):

    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    ax1.minorticks_on()
    ax1.tick_params('both', length=10, width=2, which='major', direction='in')
    ax1.tick_params('both', length=6, width=1.4, which='minor', direction='in')
    ax1.xaxis.set_ticks_position('both')
    ax1.yaxis.set_ticks_position('both')
    ax1.spines["top"].set_linewidth(2)
    ax1.spines["bottom"].set_linewidth(2)
    ax1.spines["left"].set_linewidth(2)
    ax1.spines["right"].set_linewidth(2)

    ax1.xaxis.label.set_color(mygray)
    ax1.yaxis.label.set_color(mygray)
    ax1.tick_params(axis='x', colors=mygray)
    ax1.tick_params(axis='y', colors=mygray)
    ax1.spines['left'].set_color(mygray)
    ax1.spines['top'].set_color(mygray)
    ax1.spines['bottom'].set_color(mygray)
    ax1.spines['right'].set_color(mygray)
    ax1.tick_params(axis='y', which='both', colors=mygray)
    ax1.tick_params(axis='x', which='both', colors=mygray)
    ax1.tick_params(axis='x', pad=10)

    ax1.legend(frameon=False, fontsize=fontsize, 
            loc='best', handletextpad=0.5,
            handlelength = 0.2, borderpad = 0.3, 
            labelspacing=0.3, labelcolor=mygray) 

    if xlabel is not None:
        ax1.set_xlabel(xlabel, fontdict=font, color=mygray)
    if ylabel is not None:
        ax1.set_ylabel(ylabel, fontdict=font, color=mygray)
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)

    #ax1.set_xticks([0, 5, 10, 15, 20])
    #ax1.set_yticks([0, 0.4, 0.8, 1.2, 1.6])
    #minor_locator_y = AutoMinorLocator(2)
    #ax1.yaxis.set_minor_locator(minor_locator_y)
    #minor_locator_x = AutoMinorLocator(2)
    #ax1.xaxis.set_minor_locator(minor_locator_x)

def save_figure(plt, fig, mode, name):
    fig.tight_layout()        
    if mode == 'light':
        file = "figures/" + name + "-light.png"
    else:
        file = "figures/" + name + "-dark.png"
    plt.savefig(file, bbox_inches = 'tight', pad_inches = 0.057, transparent=True)
    plt.show()    
