#
#   plot_contours
#                               
#   Script to plot a converged flowfield from the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_contours.py casename"

# Import modules and functions
from routines import *

def main():

    # Construct full filenames to read the run data
    inname = 'input_' + sys.argv[-1] + '.txt'
    outname = 'out_final_' + sys.argv[-1] + '.bin'

    # Read the settings and the case from file
    av = read_settings(inname)
    g = read_case(outname)

    # When presenting results all values should be non-dimensionalised. Two
    # variables of interest might be:
    #    1. Static pressure coefficient, (p - p_ref) / (pstag_ref - p_ref)
    #    2. Mach number, v / (ga * rgas * t)**0.5

    # First complete the "calc_secondary" function within "routines.py" to
    # calculate static pressure and Mach number, and any others you want!
    g = calc_secondary(av,g)    

    # Use the "cut_i", "mass_av" AND "area_av" functions to calculate the
    # reference pressures at the inlet plane and therefore the static pressure
    # coefficient
    # INSERT
    #########################
    i_inlet = 0
    cut = cut_i(g,i_inlet)
    p0_ref, m= mass_av(cut, 'pstag')
    p_ref, l= area_av(cut, 'p')

    g['cp'] = (g['p'] - p_ref)/(p0_ref - p_ref)
    pstag = av['pstag']
    p = av['p'][0]
    print(pstag, p)
    g['cp0'] = (g['pstag'] - pstag)/(pstag - p)
    
    #########################
    # Specify the parameters to plot
    fieldnames = ['cp', 'mach', 'cp0']; 
    colnames = ['Static pressure coefficient','Mach number', 'Stagnation pressure coefficient']

    # Plot the calculated non-dimensional parameters to show the flow solution
    for n,name in enumerate(fieldnames):

        # Open figure window
        fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();
    
        # Set aspect ratio as equal and remove axes labels
        ax.set_aspect('equal',adjustable='box'); ax.axis('off')
 
        # Plot filled contour levels
        hc = ax.pcolormesh(g['x'],g['y'],g[name],shading='gouraud')

        # Add colorbar with variable name
        colorbar(hc,colnames[n])

        # Add Mach = 1 contours
        if name == 'mach':
            ax.contour(g['x'],g['y'],g['mach'],[0.7],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['mach'],[0.75],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['mach'],[0.8],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['mach'],[0.85],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['mach'],[0.9],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['mach'],[0.95],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['mach'],[1.00],colors='w',
                linewidths=0.5)
 #           ax.contour(g['x'],g['y'],g['mach'],[0.65],colors='w',
 #               linewidths=0.5)
        elif name == 'cp':
            ax.contour(g['x'],g['y'],g['cp'],[-0.5],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['cp'],[-0.4],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['cp'],[-0.3],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['cp'],[-0.2],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['cp'],[-0.1],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['cp'],[0.0],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['cp'],[0.1],colors='w',
                linewidths=0.5)
        elif name == 'cp0':
            ax.contour(g['x'],g['y'],g['cp0'],[-0.2],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['cp0'],[-0.15],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['cp0'],[-0.125],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['cp0'],[-0.10],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['cp0'],[-0.075],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['cp0'],[-0.05],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['cp0'],[0.05],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['cp0'],[0.0],colors='w',
                linewidths=0.5)
            ax.contour(g['x'],g['y'],g['cp0'],[0.025],colors='w',
                linewidths=0.5)
        

        # Draw the walls of the block
        plot_wall(ax,g)

    # Show all the plots
    plt.show()

    
main()


