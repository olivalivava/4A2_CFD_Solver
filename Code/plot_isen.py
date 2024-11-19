#
#   plot_conv                      
#                               
#   Script to plot mass and tstag conservation history of a run executed using the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_conv.py casename"

# Import modules and functions
from routines import *

def main():

    # Construct full filenames to read the run data
    inname = 'input_' + sys.argv[-1] + '.txt'
    outname = 'out_final_' + sys.argv[-1] + '.bin'

    # Read the settings and the case from file
    av = read_settings(inname)
    g = read_case(outname)

    g = calc_secondary(av,g)

    i_inlet = 0
    
    cut_in = cut_i(g,i_inlet)

    ni = av['ni']
    nj = av['nj']
    j_mid = round(nj/2)
    i=0
    j=j_mid
    
    t0_in = cut_in['tstag']
    mass_in = cut_in['rovx'][j]*cut_in['lx'][j] + cut_in['rovy'][j]*cut_in['ly'][j]
    rov_ratio = []
    mass = []
    t0_ratio = []

    
    
    j_index = np.linspace(0, ni-2, ni-1).tolist()
    
    for i in range(ni-1):
        #j_max = np.argmax(abs(g['tstag'][i][j]-t0_in))
        t0_ratio.append(g['tstag'][i][j_mid]/t0_in)
        t0_1d = [sublist[0] for sublist in t0_ratio]

        mass = g['rovx'][i][j]*g['lx_i'][i][j] + g['rovy'][i][j]*g['ly_i'][i][j]
        rov_ratio.append(mass/mass_in)
        rov_1d = rov_ratio
        #rov_1d = [sublist[0] for sublist in rov_ratio]
        #print(rov_1d)
    
    # Open figure window for all residual data
    plt.figure(figsize=[9.6,7.2]); ax = plt.axes(); cols = gen_cols();
    ax.set_xlabel('Cell Index'); ax.set_ylabel('Stagnation Temperature Ratio');

    # Set y-axis as a log scale and turn on the gridlines
    ax.set_yscale('linear'); ax.tick_params(direction='in',which='both');     
    ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
    ax.grid(linestyle='-',color=[0.8,0.8,0.8],linewidth=0.5,which='minor',
        axis='y')
    
    #ax.plot(j_index, t0_ratio_1d)
    ax.plot(j_index, t0_1d, label="Stagnation Temperature Ratio", linestyle='-', marker='o')

    # Open figure window for all residual data
    plt.figure(figsize=[9.6,7.2]); ax = plt.axes(); cols = gen_cols();
    ax.set_xlabel('Cell Index'); ax.set_ylabel('Mass Flow Ratio');

    # Set y-axis as a log scale and turn on the gridlines
    ax.set_yscale('linear'); ax.tick_params(direction='in',which='both');     
    ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
    ax.grid(linestyle='-',color=[0.8,0.8,0.8],linewidth=0.5,which='minor',
        axis='y')
    
    ax.plot(j_index, rov_1d, label="Mass Flow Ratio", linestyle='-', marker='o')

    ax.legend()

    plt.show()
main()


