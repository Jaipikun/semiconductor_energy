"""
Main module for calculating and plotting energy as a function of x for chosen material
"""
import matplotlib.pyplot as plt


#############################
########    CONFIG    #######
#############################

PLOT_TITLE = r'$ InAs_{x}Sb_{1-x} $'
X_LABEL = r'x value of $ As_x $ composition'
PARAMETERS = [
    # Format should be a tuple containing (B_BC,B_AC,bowing parameter)
    # where material format is: A_xB_(1-x)C
    # One tuple corresponds to one symmetry point
    (0.235,0.417,0.67) #Point Gamma
    ,(0.63,1.433,0.6) #Point X
    ,(0.93,1.133,0.6) #Point L
]
LABELS = [
    r'$ \Gamma $',
    r'X',
    r'L'
]
FONT_SIZE = 20
SAVE_PLOT_AS_PNG = True
SHOW_IMAGE = False

#############################

def calculate_energy(value_a:float, value_b:float, value_c:float):
    """
    This function calculates the energy for given parameters
    y = a + bx + cx^2
    Args:
        value_a: parameter a
        value_b: parameter b
        value_c: parameter c
    Returns:
        Tuple containing list of xs and their corresponding energy values in (x,y) format
    """
    x_list = [i/100.0 for i in range(101)]
    energy_list = []
    for x_val in x_list:
        energy_value = value_a + value_b*x_val + value_c*x_val*x_val
        energy_list.append(energy_value)
    return (x_list,energy_list)

def create_all_plots(x_list,
                    energy_lists:list,
                    labels:list = [f'Default label #{num}' for num in range(20)],
                    material_name:str = 'Default material name',
                    save_to_png:bool = False,
                    show_image:bool = True):
    """
    This function generates plot for all the given energies
    """
    figure, axes = plt.subplots(figsize=(12,8))
    for i,energy_values in enumerate(energy_lists):
        axes.plot(x_list,energy_values,label=labels[i])
    axes.grid()
    axes.set_ylabel("Energy [eV]",fontsize=FONT_SIZE)
    axes.set_xlabel(X_LABEL,fontsize=FONT_SIZE)
    axes.set_title(material_name,fontsize=FONT_SIZE)
    axes.legend(fontsize=FONT_SIZE)
    axes.tick_params(axis='both',labelsize=FONT_SIZE-6)
    if show_image:
        figure.show()
        figure.waitforbuttonpress()
    if save_to_png:
        figure.savefig("energy_plot.png")

def get_parameters(parameters:tuple = (0,0,0)):
    """
    Simple function which gets parameters needed for energy calculation

    Args:
        parameters: a tuple which contains values from which parameters are obtained
    Returns:
        A tuple (a,b,c), where a,b,c are the needed parameters
    """
    a,b_ac,c = parameters
    b = b_ac - a - c
    return (a,b,c)

def main():
    """
    Main of the whole module, invokes all the functions needed for
    the calculation and plot of the material
    """
    energy_lists = []
    for params in PARAMETERS:
        a_val,b_val,c_val = get_parameters(params)
        x_list,energy_list = calculate_energy(value_a=a_val,value_b=b_val,value_c=c_val)
        energy_lists.append(energy_list)

    create_all_plots(x_list,energy_lists,LABELS,PLOT_TITLE,SAVE_PLOT_AS_PNG,SHOW_IMAGE)
    return 0

if __name__ == '__main__':
    main()
