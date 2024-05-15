"""
Main module for calculating and plotting energy as a function of x for chosen material
"""
import matplotlib.pyplot as plt


#############################
########    CONFIG    #######
#############################

PLOT_TITLE = r'$ InAs_{x}Sb_{1-x} $'
X_LABEL = r'x value of $ As_x $ composition'
PARAMETERS_BANDGAP = [
    # Format should be a tuple containing (B_BC,B_AC,bowing parameter)
    # where material format is: A_xB_(1-x)C
    # One tuple corresponds to one symmetry point
    (0.235,0.417,0.67) #Point Gamma
    ,(0.63,1.433,0.6) #Point X
    ,(0.93,1.133,0.6) #Point L
]
PLOT_BANDS = True
PARAMETERS_BANDS = {
    # Format should be a tuple containing (a,b) where a is the value
    # of the first materials parameter (the one with 1-x as index) and
    # b is the value of the second materials parameter with the exception
    # of lattice_a_0 where only 1 value should be present
    "VBO": (0 , -0.59),
    "b":(-2.0 , -1.8),
    "a_c":(-6.94 , -5.08),
    "a_v":(-0.36 , -1.00),
    "C_11":(684.7 , 832.9),
    "C_12":(373.5 , 452.6),
    "lattice_a": (6.4794, 6.0583),
    "lattice_a_0": 6.0583,
}
TEMPERATURE_DEPENDENCE = True
PARAMETERS_TEMPERATURE = {
    # Temperatures is a list of temperatures for which the calculations will be done
    # alpha and beta are in the same format as params in previous dictionary
    "Temperatures" : [0, 10, 300],
    "alpha" : (0.32*10**-3,0.276*10**-3),
    "beta" : (170,93),
}
QUANTUM_WELL_PARAMS = {
    "Well_widths":(100,200,300),
    "Compositions":(0.25,0.5,0.75)
}
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
                    show_image:bool = True,
                    image_name:str = "energy_plot.png"):
    """
    This function generates plot for all the given energies
    """
    figure, axes = plt.subplots(figsize=(20,16))
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
        figure.savefig(image_name)

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

def interpolate(values_to_interpolate:tuple,
                x_list:list):
    temp = []
    for x in x_list:
        temp.append(values_to_interpolate[0]*(1-x) + values_to_interpolate[1]*x)
    return temp

def include_temperature(energy_gap:list,
                        x_list: list,
                        temperature: int = 0):
    alpha = interpolate(PARAMETERS_TEMPERATURE['alpha'],x_list)
    beta = interpolate(PARAMETERS_TEMPERATURE['beta'],x_list)
    new_energy = [energy_gap[i] - (alpha[i]*(temperature**2) / (temperature + beta[i])) for i in range(len(energy_gap))]
    return new_energy

def calculate_not_strained_bands(energy_gap:list,
                                x_list:list):
    valence_band = interpolate(PARAMETERS_BANDS['VBO'],x_list)
    conduction_band = [valence_band[i] + energy_gap[i] for i in range(len(x_list))]
    return (valence_band,conduction_band)

def calculate_strained_bands(valence_band:list,
                            conduction_band:list,
                            x_list:list):
    lattice_a_0 = PARAMETERS_BANDS['lattice_a_0']
    interpolated_lattice = interpolate(PARAMETERS_BANDS['lattice_a'],x_list)
    c11_interpolated = interpolate(PARAMETERS_BANDS['C_11'],x_list)
    c12_interpolated = interpolate(PARAMETERS_BANDS['C_12'],x_list)
    b_interpolated = interpolate(PARAMETERS_BANDS['b'],x_list)
    a_c_interpolated = interpolate(PARAMETERS_BANDS['a_c'],x_list)
    a_v_interpolated = interpolate(PARAMETERS_BANDS['a_v'],x_list)
    eps_x_list = [(lattice_a_0 - interpolated_lattice[i])/interpolated_lattice[i] for i in range(len(x_list))]
    eps_z_list = [eps_x_list[i]*-2*(c12_interpolated[i] / c11_interpolated[i]) for i in range(len(x_list))]
    dE_s = [-b_interpolated[i] * (eps_z_list[i]-eps_x_list[i]) for i in range(len(x_list))]
    dE_hc = [a_c_interpolated[i]*(2*eps_x_list[i]+eps_z_list[i]) for i in range(len(x_list))]
    dE_hv = [a_v_interpolated[i]*(2*eps_x_list[i]+eps_z_list[i]) for i in range(len(x_list))]
    strained_conduction_band = [conduction_band[i] + dE_hc[i] for i in range(len(x_list))]
    heavy_holes = [valence_band[i] + dE_hv[i] + dE_s[i] for i in range(len(x_list))]
    light_holes = [valence_band[i] + dE_hv[i] - dE_s[i] for i in range(len(x_list))]
    return (strained_conduction_band,heavy_holes,light_holes)

def create_quantum_well(energy_list: list,
                            valence_band:list,
                            conduction_band:list,
                            heavy_holes: list,
                            light_holes: list,
                            conduction_tension: list,
                            well_width: float,
                            composition: int,
                            temperature: int):
    """
    This function creates a plot of quantum well with the given material
    """
    composition = int(composition * 100)
    position = [i for i in range(0,1001)]
    valence_hh = [0 for i in position]
    valence_lh = [0 for i in position]
    conduction = [0 + energy_list[0] for i in position]
    conduction_no_tension = [0 + energy_list[0] for i in position]
    valence_no_tension = [0 for i in position]
    for i,x in enumerate(position):
        if x in range(int(500-well_width/2.0),int(501+well_width/2.0)):
            valence_hh[i] = heavy_holes[composition] - PARAMETERS_BANDS["VBO"][1]
            valence_lh[i] = light_holes[composition] - PARAMETERS_BANDS["VBO"][1]
            conduction[i] = conduction_tension[composition] - PARAMETERS_BANDS["VBO"][1]
            conduction_no_tension[i] = conduction_band[composition] - PARAMETERS_BANDS["VBO"][1]
            valence_no_tension[i] = valence_band[composition] - PARAMETERS_BANDS["VBO"][1]
            if x==500:
                print(f"temp:{temperature},comp:{composition},width:{well_width}")
                print(f"Naprężone: {conduction[i]-valence_hh[i]}")
                print(f'nienaprezone: {conduction_no_tension[i] - valence_no_tension[i]}')
    
    fig, axes = plt.subplots(figsize=(14,10),nrows=1,ncols=2)
    axes[0].plot(position,valence_hh,'b',label='$ E_{V-HEAVY-HOLES}$')
    axes[0].plot(position,valence_lh,'--c',label='$ E_{V-LIGHT-HOLES}$')
    axes[0].plot(position,conduction,'m',label='$ E_{C-WITH-TENSION} $')
    axes[0].legend(fontsize=12,loc='center left')
    axes[0].set_title("With tension",fontsize=20)
    axes[0].set_ylabel('Energy [eV]',fontsize=20)
    axes[0].set_xlabel('Position [Å]',fontsize=20)
    axes[0].grid()
    #fig.savefig(f"Quantum_well_with_tension_{temperature}K_{well_width}nm_composition_{composition}.png")
    #fig, axes = plt.subplots(figsize=(12,8))
    axes[1].plot(position,valence_no_tension,'g',label='$ E_V $')
    axes[1].plot(position,conduction_no_tension,'r',label='$ E_C $')
    axes[1].legend(fontsize=12,loc='center right')
    axes[1].set_title("Without tension",fontsize=20)
    #axes[1].set_ylabel('Energy [eV]',fontsize=20)
    axes[1].set_xlabel('Position [Å]',fontsize=20)
    axes[1].grid()
    fig.suptitle('$ In_{'+f"{composition/100}"+'}As_{'+f"{(100-composition)/100}"+'}Sb $ - InAs base for ' + f"T = {temperature}K and width = {int(well_width/10)} nm",fontsize=20)
    fig.savefig(f"Quantum_well_{temperature}K_{int(well_width/10)}nm_composition_{composition}.png")


def main():
    """
    Main of the whole module, invokes all the functions needed for
    the calculation and plot of the material
    """
    energy_lists = []
    for params in PARAMETERS_BANDGAP:
        a_val,b_val,c_val = get_parameters(params)
        x_list,energy_list = calculate_energy(value_a=a_val,value_b=b_val,value_c=c_val)
        energy_lists.append(energy_list)
    create_all_plots(x_list,energy_lists,LABELS,PLOT_TITLE,SAVE_PLOT_AS_PNG,SHOW_IMAGE)
    if PLOT_BANDS and not TEMPERATURE_DEPENDENCE:
        energy_gap = energy_lists[0]
        valence_band,conduction_band = calculate_not_strained_bands(energy_gap,x_list)
        strained_conduction_band, heavy_holes,light_holes = calculate_strained_bands(valence_band,
                                                                                    conduction_band,
                                                                                    x_list)
        band_list = [valence_band,conduction_band,strained_conduction_band,heavy_holes,light_holes]
        labels = ['$ E_V $', '$ E_C $', '$ E_{C-with-strain} $', '$ E_{HH} $', '$ E_{LH} $']
        create_all_plots(x_list,band_list,labels,PLOT_TITLE,
                        SAVE_PLOT_AS_PNG,SHOW_IMAGE,'Energy_bands_with_strain.png')
    if PLOT_BANDS and TEMPERATURE_DEPENDENCE:
        all_bands = []
        all_labels = []
        for i,temperature in enumerate(PARAMETERS_TEMPERATURE['Temperatures']):
            energy_gap = include_temperature(energy_lists[0],x_list,temperature)
            valence_band,conduction_band = calculate_not_strained_bands(energy_gap,x_list)
            strained_conduction_band, heavy_holes,light_holes = calculate_strained_bands(valence_band,
                                                                                        conduction_band,
                                                                                        x_list)
            # We generate valence, hh and lh band only for 1st T since they're not T dependent
            if not i:
                band_list = [valence_band,conduction_band,
                            strained_conduction_band,heavy_holes,light_holes]
                labels = ['$ E_V $', f'$ E_C ({temperature}K) $',
                         '$ E_{C-with-strain} '+ f'({temperature}K) $', '$ E_{HH} $', '$ E_{LH} $']
                single_temp_bands = band_list
                single_temp_labels = labels
            else:
                band_list = [conduction_band,strained_conduction_band]
                labels = [f'$ E_C ({temperature}K) $','$ E_{C-with-strain} '+ f'({temperature}K) $']
                single_temp_bands = [all_bands[0],all_bands[3],all_bands[4]]
                single_temp_labels = [all_labels[0],all_labels[3],all_labels[4]]
                single_temp_bands.extend(band_list)
                single_temp_labels.extend(labels)
            all_bands.extend(band_list)
            all_labels.extend(labels)
            create_all_plots(x_list,single_temp_bands,single_temp_labels,PLOT_TITLE,
                        SAVE_PLOT_AS_PNG,SHOW_IMAGE,f'Energy_bands_with_strain_for_{temperature}K.png')
            
            compositions = QUANTUM_WELL_PARAMS['Compositions']
            widths = QUANTUM_WELL_PARAMS['Well_widths']
            for width in widths:
                for composition in compositions:
                    create_quantum_well(energy_list = energy_list,
                                    valence_band=valence_band,
                                    conduction_band=conduction_band,
                                    heavy_holes = heavy_holes,
                                    light_holes = light_holes,
                                    conduction_tension = strained_conduction_band,
                                    well_width = width,
                                    composition = composition,
                                    temperature = temperature)
        create_all_plots(x_list,all_bands,all_labels,PLOT_TITLE,
                        SAVE_PLOT_AS_PNG,SHOW_IMAGE,'Energy_bands_with_strain_and_temperature.png')
    
    return 0

if __name__ == '__main__':
    main()
