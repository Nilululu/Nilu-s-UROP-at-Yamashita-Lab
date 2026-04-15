"""
Creates histograms and scatter plots using matplotlib
"""
import numpy as np
import matplotlib.pyplot as plt



def create_hist( data, title):
    """

    Parameters
    data : python list
    title : string

    Raises Error if the data or title are not of expected type

    Returns: None

    """
    if not isinstance(data, list):
        if not isinstance(data, np.ndarray):
            raise TypeError("Error: the data given is not in list or numpy array format")
    if not isinstance(title, str):
        raise TypeError("Error: title should be a string")
    fig, ax = plt.subplots(1)
    
    ax.hist(data , 100)
    ax.set_title(title)
    
    plt.tight_layout()
    plt.show()
    
    return
def create_scatter (x, y, z, title, x_label, y_label, z_label):
    """
    used for comparing introns relative position to the gene and their length

    Creates a plot and returns None

    """
    
    if not isinstance(x, list):
        if not isinstance(x, np.ndarray):
            raise TypeError("Error: the input data is not in list format")
   
    
    # if type(title) != str:
        
    #     raise TypeError("Error: the title is not a string")
    
    # if type(x_label) or type(y_label) or type(z_label) != str:
        
    #     raise TypeError("Error: the input labels are not all strings")
     
    fig, ax = plt.subplots(1)
    
    sc= ax.scatter(x, y, c=z, cmap="viridis")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    cbar = plt.colorbar(sc)
    cbar.set_label(z_label)

    plt.title(title)
   
    plt.tight_layout()
    plt.show()
        
    return 


