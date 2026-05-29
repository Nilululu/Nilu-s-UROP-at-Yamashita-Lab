"""
Creates histograms and scatter plots using matplotlib
"""
import numpy as np
import matplotlib.pyplot as plt
import random


def create_hist( data, title):
    """
    Parameters
    data : python list
    title : string
    Raises Error if the data or title are not a string or data is not a list/numpy array
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
    Parameters
    x, y, z : python lists or numoy arrays
    title and labels: string
    Returns: None
    """
    
    if not isinstance(x, list):
        if not isinstance(x, np.ndarray):
            raise TypeError("Error: the input data is not in list format")
     
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

def create_2d_scatter (x, y, title, x_label = None, y_label = None):
    """
    Parameters
    x, y: python lists or numoy arrays
    title: string
    Returns: None
    """
    
    fig, ax = plt.subplots(1)
    
    ax.scatter(x, y)
    ax.set_title(title)

    if x_label:
        ax.set_xlabel(x_label)
    if y_label:
        ax.set_ylabel(y_label)
    
    plt.tight_layout()
    plt.show()
    


def get_style (kingdom, _mapping = dict(), _used = set()):
    
    colors = ["red", "blue", "green", "purple", "orange", 
              "cyan", "magenta", "yellow", "brown", "pink"]
    
    shapes = ["o","s","^","D", "*", "h", "p", "x"]
     
#         "circle": "o", "square": "s", "triangle": "^", "diamond": "D"
#         "star": "*", "hexagon": "h", "pentagon": "p", "cross": "x"

    #returning existing info
    if kingdom in _mapping:
        return _mapping
    
    #looping until we get a color and shape that has not been used before
    while True:
        color = random.choice(colors)
        
        if color in _used:
            continue
        shape = random.choice(shapes)
        
        if shape in _used:
            continue
        
        _used.add(color)
        _used.add(shape)
        
        break 
        
    _mapping[kingdom] = {"color" : color, "marker" : shape}
    
    return _mapping
        
    
# def create_filtered_hist (genomes_dict, filter_name, filter_rank, question):
#     """

#     Parameters
#     ----------
#     genomes_dict : python dict
#         gene dic for genomes and gintrons dict of all genomes.
#     filter_name : str
#         animals, insects, etc.
#     filter_rank : str
#         kingdom, subkingdom, family, etc.
#     question: str
#         gintrons, max, mean, num_gintrons, percentile. 

#     Returns
#     -------
#     None.
    

#     """
    
#     for keys in genomes_dict:
#         if genomes_dict [keys]["taxId"][filter_rank] == filter_name:
            
#     #filteration 
    
    