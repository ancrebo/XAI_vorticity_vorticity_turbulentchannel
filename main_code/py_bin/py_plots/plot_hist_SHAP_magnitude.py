# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------------------------------------------------
plot_hist_y.py
-------------------------------------------------------------------------------------------------------------------------
Created on Thu Mar 28 12:37:35 2024

@author: Andres Cremades Botella

File to plot the shap pdf in the domain. The file contains the following functions:
    - Functions:
        - plot_hist_y : function to plot the pdf of the velocities in the structures
"""
# -----------------------------------------------------------------------------------------------------------------------
# Import packages for all the functions
# -----------------------------------------------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import os

# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
# Define functions
# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

def plot_hist_SHAP_vol(data_in={"grid_shapvol":[],"grid_vol":[],"grid_shap":[],"plot_folder":"plots","xlabel":"$x^+$",
                                "ylabel":"$y^+$","zlabel":"z^+","fontsize":18,"figsize_x":10,"figsize_y":8,
                                "colormap":"viridis","colornum":2,"fig_name":"struc3d","dpi":60,"shap_max":1,
                                "shap_min":0,"vol_max":1,"vol_min":0,"cmap_flag":True,"padtext":[0,0,0],
                                "lev_min":1e-2,"nlev":None,"linewidth":1}):
    """""
    .....................................................................................................................
    # plot_histuvw_y: Function to plot the components of a vector in the field as a function of the wall-distance
    .....................................................................................................................
    Parameters
    ----------
    data_in : dict, optional
        Data required for generating the plot. 
        The default is {"grid_shapvol":[],"grid_vol":[],"grid_shap":[],"plot_folder":"plots","xlabel":"$x^+$",
                        "ylabel":"$y^+$","zlabel":"z^+","fontsize":18,"figsize_x":10,"figsize_y":8,
                        "colormap":"viridis","colornum":2,"fig_name":"struc3d","dpi":60,"shap_max":1,
                        "shap_min":0,"vol_max":1,"vol_min":0,"padtext":[0,0,0],
                        "lev_min":1e-2,"lev_delta":None,"linewidth":1}.
        Data:
            - grid_shapvol : data of the histogram
            - grid_vol     : volume values of the histogram
            - grid_shap    : shap values of the histogram
            - plot_folder  : folder to store the figures
            - xlabel       : label of the x axis
            - ylabel       : label of the y axis
            - zlabel       : label of the z axis
            - fontsize     : font size used for the figure
            - figsize_x    : size of the figure in x
            - figsize_y    : size of the figure in y
            - colormap     : colormap used for the figure
            - colornum     : number of colors of the colormap, two curves are used. The number of levels of the 
                             colormap needs to be higher than 2 
            - fig_name     : name of the saved figure
            - dpi          : dots per inch of the saved figure
            - shap_max     : maximum value of the SHAP, if not input, default
            - shap_min     : minimum value of the SHAP, if not input, default
            - vol_max      : maximum value of the volume, if not input, default
            - vol_min      : minimum value of the volume, if not input, default
            - cmap_flag    : flag for using colormaps
            - padtext      : padding used for the labels text
            - lev_min      : minimum level of the histogram
            - nlev         : number of levels
            - linewidth    : line width of the pdf

    Returns
    -------
    None.

    """
    
    from matplotlib import ticker
    from py_bin.py_class.plot_format import plot_format
    
    # -------------------------------------------------------------------------------------------------------------------
    # Read the parameters of the plot
    # -------------------------------------------------------------------------------------------------------------------
    grid_shapvol = np.array(data_in["grid_shapvol"],dtype="float")
    grid_vol     = np.array(data_in["grid_vol"],dtype="float")
    grid_shap    = np.array(data_in["grid_shap"],dtype="float")
    plot_folder  = str(data_in["plot_folder"]) 
    xlabel       = str(data_in["xlabel"]) 
    ylabel       = str(data_in["ylabel"])
    zlabel       = str(data_in["zlabel"])
    fontsize     = int(data_in["fontsize"])
    figsize_x    = int(data_in["figsize_x"])
    figsize_y    = int(data_in["figsize_y"])
    colormap     = data_in["colormap"]
    colornum     = int(data_in["colornum"])
    fig_name     = str(data_in["fig_name"])
    dpi          = float(data_in["dpi"])
    cmapflag     = False
    padtext      = np.array(data_in["padtext"],dtype="int")
    padtext_x    = padtext[0]
    padtext_y    = padtext[1]
    padtext_z    = padtext[2]
    lev_min      = float(data_in["lev_min"])
    nlev         = data_in["nlev"]
    linewidth    = int(data_in["linewidth"])
    if "shap_max" in data_in:
        shap_max = float(data_in["shap_max"])
    else:
        shap_max = np.max(shap)
    if "shap_min" in data_in:
        shap_min = float(data_in["shap_min"])
    else:
        shap_min = np.min(shap)
    if "vol_max" in data_in:
        vol_max = float(data_in["vol_max"])
    else:
        vol_max = np.max(vol)
    if "vol_min" in data_in:
        vol_min = float(data_in["vol_min"])
    else:
        vol_min = np.min(vol)
    
    if isinstance(colormap, str):
        colormap = str(colormap)
    else:
        colormap = np.array(colormap,dtype="float")
    
      
    # ----------------------------------------------------------------------------------------------------------------
    # Calculate the number of levels
    # ----------------------------------------------------------------------------------------------------------------
    if nlev is None:
        explev0   = np.log10(lev_min)
        explevs   = np.linspace(explev0,0,5)
        levels    = 10**explevs 
        locator   = ticker.LogLocator(numticks=5)
    else:
        explev0   = np.log10(lev_min)
        explevs   = np.linspace(explev0,0,nlev+1)
        levels    = 10**explevs 
        locator   = ticker.LogLocator(numticks=nlev+1)
            
    
    # -------------------------------------------------------------------------------------------------------------------
    # Create the plot for shap and volume
    # -------------------------------------------------------------------------------------------------------------------
    titlefig   = "$\phi - V^+$"
    data_plot  = {"xlabel":xlabel,"ylabel":ylabel,"zlabel":[],"fontsize":fontsize,"figsize_x":figsize_x,
                  "figsize_y":figsize_y,"xscale":"linear","yscale":"linear","zscale":"linear","colormap":colormap,
                  "colornum":colornum,"legend":True,"fig_name":fig_name,"dpi":dpi,"plot_folder":plot_folder,
                  "xmin":vol_min,"xmax":vol_max,"ymin":shap_min,"ymax":shap_max,"zmin":None,"zmax":None}
    plot_pred = plot_format(data_in=data_plot)
    plot_pred.create_figure()
    plot_pred.add_hist2d_y(data_in={"xx":grid_vol,"yy":grid_shap,"xxyy":grid_shapvol,"levels":levels,"colormap":None,
                                    "locator":locator,"alp":0.65,"linewidth":linewidth})
    labels = ['{0:.1e}'.format(levels[ii]) for ii in np.arange(len(levels))]
    plot_pred.plot_layout_pcolor(data_in={"title":titlefig,"colorbar":True,"b_text":None,
                                          "colorticks":levels,"colorlabels":labels,"equal":False,"xticks":None,
                                          "yticks":None,"xticklabels":None,"yticklabels":None})
    try:
        os.mkdir(plot_folder)
    except:
        print("Existing folder...",flush=True)
    # plot_pred.plot_save_png()
    # plot_pred.plot_save_pdf()
    # plot_pred.close()
   
    
def plot_hist_SHAP_vol_type(data_in={"grid_shapvol":[],"grid_vol":[],"grid_shap":[],"plot_folder":"plots",
                                     "xlabel":"$x^+$","ylabel":"$y^+$","zlabel":"z^+","fontsize":18,"figsize_x":10,
                                     "figsize_y":8,"colormap":"viridis","colornum":2,"fig_name":"struc3d","dpi":60,
                                     "shap_max":1,"shap_min":0,"vol_max":1,"vol_min":0,"cmap_flag":True,
                                     "padtext":[0,0,0],"lev_min":1e-2,"nlev":None,"linewidth":1,
                                     "labels_pdf":["0","1"],"colors_pdf":["#FF0000","#00FF00"],"typescale":"linear"}):
    """""
    .....................................................................................................................
    # plot_histuvw_y: Function to plot the components of a vector in the field as a function of the wall-distance
    .....................................................................................................................
    Parameters
    ----------
    data_in : dict, optional
        Data required for generating the plot. 
        The default is {"grid_shapvol":[],"grid_vol":[],"grid_shap":[],"plot_folder":"plots","xlabel":"$x^+$",
                        "ylabel":"$y^+$","zlabel":"z^+","fontsize":18,"figsize_x":10,"figsize_y":8,
                        "colormap":"viridis","colornum":2,"fig_name":"struc3d","dpi":60,"shap_max":1,
                        "shap_min":0,"vol_max":1,"vol_min":0,"padtext":[0,0,0],
                        "lev_min":1e-2,"lev_delta":None,"linewidth":1,"labels_pdf":["0","1"],
                        "colors_pdf":["#FF0000","#00FF00"]}.
        Data:
            - grid_shapvol : data of the histogram
            - grid_vol     : volume values of the histogram
            - grid_shap    : shap values of the histogram
            - plot_folder  : folder to store the figures
            - xlabel       : label of the x axis
            - ylabel       : label of the y axis
            - zlabel       : label of the z axis
            - fontsize     : font size used for the figure
            - figsize_x    : size of the figure in x
            - figsize_y    : size of the figure in y
            - colormap     : colormap used for the figure
            - colornum     : number of colors of the colormap, two curves are used. The number of levels of the 
                             colormap needs to be higher than 2 
            - fig_name     : name of the saved figure
            - dpi          : dots per inch of the saved figure
            - shap_max     : maximum value of the SHAP, if not input, default
            - shap_min     : minimum value of the SHAP, if not input, default
            - vol_max      : maximum value of the volume, if not input, default
            - vol_min      : minimum value of the volume, if not input, default
            - cmap_flag    : flag for using colormaps
            - padtext      : padding used for the labels text
            - lev_min      : minimum level of the histogram
            - nlev         : number of levels
            - linewidth    : line width of the pdf
            - labels_pdf   : labels to add to the colorbar
            - colors_pdf   : colors to use in the pdf

    Returns
    -------
    None.

    """
    
    from matplotlib import ticker
    from py_bin.py_class.plot_format import plot_format
    
    # -------------------------------------------------------------------------------------------------------------------
    # Read the parameters of the plot
    # -------------------------------------------------------------------------------------------------------------------
    grid_shapvol = data_in["grid_shapvol"]
    grid_vol     = data_in["grid_vol"]
    grid_shap    = data_in["grid_shap"]
    plot_folder  = str(data_in["plot_folder"]) 
    xlabel       = str(data_in["xlabel"]) 
    ylabel       = str(data_in["ylabel"])
    zlabel       = str(data_in["zlabel"])
    fontsize     = int(data_in["fontsize"])
    figsize_x    = int(data_in["figsize_x"])
    figsize_y    = int(data_in["figsize_y"])
    colormap     = data_in["colormap"]
    colornum     = int(data_in["colornum"])
    fig_name     = str(data_in["fig_name"])
    dpi          = float(data_in["dpi"])
    if "typescale" in data_in:
        typescale = str(data_in["typescale"])
    else:
        typescale = "linear"
    cmapflag     = False
    padtext      = np.array(data_in["padtext"],dtype="int")
    padtext_x    = padtext[0]
    padtext_y    = padtext[1]
    padtext_z    = padtext[2]
    lev_min      = float(data_in["lev_min"])
    nlev         = data_in["nlev"]
    linewidth    = int(data_in["linewidth"])
    labels_pdf   = data_in["labels_pdf"]
    colors_pdf   = data_in["colors_pdf"]
    if "shap_max" in data_in:
        shap_max = float(data_in["shap_max"])
    else:
        shap_max = np.max(shap)
    if "shap_min" in data_in:
        shap_min = float(data_in["shap_min"])
    else:
        shap_min = np.min(shap)
    if "vol_max" in data_in:
        vol_max = float(data_in["vol_max"])
    else:
        vol_max = np.max(vol)
    if "vol_min" in data_in:
        vol_min = float(data_in["vol_min"])
    else:
        vol_min = np.min(vol)
    
    if isinstance(colormap, str):
        colormap = str(colormap)
    else:
        colormap = np.array(colormap,dtype="float")
       
    # ----------------------------------------------------------------------------------------------------------------
    # Separe the structures
    # ----------------------------------------------------------------------------------------------------------------
    grid_shapvolQ2 = np.array(grid_shapvol[0],dtype="float")
    grid_volQ2     = np.array(grid_vol[0],dtype="float")
    grid_shapQ2    = np.array(grid_shap[0],dtype="float")
    labels_pdfQ2   = str(labels_pdf[0])
    colors_pdfQ2   = str(colors_pdf[0])
    grid_shapvolQ4 = np.array(grid_shapvol[1],dtype="float")
    grid_volQ4     = np.array(grid_vol[1],dtype="float")
    grid_shapQ4    = np.array(grid_shap[1],dtype="float")
    labels_pdfQ4   = str(labels_pdf[1])
    colors_pdfQ4   = str(colors_pdf[1])
    grid_shapvol_streakl = np.array(grid_shapvol[2],dtype="float")
    grid_vol_streakl     = np.array(grid_vol[2],dtype="float")
    grid_shap_streakl    = np.array(grid_shap[2],dtype="float")
    labels_pdf_streakl   = str(labels_pdf[2])
    colors_pdf_streakl   = str(colors_pdf[2])
    grid_shapvol_streakh = np.array(grid_shapvol[3],dtype="float")
    grid_vol_streakh     = np.array(grid_vol[3],dtype="float")
    grid_shap_streakh    = np.array(grid_shap[3],dtype="float")
    labels_pdf_streakh   = str(labels_pdf[3])
    colors_pdf_streakh   = str(colors_pdf[3])
    grid_shapvol_chong   = np.array(grid_shapvol[4],dtype="float")
    grid_vol_chong       = np.array(grid_vol[4],dtype="float")
    grid_shap_chong      = np.array(grid_shap[4],dtype="float")
    labels_pdf_chong     = str(labels_pdf[4])
    colors_pdf_chong     = str(colors_pdf[4])
    
    # ----------------------------------------------------------------------------------------------------------------
    # Calculate the number of levels
    # ----------------------------------------------------------------------------------------------------------------
    lev_max_Q2      = np.max(grid_shapvolQ2)
    exp_max_Q2      = np.log10(lev_max_Q2)
    lev_max_Q4      = np.max(grid_shapvolQ4)
    exp_max_Q4      = np.log10(lev_max_Q4)
    lev_max_streakl = np.max(grid_shapvol_streakl)
    exp_max_streakl = np.log10(lev_max_streakl)
    lev_max_streakh = np.max(grid_shapvol_streakh)
    exp_max_streakh = np.log10(lev_max_streakh)
    lev_max_chong   = np.max(grid_shapvol_chong)
    exp_max_chong   = np.log10(lev_max_chong)
    if nlev is None:
        explev0         = np.log10(lev_min)
        explevs_Q2      = np.linspace(explev0,exp_max_Q2,5)
        levels_Q2       = 10**explevs_Q2
        explevs_Q4      = np.linspace(explev0,exp_max_Q4,5)
        levels_Q4       = 10**explevs_Q4
        explevs_streakl = np.linspace(explev0,exp_max_streakl,5)
        levels_streakl  = 10**explevs_streakl
        explevs_streakh = np.linspace(explev0,exp_max_streakh,5)
        levels_streakh  = 10**explevs_streakh
        explevs_chong   = np.linspace(explev0,exp_max_chong,5)
        levels_chong    = 10**explevs_chong
        locator         = ticker.LogLocator(numticks=5)
    else:
        explev0         = np.log10(lev_min)
        explevs_Q2      = np.linspace(explev0,exp_max_Q2,nlev+1)
        levels_Q2       = 10**explevs_Q2
        explevs_Q4      = np.linspace(explev0,exp_max_Q4,nlev+1)
        levels_Q4       = 10**explevs_Q4
        explevs_streakl = np.linspace(explev0,exp_max_streakl,nlev+1)
        levels_streakl  = 10**explevs_streakl
        explevs_streakh = np.linspace(explev0,exp_max_streakh,nlev+1)
        levels_streakh  = 10**explevs_streakh
        explevs_chong   = np.linspace(explev0,exp_max_chong,nlev+1)
        levels_chong    = 10**explevs_chong
        locator         = ticker.LogLocator(numticks=nlev+1)
            
    
    # -------------------------------------------------------------------------------------------------------------------
    # Create the plot for shap and volume
    # -------------------------------------------------------------------------------------------------------------------
    titlefig   = " "
    data_plot  = {"xlabel":xlabel,"ylabel":ylabel,"zlabel":[],"fontsize":fontsize,"figsize_x":figsize_x,
                  "figsize_y":figsize_y,"xscale":typescale,"yscale":"linear","zscale":"linear","colormap":colormap,
                  "colornum":colornum,"legend":True,"fig_name":fig_name,"dpi":dpi,"plot_folder":plot_folder,
                  "xmin":vol_min,"xmax":vol_max,"ymin":shap_min,"ymax":shap_max,"zmin":None,"zmax":None}
    plot_pred = plot_format(data_in=data_plot)
    plot_pred.create_figure()
    plot_pred.add_hist2d_y(data_in={"xx":grid_volQ2,"yy":grid_shapQ2,"xxyy":grid_shapvolQ2,"levels":levels_Q2,
                                    "colormap":colors_pdfQ2,
                                    "locator":locator,"alp":0.5,"linewidth":linewidth,"cmapflag":False})
    
    plot_pred.add_hist2d_y(data_in={"xx":grid_volQ4,"yy":grid_shapQ4,"xxyy":grid_shapvolQ4,"levels":levels_Q4,
                                    "colormap":colors_pdfQ4,
                                    "locator":locator,"alp":0.5,"linewidth":linewidth,"cmapflag":False})
    plot_pred.add_hist2d_y(data_in={"xx":grid_vol_streakl,"yy":grid_shap_streakl,"xxyy":grid_shapvol_streakl,
                                    "levels":levels_streakl,
                                    "colormap":colors_pdf_streakl,
                                    "locator":locator,"alp":0.5,"linewidth":linewidth,"cmapflag":False})
    plot_pred.add_hist2d_y(data_in={"xx":grid_vol_streakh,"yy":grid_shap_streakh,"xxyy":grid_shapvol_streakh,
                                    "levels":levels_streakh,
                                    "colormap":colors_pdf_streakh,
                                    "locator":locator,"alp":0.5,"linewidth":linewidth,"cmapflag":False})
    plot_pred.add_hist2d_y(data_in={"xx":grid_vol_chong,"yy":grid_shap_chong,"xxyy":grid_shapvol_chong,
                                    "levels":levels_chong,
                                    "colormap":colors_pdf_chong,
                                    "locator":locator,"alp":0.5,"linewidth":linewidth,"cmapflag":False})
    labels = ['{0:.1e}'.format(levels_Q2[ii]) for ii in np.arange(len(levels_Q2))]
    plot_pred.plot_layout_pcolor(data_in={"title":titlefig,"colorbar":True,"b_text":None,
                                          "colorticks":levels_Q2,"colorlabels":labels,"equal":False,"xticks":None,
                                          "yticks":None,"xticklabels":None,"yticklabels":None,
                                          "legend_color":[colors_pdfQ2,colors_pdfQ4,colors_pdf_streakl,colors_pdf_streakh,colors_pdf_chong],
                                          "flag_legend_types":True,"legend_label_types":[labels_pdfQ2,labels_pdfQ4,
                                                                                         labels_pdf_streakl,
                                                                                         labels_pdf_streakh,labels_pdf_chong]})
    try:
        os.mkdir(plot_folder)
    except:
        print("Existing folder...",flush=True)
    plot_pred.plot_save_png()
    plot_pred.plot_save_pdf()
    plot_pred.close()
    