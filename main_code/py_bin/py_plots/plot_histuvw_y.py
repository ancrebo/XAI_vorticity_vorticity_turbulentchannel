# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------------------------------------------------
plot_histuvw_y.py
-------------------------------------------------------------------------------------------------------------------------
Created on Thu Mar 28 12:37:35 2024

@author: Andres Cremades Botella

File to plot the velocity pdf in a type of structure. The file contains the following functions:
    - Functions:
        - plot_histuvw_y : function to plot the pdf of the velocities in the structures
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

def plot_histuvw_y(data_in={"plot_folder":"plots","plot_fileu":"fileu","plot_filev":"filev","plot_filew":"filew",
                            "ylabel":"y","xlabelu":"u","xlabelv":"v","xlabelw":"w","fontsize":18,"figsize_x":10,
                            "figsize_y":8,"colormap":"viridis","colornum":2,"dpi":60,"uu_struc":[],"vv_struc":[],
                            "ww_struc":[],"yplus_struc":[],"yplusmesh":[],"bins":1,"lev_min":1e-2,"lev_delta":None,
                            "linewidth":2,"umin":None,"umax":None,"vmin":None,"vmax":None,"wmin":None,"wmax":None}):
    """""
    .....................................................................................................................
    # plot_histuvw_y: Function to plot the pdf of the velocities in the structures
    .....................................................................................................................
    Parameters
    ----------
    data_in : dict, optional
        Data required for generating the plot. 
        The default is {"plot_folder":"plots","plot_fileu":"fileu","plot_filev":"filev","plot_filew":"filew",
                        "xlabel":"y","ylabelu":"u","ylabelv":"v","ylabelw":"w","fontsize":18,"figsize_x":10,
                        "figsize_y":8,"colormap":"viridis","colornum":2,"dpi":60,"uu_struc":[],"vv_struc":[],
                        "ww_struc":[],"yplus_struc":[],"yplusmesh":[],"bins":1,"lev_min":1e-2,"lev_delta":None,
                        "linewidth":2,"umin":None,"umax":None,"vmin":None,"vmax":None,"wmin":None,"wmax":None}.
        Data:
            - plot_folder      : folder to store the plots
            - plot_fileu       : file to save the pdf of the streamwise velocity
            - plot_filev       : file to save the pdf of the wall-normal velocity
            - plot_filew       : file to save the pdf of the spanwise velocity
            - ylabel           : label of the y axis
            - xlabelu          : label of the x axis for the streamwise velocity
            - xlabelv          : label of the x axis for the wall-normal velocity
            - xlabelw          : label of the x axis for the spanwise velocity
            - fontsize         : font size used for the figure
            - figsize_x        : size of the figure in x
            - figsize_y        : size of the figure in y
            - colormap         : colormap used for the figure
            - colornum         : number of colors of the colormap, two curves are used. The number of levels of the 
                                 colormap needs to be higher than 2 
            - dpi              : dots per inch of the saved figure
            - uu_struc         : array of the streamwise velocity
            - vv_struc         : array of the wall-normal velocity
            - ww_struc         : array of the spanwise velocity
            - yplus_struc      : array of the wall-distance 
            - yplus_mesh       : positions of the mesh in the wall-normal direction
            - bins             : bins in the pdf
            - lev_min          : minimum value of the levels of the pdf
            - lev_delta        : distance between levels (None if distance needs to be calculated by the code)
            - linewidth        : width of the line
            - umin             : minimum value of the streamwise velocity in the histogram
            - umax             : maximum value of the streamwise velocity in the histogram
            - vmin             : minimum value of the wall-normal velocity in the histogram
            - vmax             : maximum value of the wall-normal velocity in the histogram
            - wmin             : minimum value of the spanwise velocity in the histogram
            - wmax             : maximum value of the spanwise velocity in the histogram

    Returns
    -------
    None.

    """
    
    
    # -------------------------------------------------------------------------------------------------------------------
    # Import packages
    # -------------------------------------------------------------------------------------------------------------------
    from py_bin.py_class.plot_format import plot_format
    from matplotlib import ticker
    
    # -------------------------------------------------------------------------------------------------------------------
    # Read the parameters of the plot
    # -------------------------------------------------------------------------------------------------------------------
    plot_folder           = str(data_in["plot_folder"])
    plot_fileu            = str(data_in["plot_fileu"])
    plot_filev            = str(data_in["plot_filev"])
    plot_filew            = str(data_in["plot_filew"])
    ylabel                = str(data_in["ylabel"])
    xlabelu               = str(data_in["xlabelu"])
    xlabelv               = str(data_in["xlabelv"])
    xlabelw               = str(data_in["xlabelw"])
    fontsize              = int(data_in["fontsize"])
    figsize_x             = int(data_in["figsize_x"])
    figsize_y             = int(data_in["figsize_y"])
    colormap              = str(data_in["colormap"])
    colornum              = int(data_in["colornum"])
    dpi                   = float(data_in["dpi"])
    uu_struc              = np.array(data_in["uu_struc"],dtype="float")
    vv_struc              = np.array(data_in["vv_struc"],dtype="float")
    ww_struc              = np.array(data_in["ww_struc"],dtype="float")
    yplus_struc           = np.array(data_in["yplus_struc"],dtype="float")
    yplusmesh             = np.array(data_in["yplusmesh"],dtype="float")
    bins                  = int(data_in["bins"])
    lev_min               = float(data_in["lev_min"])
    lev_delta             = data_in["lev_delta"]
    linewidth             = int(data_in["linewidth"])
    umin_in               = data_in["umin"]
    umax_in               = data_in["umax"]
    vmin_in               = data_in["vmin"]
    vmax_in               = data_in["vmax"]
    wmin_in               = data_in["wmin"]
    wmax_in               = data_in["wmax"]
    
    # -------------------------------------------------------------------------------------------------------------------
    # Create the mesh for the pdf
    # -------------------------------------------------------------------------------------------------------------------
    diffy = np.diff(yplusmesh)/2
    binsy = np.concatenate(([-diffy[0]],yplusmesh[:-1]+diffy,[yplusmesh[-1]+diffy[-1]]))
      
    # ----------------------------------------------------------------------------------------------------------------
    # Calculate the number of levels
    # ----------------------------------------------------------------------------------------------------------------
    if lev_delta is None:
        explev0   = np.log10(lev_min)
        explevs   = np.linspace(explev0,0,colornum)
        levels    = 10**explevs 
        locator   = ticker.LogLocator(numticks=colornum)
    else:
        lev_delta = int(lev_delta)
        levels    = [lev_min,lev_min+lev_delta,lev_min+2*lev_delta,lev_min+3*lev_delta,lev_min+4*lev_delta]
        locator   = ticker.LinearLocator()
        
    # ----------------------------------------------------------------------------------------------------------------
    # Calculate the pdf for u
    # ----------------------------------------------------------------------------------------------------------------
    hist_uy,hist_u,hist_y = np.histogram2d(uu_struc,yplus_struc,bins=(bins,binsy))
    hist_y                = hist_y[:-1]+np.diff(hist_y)/2
    hist_u                = hist_u[:-1]+np.diff(hist_u)/2
    grid_u,grid_y         = np.meshgrid(hist_u,hist_y)
    grid_uy               = hist_uy.T.copy()
    grid_uy              /= np.max(grid_uy)
    ucontent              = grid_u[np.where(grid_uy>=lev_min)]
    if umin_in is None:
        umin              = np.min(ucontent)
    else:
        umin              = float(umin_in)
    if umax_in is None:
        umax              = np.max(ucontent) 
    else:
        umax              = float(umax_in)
    
    
    # -------------------------------------------------------------------------------------------------------------------
    # Create the plot for u
    # -------------------------------------------------------------------------------------------------------------------
    titlefig   = "$u^+ - y^+$"
    data_plot  = {"xlabel":xlabelu,"ylabel":ylabel,"zlabel":[],"fontsize":fontsize,"figsize_x":figsize_x,
                  "figsize_y":figsize_y,"xscale":"linear","yscale":"log","zscale":"linear","colormap":colormap,
                  "colornum":colornum,"legend":True,"fig_name":plot_fileu,"dpi":dpi,"plot_folder":plot_folder,
                  "xmin":umin,"xmax":umax,"ymin":1,"ymax":np.max(grid_y),"zmin":None,"zmax":None}
    plot_pred = plot_format(data_in=data_plot)
    plot_pred.create_figure()
    try:
        plot_pred.add_hist2d_y(data_in={"xx":grid_u,"yy":grid_y,"xxyy":grid_uy,"levels":levels,"colormap":None,
                                        "locator":locator,"alp":0.65,"linewidth":linewidth})
        labels = ['{0:.1e}'.format(levels[ii]) for ii in np.arange(len(levels))]
        plot_pred.plot_layout_pcolor(data_in={"title":None,"colorbar":True,"b_text":None,
                                              "colorticks":levels,"colorlabels":labels,"equal":False,
                                              "xticks":None,"yticks":None,"xticklabels":None,"yticklabels":None})
    except:
        print("Error in contours",flush=True)
    try:
        os.mkdir(plot_folder)
    except:
        print("Existing folder...",flush=True)
    plot_pred.plot_save_png()
    plot_pred.plot_save_pdf()
    plot_pred.close()
   
    # ----------------------------------------------------------------------------------------------------------------
    # Calculate the pdf for v
    # ----------------------------------------------------------------------------------------------------------------
    hist_vy,hist_v,hist_y = np.histogram2d(vv_struc,yplus_struc,bins=(bins,binsy))
    hist_y                = hist_y[:-1]+np.diff(hist_y)/2
    hist_v                = hist_v[:-1]+np.diff(hist_v)/2
    grid_v,grid_y         = np.meshgrid(hist_v,hist_y)
    grid_vy               = hist_vy.T.copy()
    grid_vy              /= np.max(grid_vy)
    vcontent              = grid_v[np.where(grid_vy>=lev_min)]
    if vmin_in is None:
        vmin              = np.min(vcontent)
    else:
        vmin              = float(vmin_in)
    if vmax_in is None:
        vmax              = np.max(vcontent) 
    else:
        vmax              = float(vmax_in) 
    
    # -------------------------------------------------------------------------------------------------------------------
    # Create the plot for v
    # -------------------------------------------------------------------------------------------------------------------
    titlefig   = "$v^+ - y^+$"
    data_plot  = {"xlabel":xlabelv,"ylabel":ylabel,"zlabel":[],"fontsize":fontsize,"figsize_x":figsize_x,
                  "figsize_y":figsize_y,"xscale":"linear","yscale":"log","zscale":"linear","colormap":colormap,
                  "colornum":colornum,"legend":True,"fig_name":plot_filev,"dpi":dpi,"plot_folder":plot_folder,
                  "xmin":vmin,"xmax":vmax,"ymin":1,"ymax":np.max(grid_y),"zmin":None,"zmax":None}
    plot_pred = plot_format(data_in=data_plot)
    plot_pred.create_figure()
    try:
        plot_pred.add_hist2d_y(data_in={"xx":grid_v,"yy":grid_y,"xxyy":grid_vy,"levels":levels,"colormap":None,
                                      "locator":locator,"alp":0.65,"linewidth":linewidth})
        labels = ['{0:.1e}'.format(levels[ii]) for ii in np.arange(len(levels))]
        plot_pred.plot_layout_pcolor(data_in={"title":None,"colorbar":True,"b_text":None,
                                              "colorticks":levels,"colorlabels":labels,"equal":False,
                                              "xticks":None,"yticks":None,"xticklabels":None,"yticklabels":None})
    except:
        print("Error in contours",flush=True)
    try:
        os.mkdir(plot_folder)
    except:
        print("Existing folder...",flush=True)
    plot_pred.plot_save_png()
    plot_pred.plot_save_pdf()
    plot_pred.close()
    
    # ----------------------------------------------------------------------------------------------------------------
    # Calculate the pdf for u
    # ----------------------------------------------------------------------------------------------------------------
    hist_wy,hist_w,hist_y = np.histogram2d(ww_struc,yplus_struc,bins=(bins,binsy))
    hist_y                = hist_y[:-1]+np.diff(hist_y)/2
    hist_w                = hist_w[:-1]+np.diff(hist_w)/2
    grid_w,grid_y         = np.meshgrid(hist_w,hist_y)
    grid_wy               = hist_wy.T.copy()
    grid_wy              /= np.max(grid_wy)
    wcontent              = grid_w[np.where(grid_wy>=lev_min)]
    if wmin_in is None:
        wmin              = np.min(wcontent)
    else:
        wmin              = float(wmin_in)
    if wmax_in is None:
        wmax              = np.max(wcontent) 
    else:
        wmax              = float(wmax_in)
    
    # -------------------------------------------------------------------------------------------------------------------
    # Create the plot for w
    # -------------------------------------------------------------------------------------------------------------------
    titlefig   = "$v^+ - y^+$"
    data_plot  = {"xlabel":xlabelw,"ylabel":ylabel,"zlabel":[],"fontsize":fontsize,"figsize_x":figsize_x,
                  "figsize_y":figsize_y,"xscale":"linear","yscale":"log","zscale":"linear","colormap":colormap,
                  "colornum":colornum,"legend":True,"fig_name":plot_filew,"dpi":dpi,"plot_folder":plot_folder,
                  "xmin":wmin,"xmax":wmax,"ymin":1,"ymax":np.max(grid_y),"zmin":None,"zmax":None}
    plot_pred = plot_format(data_in=data_plot)
    plot_pred.create_figure()
    try:
        plot_pred.add_hist2d_y(data_in={"xx":grid_w,"yy":grid_y,"xxyy":grid_wy,"levels":levels,"colormap":None,
                                      "locator":locator,"alp":0.65,"linewidth":linewidth})
        labels = ['{0:.1e}'.format(levels[ii]) for ii in np.arange(len(levels))]
        plot_pred.plot_layout_pcolor(data_in={"title":None,"colorbar":True,"b_text":None,
                                              "colorticks":levels,"colorlabels":labels,"equal":False,
                                              "xticks":None,"yticks":None,"xticklabels":None,"yticklabels":None})
    except:
        print("Error in contours",flush=True)
    try:
        os.mkdir(plot_folder)
    except:
        print("Existing folder...",flush=True)
    plot_pred.plot_save_png()
    plot_pred.plot_save_pdf()
    plot_pred.close()
    

    
def plot_histuvw_y_lowmem(data_in={"plot_folder":"plots","plot_fileu":"fileu","plot_filev":"filev","plot_filew":"filew",
                                   "ylabel":"y","xlabelu":"u","xlabelv":"v","xlabelw":"w","fontsize":18,"figsize_x":10,
                                   "figsize_y":8,"colormap":"viridis","colornum":2,"dpi":60,"grid_uy":[],
                                   "grid_vy":[],"grid_wy":[],"grid_y":[],"grid_u":[],"grid_v":[],
                                   "grid_w":[],"lev_min":1e-2,"lev_delta":None,"linewidth":2,"umin":None,
                                   "umax":None,"vmin":None,"vmax":None,"wmin":None,"wmax":None,"ymax":125}):
    """""
    .....................................................................................................................
    # plot_histuvw_y_lowmem: Function to plot the pdf of the velocities in the structures for low memory consumption
    .....................................................................................................................
    Parameters
    ----------
    data_in : dict, optional
        Data required for generating the plot. 
        The default is {"plot_folder":"plots","plot_fileu":"fileu","plot_filev":"filev","plot_filew":"filew",
                        "ylabel":"y","xlabelu":"u","xlabelv":"v","xlabelw":"w","fontsize":18,"figsize_x":10,
                        "figsize_y":8,"colormap":"viridis","colornum":2,"dpi":60,"grid_uy":grid_uy,
                        "grid_vy":grid_vy,"grid_wy":grid_wy,"grid_y":grid_y,"grid_u":grid_u,"grid_v":grid_v,
                        "grid_w":grid_w,"lev_min":1e-2,"lev_delta":None,"linewidth":2,"umin":None,
                        "umax":None,"vmin":None,"vmax":None,"wmin":None,"wmax":None,"ymax":125}.
        Data:
            - plot_folder      : folder to store the plots
            - plot_fileu       : file to save the pdf of the streamwise velocity
            - plot_filev       : file to save the pdf of the wall-normal velocity
            - plot_filew       : file to save the pdf of the spanwise velocity
            - ylabel           : label of the y axis
            - xlabelu          : label of the x axis for the streamwise velocity
            - xlabelv          : label of the x axis for the wall-normal velocity
            - xlabelw          : label of the x axis for the spanwise velocity
            - fontsize         : font size used for the figure
            - figsize_x        : size of the figure in x
            - figsize_y        : size of the figure in y
            - colormap         : colormap used for the figure
            - colornum         : number of colors of the colormap, two curves are used. The number of levels of the 
                                 colormap needs to be higher than 2 
            - dpi              : dots per inch of the saved figure
            - grid_uy          : grid of the repetitions in u and y
            - grid_vy          : grid of the repetitions in v and y
            - grid_wy          : grid of the repetitions in w and y
            - grid_y           : grid of the repetitions in y
            - grid_u           : grid of the repetitions in u
            - grid_v           : grid of the repetitions in v
            - grid_w           : grid of the repetitions in w
            - lev_min          : minimum value of the levels of the pdf
            - lev_delta        : distance between levels (None if distance needs to be calculated by the code)
            - linewidth        : width of the line
            - umin             : minimum value of the streamwise velocity in the histogram
            - umax             : maximum value of the streamwise velocity in the histogram
            - vmin             : minimum value of the wall-normal velocity in the histogram
            - vmax             : maximum value of the wall-normal velocity in the histogram
            - wmin             : minimum value of the spanwise velocity in the histogram
            - wmax             : maximum value of the spanwise velocity in the histogram
            - ymax             : maximum y coordinate

    Returns
    -------
    None.

    """
    
    
    # -------------------------------------------------------------------------------------------------------------------
    # Import packages
    # -------------------------------------------------------------------------------------------------------------------
    from py_bin.py_class.plot_format import plot_format
    from matplotlib import ticker
    
    # -------------------------------------------------------------------------------------------------------------------
    # Read the parameters of the plot
    # -------------------------------------------------------------------------------------------------------------------
    plot_folder           = str(data_in["plot_folder"])
    plot_fileu            = str(data_in["plot_fileu"])
    plot_filev            = str(data_in["plot_filev"])
    plot_filew            = str(data_in["plot_filew"])
    ylabel                = str(data_in["ylabel"])
    xlabelu               = str(data_in["xlabelu"])
    xlabelv               = str(data_in["xlabelv"])
    xlabelw               = str(data_in["xlabelw"])
    fontsize              = int(data_in["fontsize"])
    figsize_x             = int(data_in["figsize_x"])
    figsize_y             = int(data_in["figsize_y"])
    colormap              = str(data_in["colormap"])
    colornum              = int(data_in["colornum"])
    dpi                   = float(data_in["dpi"])
    grid_uy               = np.array(data_in["grid_uy"],dtype="float")
    grid_vy               = np.array(data_in["grid_vy"],dtype="float")
    grid_wy               = np.array(data_in["grid_wy"],dtype="float")
    grid_y                = np.array(data_in["grid_y"],dtype="float")
    grid_u                = np.array(data_in["grid_u"],dtype="float")
    grid_v                = np.array(data_in["grid_v"],dtype="float")
    grid_w                = np.array(data_in["grid_w"],dtype="float")
    lev_min               = float(data_in["lev_min"])
    lev_delta             = data_in["lev_delta"]
    linewidth             = int(data_in["linewidth"])
    umin_in               = data_in["umin"]
    umax_in               = data_in["umax"]
    vmin_in               = data_in["vmin"]
    vmax_in               = data_in["vmax"]
    wmin_in               = data_in["wmin"]
    wmax_in               = data_in["wmax"]
    if "ymax" in data_in:
        ymax              = float(data_in["ymax"])
    else:
        ymax              = np.max(grid_y)
    
      
    # ----------------------------------------------------------------------------------------------------------------
    # Calculate the number of levels
    # ----------------------------------------------------------------------------------------------------------------
    if lev_delta is None:
        lev_delta = 5
    #     explev0   = np.log10(lev_min)
    #     explevs   = np.linspace(explev0,0,5)
    #     levels    = 10**explevs 
    #     locator   = ticker.LogLocator(numticks=5)
    # else:
    #     lev_delta = int(lev_delta)
    #     levels    = [lev_min,lev_min+lev_delta,lev_min+2*lev_delta,lev_min+3*lev_delta,lev_min+4*lev_delta]
    #     locator   = ticker.LinearLocator()
    explev0   = np.log10(lev_min)
    explevs   = np.linspace(explev0,0,lev_delta)
    levels    = 10**explevs 
    locator   = ticker.LogLocator(numticks=lev_delta)
        
    # ----------------------------------------------------------------------------------------------------------------
    # Calculate the pdf for u
    # ----------------------------------------------------------------------------------------------------------------
    grid_uy              /= np.max(grid_uy[grid_y>1])
    ucontent              = grid_u[np.where(grid_uy>=lev_min)]
    if umin_in is None:
        umin              = np.min(ucontent)
    else:
        umin              = float(umin_in)
    if umax_in is None:
        umax              = np.max(ucontent) 
    else:
        umax              = float(umax_in)
    
    
    # -------------------------------------------------------------------------------------------------------------------
    # Create the plot for u
    # -------------------------------------------------------------------------------------------------------------------
    titlefig   = "$u^+ - y^+$"
    data_plot  = {"xlabel":xlabelu,"ylabel":ylabel,"zlabel":[],"fontsize":fontsize,"figsize_x":figsize_x,
                  "figsize_y":figsize_y,"xscale":"linear","yscale":"log","zscale":"linear","colormap":colormap,
                  "colornum":colornum,"legend":True,"fig_name":plot_fileu,"dpi":dpi,"plot_folder":plot_folder,
                  "xmin":umin,"xmax":umax,"ymin":1,"ymax":ymax,"zmin":None,"zmax":None}
    plot_pred = plot_format(data_in=data_plot)
    plot_pred.create_figure()
    try:
        plot_pred.add_hist2d_y(data_in={"xx":grid_u,"yy":grid_y,"xxyy":grid_uy,"levels":levels,"colormap":None,
                                        "locator":locator,"alp":0.65,"linewidth":linewidth})
        labels = ['{0:.1e}'.format(levels[ii]) for ii in np.arange(len(levels))]
        plot_pred.plot_layout_pcolor(data_in={"title":None,"colorbar":True,"b_text":None,
                                              "colorticks":levels,"colorlabels":labels,"equal":False,
                                              "xticks":None,"yticks":None,"xticklabels":None,"yticklabels":None})
    except:
        print("Error in contours",flush=True)
    try:
        os.mkdir(plot_folder)
    except:
        print("Existing folder...",flush=True)
    plot_pred.plot_save_png()
    plot_pred.plot_save_pdf()
    plot_pred.close()
   
    # ----------------------------------------------------------------------------------------------------------------
    # Calculate the pdf for v
    # ----------------------------------------------------------------------------------------------------------------
    grid_vy              /= np.max(grid_vy[grid_y>1])
    vcontent              = grid_v[np.where(grid_vy>=lev_min)]
    if vmin_in is None:
        vmin              = np.min(vcontent)
    else:
        vmin              = float(vmin_in)
    if vmax_in is None:
        vmax              = np.max(vcontent) 
    else:
        vmax              = float(vmax_in) 
    
    # -------------------------------------------------------------------------------------------------------------------
    # Create the plot for v
    # -------------------------------------------------------------------------------------------------------------------
    titlefig   = "$v^+ - y^+$"
    data_plot  = {"xlabel":xlabelv,"ylabel":ylabel,"zlabel":[],"fontsize":fontsize,"figsize_x":figsize_x,
                  "figsize_y":figsize_y,"xscale":"linear","yscale":"log","zscale":"linear","colormap":colormap,
                  "colornum":colornum,"legend":True,"fig_name":plot_filev,"dpi":dpi,"plot_folder":plot_folder,
                  "xmin":vmin,"xmax":vmax,"ymin":1,"ymax":np.max(grid_y),"zmin":None,"zmax":None}
    plot_pred = plot_format(data_in=data_plot)
    plot_pred.create_figure()
    try:
        plot_pred.add_hist2d_y(data_in={"xx":grid_v,"yy":grid_y,"xxyy":grid_vy,"levels":levels,"colormap":None,
                                      "locator":locator,"alp":0.65,"linewidth":linewidth})
        labels = ['{0:.1e}'.format(levels[ii]) for ii in np.arange(len(levels))]
        plot_pred.plot_layout_pcolor(data_in={"title":None,"colorbar":True,"b_text":None,
                                              "colorticks":levels,"colorlabels":labels,"equal":False,
                                              "xticks":None,"yticks":None,"xticklabels":None,"yticklabels":None})
    except:
        print("Error in contours",flush=True)
    try:
        os.mkdir(plot_folder)
    except:
        print("Existing folder...",flush=True)
    plot_pred.plot_save_png()
    plot_pred.plot_save_pdf()
    plot_pred.close()
    
    # ----------------------------------------------------------------------------------------------------------------
    # Calculate the pdf for u
    # ----------------------------------------------------------------------------------------------------------------
    grid_wy              /= np.max(grid_wy[grid_y>1])
    wcontent              = grid_w[np.where(grid_wy>=lev_min)]
    if wmin_in is None:
        wmin              = np.min(wcontent)
    else:
        wmin              = float(wmin_in)
    if wmax_in is None:
        wmax              = np.max(wcontent) 
    else:
        wmax              = float(wmax_in)
    
    # -------------------------------------------------------------------------------------------------------------------
    # Create the plot for w
    # -------------------------------------------------------------------------------------------------------------------
    titlefig   = "$v^+ - y^+$"
    data_plot  = {"xlabel":xlabelw,"ylabel":ylabel,"zlabel":[],"fontsize":fontsize,"figsize_x":figsize_x,
                  "figsize_y":figsize_y,"xscale":"linear","yscale":"log","zscale":"linear","colormap":colormap,
                  "colornum":colornum,"legend":True,"fig_name":plot_filew,"dpi":dpi,"plot_folder":plot_folder,
                  "xmin":wmin,"xmax":wmax,"ymin":1,"ymax":np.max(grid_y),"zmin":None,"zmax":None}
    plot_pred = plot_format(data_in=data_plot)
    plot_pred.create_figure()
    try:
        plot_pred.add_hist2d_y(data_in={"xx":grid_w,"yy":grid_y,"xxyy":grid_wy,"levels":levels,"colormap":None,
                                      "locator":locator,"alp":0.65,"linewidth":linewidth})
        labels = ['{0:.1e}'.format(levels[ii]) for ii in np.arange(len(levels))]
        plot_pred.plot_layout_pcolor(data_in={"title":None,"colorbar":True,"b_text":None,
                                              "colorticks":levels,"colorlabels":labels,"equal":False,
                                              "xticks":None,"yticks":None,"xticklabels":None,"yticklabels":None})
    except:
        print("Error in contours",flush=True)
    try:
        os.mkdir(plot_folder)
    except:
        print("Existing folder...",flush=True)
    plot_pred.plot_save_png()
    plot_pred.plot_save_pdf()
    plot_pred.close()
    
def plot_histuvw_y_lowmem_combined(data_in={"plot_folder":"plots","plot_fileu":"fileu","plot_filev":"filev","plot_filew":"filew",
                                            "ylabel":"y","xlabelu":"u","xlabelv":"v","xlabelw":"w","fontsize":18,"figsize_x":10,
                                            "figsize_y":8,"colormap":"viridis","colornum":2,"dpi":60,"grid_uy1":[],
                                            "grid_vy1":[],"grid_wy1":[],"grid_y1":[],"grid_u1":[],"grid_v1":[],
                                            "grid_w1":[],"grid_uy2":[],"grid_vy2":[],"grid_wy2":[],"grid_y2":[],
                                            "grid_u2":[],"grid_v2":[],"grid_w2":[],"lev_min":1e-2,"lev_delta":None,
                                            "linewidth":2,"umin":None,"umax":None,"vmin":None,"vmax":None,"wmin":None,
                                            "wmax":None,"ymax":125}):
    """""
    .....................................................................................................................
    # plot_histuvw_y_lowmem_combined: Function to plot the pdf of the velocities in the structures for low memory consumption
    .....................................................................................................................
    Parameters
    ----------
    data_in : dict, optional
        Data required for generating the plot. 
        The default is {"plot_folder":"plots","plot_fileu":"fileu","plot_filev":"filev","plot_filew":"filew",
                        "ylabel":"y","xlabelu":"u","xlabelv":"v","xlabelw":"w","fontsize":18,"figsize_x":10,
                        "figsize_y":8,"colormap":"viridis","colornum":2,"dpi":60,"grid_uy1":[],
                        "grid_vy1":[],"grid_wy1":[],"grid_y1":[],"grid_u1":[],"grid_v1":[],
                        "grid_w1":[],"grid_uy2":[],"grid_vy2":[],"grid_wy2":[],"grid_y2":[],
                        "grid_u2":[],"grid_v2":[],"grid_w2":[],"lev_min":1e-2,"lev_delta":None,
                        "linewidth":2,"umin":None,"umax":None,"vmin":None,"vmax":None,"wmin":None,
                        "wmax":None,"ymax":125}.
        Data:
            - plot_folder      : folder to store the plots
            - plot_fileu       : file to save the pdf of the streamwise velocity
            - plot_filev       : file to save the pdf of the wall-normal velocity
            - plot_filew       : file to save the pdf of the spanwise velocity
            - ylabel           : label of the y axis
            - xlabelu          : label of the x axis for the streamwise velocity
            - xlabelv          : label of the x axis for the wall-normal velocity
            - xlabelw          : label of the x axis for the spanwise velocity
            - fontsize         : font size used for the figure
            - figsize_x        : size of the figure in x
            - figsize_y        : size of the figure in y
            - colormap         : colormap used for the figure
            - colornum         : number of colors of the colormap, two curves are used. The number of levels of the 
                                 colormap needs to be higher than 2 
            - dpi              : dots per inch of the saved figure
            - grid_uy1         : grid of the repetitions in u and y for case 1
            - grid_vy1         : grid of the repetitions in v and y for case 1
            - grid_wy1         : grid of the repetitions in w and y for case 1
            - grid_y1          : grid of the repetitions in y for case 1
            - grid_u1          : grid of the repetitions in u for case 1
            - grid_v1          : grid of the repetitions in v for case 1
            - grid_w1          : grid of the repetitions in w for case 1
            - grid_uy2         : grid of the repetitions in u and y for case 2
            - grid_vy2         : grid of the repetitions in v and y for case 2
            - grid_wy2         : grid of the repetitions in w and y for case 2
            - grid_y2          : grid of the repetitions in y for case 2
            - grid_u2          : grid of the repetitions in u for case 2
            - grid_v2          : grid of the repetitions in v for case 2
            - grid_w2          : grid of the repetitions in w for case 2
            - lev_min          : minimum value of the levels of the pdf
            - lev_delta        : distance between levels (None if distance needs to be calculated by the code)
            - linewidth        : width of the line
            - umin             : minimum value of the streamwise velocity in the histogram
            - umax             : maximum value of the streamwise velocity in the histogram
            - vmin             : minimum value of the wall-normal velocity in the histogram
            - vmax             : maximum value of the wall-normal velocity in the histogram
            - wmin             : minimum value of the spanwise velocity in the histogram
            - wmax             : maximum value of the spanwise velocity in the histogram
            - ymax             : maximum y coordinate

    Returns
    -------
    None.

    """
    
    
    # -------------------------------------------------------------------------------------------------------------------
    # Import packages
    # -------------------------------------------------------------------------------------------------------------------
    from py_bin.py_class.plot_format import plot_format
    from matplotlib import ticker
    
    # -------------------------------------------------------------------------------------------------------------------
    # Read the parameters of the plot
    # -------------------------------------------------------------------------------------------------------------------
    plot_folder           = str(data_in["plot_folder"])
    plot_fileu            = str(data_in["plot_fileu"])
    plot_filev            = str(data_in["plot_filev"])
    plot_filew            = str(data_in["plot_filew"])
    ylabel                = str(data_in["ylabel"])
    xlabelu               = str(data_in["xlabelu"])
    xlabelv               = str(data_in["xlabelv"])
    xlabelw               = str(data_in["xlabelw"])
    fontsize              = int(data_in["fontsize"])
    figsize_x             = int(data_in["figsize_x"])
    figsize_y             = int(data_in["figsize_y"])
    colormap              = str(data_in["colormap"])
    colornum              = int(data_in["colornum"])
    dpi                   = float(data_in["dpi"])
    grid_uy1              = np.array(data_in["grid_uy1"],dtype="float")
    grid_vy1              = np.array(data_in["grid_vy1"],dtype="float")
    grid_wy1              = np.array(data_in["grid_wy1"],dtype="float")
    grid_y1               = np.array(data_in["grid_y1"],dtype="float")
    grid_u1               = np.array(data_in["grid_u1"],dtype="float")
    grid_v1               = np.array(data_in["grid_v1"],dtype="float")
    grid_w1               = np.array(data_in["grid_w1"],dtype="float")
    grid_uy2              = np.array(data_in["grid_uy2"],dtype="float")
    grid_vy2              = np.array(data_in["grid_vy2"],dtype="float")
    grid_wy2              = np.array(data_in["grid_wy2"],dtype="float")
    grid_y2               = np.array(data_in["grid_y2"],dtype="float")
    grid_u2               = np.array(data_in["grid_u2"],dtype="float")
    grid_v2               = np.array(data_in["grid_v2"],dtype="float")
    grid_w2               = np.array(data_in["grid_w2"],dtype="float")
    lev_min               = float(data_in["lev_min"])
    lev_delta             = data_in["lev_delta"]
    linewidth             = int(data_in["linewidth"])
    umin_in               = data_in["umin"]
    umax_in               = data_in["umax"]
    vmin_in               = data_in["vmin"]
    vmax_in               = data_in["vmax"]
    wmin_in               = data_in["wmin"]
    wmax_in               = data_in["wmax"]
    if "ymax" in data_in:
        ymax              = float(data_in["ymax"])
    else:
        ymax              = np.max(grid_y)
    
      
    # ----------------------------------------------------------------------------------------------------------------
    # Calculate the number of levels
    # ----------------------------------------------------------------------------------------------------------------
    if lev_delta is None:
        lev_delta = 5
    #     explev0   = np.log10(lev_min)
    #     explevs   = np.linspace(explev0,0,5)
    #     levels    = 10**explevs 
    #     locator   = ticker.LogLocator(numticks=5)
    # else:
    #     lev_delta = int(lev_delta)
    #     levels    = [lev_min,lev_min+lev_delta,lev_min+2*lev_delta,lev_min+3*lev_delta,lev_min+4*lev_delta]
    #     locator   = ticker.LinearLocator()
    explev0   = np.log10(lev_min)
    explevs   = np.linspace(explev0,0,lev_delta)
    levels    = 10**explevs 
    locator   = ticker.LogLocator(numticks=lev_delta)
        
    # ----------------------------------------------------------------------------------------------------------------
    # Calculate the pdf for u
    # ----------------------------------------------------------------------------------------------------------------
    grid_uy1             /= np.max(grid_uy1[grid_y1>1])
    grid_uy2             /= np.max(grid_uy2[grid_y2>1])
    ucontent1             = grid_u1[np.where(grid_uy1>=lev_min)]
    if umin_in is None:
        umin              = np.min(ucontent1)
    else:
        umin              = float(umin_in)
    if umax_in is None:
        umax              = np.max(ucontent1) 
    else:
        umax              = float(umax_in)
    
    
    # -------------------------------------------------------------------------------------------------------------------
    # Create the plot for u
    # -------------------------------------------------------------------------------------------------------------------
    titlefig   = "$u^+ - y^+$"
    data_plot  = {"xlabel":xlabelu,"ylabel":ylabel,"zlabel":[],"fontsize":fontsize,"figsize_x":figsize_x,
                  "figsize_y":figsize_y,"xscale":"linear","yscale":"log","zscale":"linear","colormap":colormap,
                  "colornum":colornum,"legend":True,"fig_name":plot_fileu,"dpi":dpi,"plot_folder":plot_folder,
                  "xmin":umin,"xmax":umax,"ymin":1,"ymax":ymax,"zmin":None,"zmax":None}
    plot_pred = plot_format(data_in=data_plot)
    plot_pred.create_figure()
    try:
        plot_pred.add_hist2d_y_cont(data_in={"xx":grid_u1,"yy":grid_y1,"xxyy":grid_uy1,"levels":levels,"colormap":colormap,
                                             "locator":locator,"alp":0.65,"linewidth":linewidth,"cmapflag":True})
        plot_pred.add_hist2d_y_cont(data_in={"xx":grid_u2,"yy":grid_y2,"xxyy":grid_uy2,"levels":levels,"colormap":colormap,
                                             "locator":locator,"alp":0.65,"linewidth":linewidth,"cmapflag":True,"linestyle":"--"})
        labels = ['{0:.1e}'.format(levels[ii]) for ii in np.arange(len(levels))]
        plot_pred.plot_layout_pcolor(data_in={"title":None,"colorbar":True,"b_text":None,
                                              "colorticks":levels,"colorlabels":labels,"equal":False,
                                              "xticks":None,"yticks":None,"xticklabels":None,"yticklabels":None})
    except:
        print("Error in contours",flush=True)
    try:
        os.mkdir(plot_folder)
    except:
        print("Existing folder...",flush=True)
    plot_pred.plot_save_png()
    plot_pred.plot_save_pdf()
    plot_pred.close()
   
    # ----------------------------------------------------------------------------------------------------------------
    # Calculate the pdf for v
    # ----------------------------------------------------------------------------------------------------------------
    grid_vy1             /= np.max(grid_vy1[grid_y1>1])
    grid_vy2             /= np.max(grid_vy2[grid_y2>1])
    vcontent1             = grid_v1[np.where(grid_vy1>=lev_min)]
    if vmin_in is None:
        vmin              = np.min(vcontent1)
    else:
        vmin              = float(vmin_in)
    if vmax_in is None:
        vmax              = np.max(vcontent1) 
    else:
        vmax              = float(vmax_in) 
    
    # -------------------------------------------------------------------------------------------------------------------
    # Create the plot for v
    # -------------------------------------------------------------------------------------------------------------------
    titlefig   = "$v^+ - y^+$"
    data_plot  = {"xlabel":xlabelv,"ylabel":ylabel,"zlabel":[],"fontsize":fontsize,"figsize_x":figsize_x,
                  "figsize_y":figsize_y,"xscale":"linear","yscale":"log","zscale":"linear","colormap":colormap,
                  "colornum":colornum,"legend":True,"fig_name":plot_filev,"dpi":dpi,"plot_folder":plot_folder,
                  "xmin":vmin,"xmax":vmax,"ymin":1,"ymax":np.max(grid_y1),"zmin":None,"zmax":None}
    plot_pred = plot_format(data_in=data_plot)
    plot_pred.create_figure()
    try:
        plot_pred.add_hist2d_y_cont(data_in={"xx":grid_u1,"yy":grid_y1,"xxyy":grid_vy1,"levels":levels,"colormap":colormap,
                                             "locator":locator,"alp":0.65,"linewidth":linewidth,"cmapflag":True})
        plot_pred.add_hist2d_y_cont(data_in={"xx":grid_u2,"yy":grid_y2,"xxyy":grid_vy2,"levels":levels,"colormap":colormap,
                                             "locator":locator,"alp":0.65,"linewidth":linewidth,"cmapflag":True,"linestyle":"--"})
        labels = ['{0:.1e}'.format(levels[ii]) for ii in np.arange(len(levels))]
        plot_pred.plot_layout_pcolor(data_in={"title":None,"colorbar":True,"b_text":None,
                                              "colorticks":levels,"colorlabels":labels,"equal":False,
                                              "xticks":None,"yticks":None,"xticklabels":None,"yticklabels":None})
    except:
        print("Error in contours",flush=True)
    try:
        os.mkdir(plot_folder)
    except:
        print("Existing folder...",flush=True)
    plot_pred.plot_save_png()
    plot_pred.plot_save_pdf()
    plot_pred.close()
    
    # ----------------------------------------------------------------------------------------------------------------
    # Calculate the pdf for u
    # ----------------------------------------------------------------------------------------------------------------
    grid_wy1             /= np.max(grid_wy1[grid_y1>1])
    grid_wy2             /= np.max(grid_wy2[grid_y2>1])
    wcontent1             = grid_w1[np.where(grid_wy1>=lev_min)]
    if wmin_in is None:
        wmin              = np.min(wcontent1)
    else:
        wmin              = float(wmin_in)
    if wmax_in is None:
        wmax              = np.max(wcontent1) 
    else:
        wmax              = float(wmax_in)
    
    # -------------------------------------------------------------------------------------------------------------------
    # Create the plot for w
    # -------------------------------------------------------------------------------------------------------------------
    titlefig   = "$v^+ - y^+$"
    data_plot  = {"xlabel":xlabelw,"ylabel":ylabel,"zlabel":[],"fontsize":fontsize,"figsize_x":figsize_x,
                  "figsize_y":figsize_y,"xscale":"linear","yscale":"log","zscale":"linear","colormap":colormap,
                  "colornum":colornum,"legend":True,"fig_name":plot_filew,"dpi":dpi,"plot_folder":plot_folder,
                  "xmin":wmin,"xmax":wmax,"ymin":1,"ymax":np.max(grid_y1),"zmin":None,"zmax":None}
    plot_pred = plot_format(data_in=data_plot)
    plot_pred.create_figure()
    try:
        plot_pred.add_hist2d_y_cont(data_in={"xx":grid_u1,"yy":grid_y1,"xxyy":grid_wy1,"levels":levels,"colormap":colormap,
                                             "locator":locator,"alp":0.65,"linewidth":linewidth,"cmapflag":True})
        plot_pred.add_hist2d_y_cont(data_in={"xx":grid_u2,"yy":grid_y2,"xxyy":grid_wy2,"levels":levels,"colormap":colormap,
                                             "locator":locator,"alp":0.65,"linewidth":linewidth,"cmapflag":True,"linestyle":"--"})
        labels = ['{0:.1e}'.format(levels[ii]) for ii in np.arange(len(levels))]
        plot_pred.plot_layout_pcolor(data_in={"title":None,"colorbar":True,"b_text":None,
                                              "colorticks":levels,"colorlabels":labels,"equal":False,
                                              "xticks":None,"yticks":None,"xticklabels":None,"yticklabels":None})
    except:
        print("Error in contours",flush=True)
    try:
        os.mkdir(plot_folder)
    except:
        print("Existing folder...",flush=True)
    plot_pred.plot_save_png()
    plot_pred.plot_save_pdf()
    plot_pred.close()