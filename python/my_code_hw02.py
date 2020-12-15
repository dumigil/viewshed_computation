#-- my_code_hw02.py
#-- Assignment 02 GEO1015.2020
#-- [YOUR NAME] 
#-- [YOUR STUDENT NUMBER] 
#-- [YOUR NAME] 
#-- [YOUR STUDENT NUMBER] 

import sys
import math
import numpy
import rasterio
from rasterio import features


def output_viewshed(d, viewpoints, maxdistance, output_file):
    """
    !!! TO BE COMPLETED !!!
     
    Function that writes the output raster
     
    Input:
        d:            the input datasets (rasterio format)  
        viewpoints:   a list of the viewpoints (x, y, height)
        maxdistance:  max distance one can see
        output_file:  path of the file to write as output
        
    Output:
        none (but output GeoTIFF file written to 'output-file')
    """  
    
    # [this code can and should be removed/modified/reutilised]
    # [it's just there to help you]

    #-- numpy of input
    npi  = d.read(1)
    #-- fetch the 1st viewpoint
    v = viewpoints[0]
    #-- index of this point in the numpy raster
    vrow, vcol = d.index(v[0], v[1])
    #-- the results of the viewshed in npvs, all values=0
    npvs = numpy.full(d.shape,3 ,dtype=numpy.int8)
    #-- put that pixel with value 2
    npvs[vrow , vcol] = 2
    vind = (vrow, vcol)
    print(npi[300][300])
    test_row, test_col = d.index((v[0]+maxdistance), v[1])
    print(test_row, test_col)
    r = test_col - vcol
    circle_cells = []
    for angle in numpy.arange(0, 360, 0.5):
        x = r * math.sin(math.radians(angle)) + vcol
        y = r * math.cos(math.radians(angle)) + vrow
        npvs[int(round(y))][int(round(x))] = 2
        circle_cells.append((int(round(x)),int(round(y))))
    circle_cell_final = list(set(circle_cells))

    '''for y in range(d.shape[1]):
        for x in range(d.shape[0]):
            if abs((x-vrow)**2 + (y-vcol)**2 - r**2) == EPSILON**2:
                circle_cells.append((x,y))
                npvs[x,y] = 2'''
    
    for cell in circle_cell_final:
        line_pq = Bresenham_with_rasterio(d, vind, cell)
        walkline = []
        x, y = ((numpy.where(line_pq)))
        XY = [i for i in zip(x, y)]
        if XY[0] != vind:
            XY = numpy.flipud(XY)
        for pixel in XY:
            walkline.append((npi[pixel[0]][pixel[1]]))
        print((XY[1]))
                
    #-- write this to disk
    with rasterio.open(output_file, 'w', 
                       driver='GTiff', 
                       height=npi.shape[0],
                       width=npi.shape[1], 
                       count=1, 
                       dtype=rasterio.uint8,
                       crs=d.crs, 
                       transform=d.transform) as dst:
        dst.write(npvs.astype(rasterio.uint8), 1)

    print("Viewshed file written to '%s'" % output_file)

def viewshed_circle(d,v,maxdistance):
    center = v
    radius = maxdistance / d.shape[0]
    height = d.shape[0]
    width = d.shape[1]
    Y, X = numpy.ogrid[:height,:width]
    dist_from_center = numpy.sqrt((X - center[0])**2 + (Y - center[1])**2)
    mask = dist_from_center <= radius
    numpy.savetxt('mask.txt', mask)
    return mask



def Bresenham_with_rasterio(d, v, q):
    a = v
    b = q
    #-- create in-memory a simple GeoJSON LineString
    v = {}
    v["type"] = "LineString"
    v["coordinates"] = []
    v["coordinates"].append(d.xy(a[0], a[1]))
    v["coordinates"].append(d.xy(b[0], b[1]))
    shapes = [(v, 1)]
    line_list = []
    re = features.rasterize(shapes, 
                             out_shape=d.shape, 
                             all_touched=True,
                             transform=d.transform)
    # re is a numpy with d.shape where the line is rasterised (values != 0)
    return re

