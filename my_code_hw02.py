#-- my_code_hw02.py
#-- Assignment 02 GEO1015.2020
#-- Michiel de Jong  
#-- 4376978 
 

import sys
import math
import numpy
import rasterio
from rasterio import features

def distance(point1, point2):
	    return math.sqrt(((point2[0]-point1[0])*(point2[0]-point1[0]))+((point2[1]-point1[1])*(point2[1]-point1[1])))

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

    #-- numpy of input
    npi  = d.read(1)
    npvs = numpy.full(d.shape,3 ,dtype=numpy.int8)

    for num, each in enumerate(viewpoints):
        #-- fetch the 1st viewpoint
        v = viewpoints[num]
        #-- index of this point in the numpy raster
        vrow, vcol = d.index(v[0], v[1])
        #-- the results of the viewshed in npvs, all values=0
        #-- put that pixel with value 2
        vind = (vrow, vcol)
        #calculate raster in index space
        test_row, test_col = d.index((v[0]+maxdistance), v[1])
        r = test_col - vcol
        circle_cells = []

        #construct the circle
        for angle in numpy.arange(0, 360, 1/len(npi)):
            x = r * math.sin(math.radians(angle)) + vrow
            y = r * math.cos(math.radians(angle)) + vcol
            circle_cells.append((int(round(x)),int(round(y))))
        
        #remove duplicates
        circle_cell_final = list(set(circle_cells))

        #perform LoS queries on all cells in circle
        for cell in circle_cell_final:
            XY = Bresenham_with_rasterio(d, vind, cell)
            
            #initial height is height in DEM + eyeheight from .json
            init_height = npi[vind[0],vind[1]] +v[2] 
            
            #calculate delta x for initial tangent
            cell_1 = d.xy(XY[1][0],XY[1][1])
            dis_init = distance(v, cell_1)

            #calculate initial slope and initial tangent
            slope_init = ((npi[XY[1]])-init_height) / dis_init
            tcur = (slope_init * dis_init) + init_height
            
            #loop through all points on line, except those that have already been seen
            for i, pos in enumerate(XY):
                if npvs[XY[i]] != 1 and npvs[XY[i]] != 2:
                    if [XY[i]] != vind:
                        #get xy coords
                        point = d.xy(XY[i][0],XY[i][1])
                        #check tcur value for this pixel
                        tcur = (slope_init * distance(v, point)) + init_height                       
                        if tcur <= npi[XY[i]]:
                            #visible! yes!
                            npvs[XY[i]] = 1
                            #update slope and tangent
                            slope_init = (npi[XY[i]]- init_height) / distance(v,point)
                            tcur = (slope_init * distance(v, point)) + init_height                       
                        else:
                            #invisible
                            npvs[XY[i]] = 0
                    else:
                        continue
                else:
                    continue
        #make sure we can see the viewpoints
        npvs[vrow , vcol] = 2
        
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
    
    re = features.rasterize(shapes, 
                             out_shape=d.shape, 
                             all_touched=False,
                             transform=d.transform)
    #find where array is true
    where = numpy.argwhere(re==1)
    line = []

    #sort according to which quarter of circle
    for i in where:
        line.append(tuple(i))
    if a[0]<b[0] and a[1]<=b[1]:
        line = sorted(line, key = lambda k: (k[0], k[1]))
    elif a[0]>b[0] and a[1]<=b[1]:
        line = sorted(line, key = lambda k: (-k[0], k[1]))
    elif a[0]<=b[0] and a[1]>b[1]:
        line = sorted(line, key = lambda k: (k[0], -k[1]))
    elif a[0]>b[0] and a[1]>=b[1]:
        line = sorted(line, key = lambda k: (-k[0], -k[1]))
    return line

