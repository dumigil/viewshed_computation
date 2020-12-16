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
    
    # [this code can and should be removed/modified/reutilised]
    # [it's just there to help you]

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
        print([vind])
        test_row, test_col = d.index((v[0]+maxdistance), v[1])
        r = test_col - vcol
        circle_cells = []
        for angle in numpy.arange(0, 360, 0.1):
            x = r * math.sin(math.radians(angle)) + vrow
            y = r * math.cos(math.radians(angle)) + vcol

            circle_cells.append((int(round(x)),int(round(y))))
        circle_cell_final = list(set(circle_cells))
        print(len(circle_cells))
        print(len(circle_cell_final))
        
        
        for cell in circle_cell_final:
            XY = Bresenham_with_rasterio(d, vind, cell)
            init_height = npi[vind[0],vind[1]] +v[2] 
            #x_new, y_new = ((numpy.where(line_pq)))
            #print(x_new,y_new)
            #XY = [i for i in zip(x_new, y_new)]
                      
            
            cell_1 = d.xy(XY[1][0],XY[1][1])
            dis_init = distance(v, cell_1)

            slope_init = ((npi[XY[1]])-init_height) / dis_init
            tcur = (slope_init * dis_init) + init_height
            print(npi[XY[1]])
            print(dis_init)
            print(slope_init)
            #print(XY[0],XY[-1])

            for i, pos in enumerate(XY):
                if npvs[XY[i]] != 1 and npvs[XY[i]] != 2:
                    if [XY[i]] != vind:
                        point = d.xy(XY[i][0],XY[i][1])
                        prev_point = d.xy(XY[i-1][0],XY[i-1][1])
                        dx = distance(point,prev_point)

                        tcur = (slope_init * distance(v, point)) + init_height                       
                        #print(tcur)
                        if tcur <= npi[XY[i]]:
                            npvs[XY[i]] = 1
                            slope_init = (npi[XY[i]]- init_height) / distance(v,point)
                            #tcur = (slope_init * distance(v, point)) + init_height                       
                        else:
                            npvs[XY[i]] = 0
                    else:
                        continue
                else:
                    continue
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
    #print(shapes)

    re = features.rasterize(shapes, 
                             out_shape=d.shape, 
                             all_touched=True,
                             transform=d.transform)
    # re is a numpy with d.shape where the line is rasterised (values != 0)

    where = numpy.argwhere(re==1)
    out = []
    for i in where:
        out.append(tuple(i))
    if a[0]>b[0] and a[1]<=b[1]:
        out = sorted(out, key = lambda k: (-k[0], k[1]))
    elif a[0]>b[0] and a[1]>=b[1]:
        out = sorted(out, key = lambda k: (-k[0], -k[1]))
    elif a[0]<=b[0] and a[1]>b[1]:
        out = sorted(out, key = lambda k: (k[0], -k[1]))        
    return out

