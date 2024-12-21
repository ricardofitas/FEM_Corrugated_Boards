# Created on 26.07.2024 at TU Darmstadt, Papierfabrikation department; Author @akram.metar@stud.tu-darmstadt.de
# Salome mecca 2022 version; 2d model; complete version
# Simulation of general geometry in place of flute in corrugated cardboards in Flat crush Test

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ========================================================
# ======== Geometry parameterisation =====================
# ========================================================
x_cordinates = [0,1,2,3,4,5]
y_cordinates = [0,0,1,1,0,0]
amplitude = 1.0
flute_thickness = 0.2
liner_thickness = 0.2



#----------------------------------------------
#----------------------------------------------

# Importing important modules for calculations
from math import sin,radians,floor,ceil,sqrt,atan2,cos
import numpy as np




# creating the original values for x and y cordinates
x_cordinates = np.array(x_cordinates)*amplitude
y_cordinates = np.array(y_cordinates)*amplitude
n_data_points = len(x_cordinates)


# creating the list for new x_cordinates and y_cordinates (offset from the original one)
new_x_cordinates = []
new_y_cordinates = []
# creating a new list of cordinates for thickned line

def calculate_offset_point(Ax, Ay, Bx, By, Cx, Cy, offset_distance):
    angle1 = atan2((By-Ay),(Bx-Ax))
    angle2 = atan2((Cy-By),(Cx-Bx))
    half_angle = (np.pi)+angle1-angle2
    mid_vector_length = offset_distance*cos(half_angle)
    other_angle = half_angle+angle1
    dx = mid_vector_length*cos(other_angle)
    dy = mid_vector_length*sin(other_angle)
    #print(np.rad2deg(half_angle))
    return (dx,dy)



if (n_data_points%2)==0:
    for count,value in enumerate (x_cordinates):
        if count==int(n_data_points-1) :
            x2 = x_cordinates[count]
            y2 = y_cordinates[count]
            new_x_cordinates.append(x2)
            new_y_cordinates.append(y2+flute_thickness)


        elif count==0 :
            x2 = x_cordinates[count]
            y2 = y_cordinates[count]
            new_x_cordinates.append(x2)
            new_y_cordinates.append(y2+flute_thickness)

        
        else:
            x1 = x_cordinates[count]
            x2 = x_cordinates[count+1]
            y1 = y_cordinates[count]
            y2 = y_cordinates[count+1]
            x0 = x_cordinates[count-1]
            y0 = y_cordinates[count-1]


            final_dx,final_dy = calculate_offset_point(x0,y0,x1,y1,x2,y2,flute_thickness)
            #final_dy = abs(final_dy)
            



            new_x_cordinates.append(x1-final_dy)
            new_y_cordinates.append(y1+final_dx)

else :
    for count,value in enumerate (x_cordinates):
        if count==int(n_data_points-1) :
            x2 = x_cordinates[count]
            y2 = y_cordinates[count]
            new_x_cordinates.append(x2)
            new_y_cordinates.append(y2+flute_thickness)


        elif count==0 :
            x2 = x_cordinates[count]
            y2 = y_cordinates[count]
            new_x_cordinates.append(x2)
            new_y_cordinates.append(y2+flute_thickness)

        elif count == int((n_data_points-1)/2):
            x2 = x_cordinates[count]
            y2 = y_cordinates[count]
            new_x_cordinates.append(x2)
            new_y_cordinates.append(y2+flute_thickness)
            
        else:
            x1 = x_cordinates[count]
            x2 = x_cordinates[count+1]
            y1 = y_cordinates[count]
            y2 = y_cordinates[count+1]
            x0 = x_cordinates[count-1]
            y0 = y_cordinates[count-1]


            tan0 = x2-x1
            tan1 = -(y2-y1)
            magnitude = np.sqrt(tan1**2+tan0**2)
            dx = tan1*flute_thickness/magnitude
            dy = tan0*flute_thickness/magnitude

            tan2 = x0-x1
            tan3 = (y0-y1)
            magnitude1 = np.sqrt(tan2**2+tan3**2)
            dx1 = tan3*flute_thickness/magnitude1
            dy1 = tan2*flute_thickness/magnitude1

            #final_dx = (dx+dx1)/2.0
            #final_dy = (dy+dy1)/2.0

            final_dx,final_dy = calculate_offset_point(x0,y0,x1,y1,x2,y2,flute_thickness)
            #final_dy = abs(final_dy)




            new_x_cordinates.append(x1+final_dx)
            new_y_cordinates.append(y1+final_dy)
        


# round the cordinate points to 6 decimal places
#x_cordinates = np.round(x_cordinates,6)
#y_cordinates = np.round(y_cordinates,6)
#new_x_cordinates = np.round(new_x_cordinates,6)
#new_y_cordinates = np.round(new_y_cordinates,6)

import sys
import salome
salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()

#=======================================================================
#============ Geometry module ==========================================
#=======================================================================

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS

geompy = geomBuilder.New()
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )

# creating the the points for the flute geometry from x and y_cordinate lists
for count,value in enumerate (x_cordinates):
    globals()[f'Vertex_{int(count+1)}'] = geompy.MakeVertex(x_cordinates[count],y_cordinates[count],0.0)
    #geompy.addToStudy(globals()[f'Vertex_{int(count+1)}'],'Vertex_'+str(count+1))

# creating the lines and joining them
all_lines = []
for count in range(1,len(x_cordinates)):
    globals()[f'Line_{int(count)}'] = geompy.MakeLineTwoPnt(globals()[f'Vertex_{int(count)}'],globals()[f'Vertex_{int(count+1)}'])
    #geompy.addToStudy(globals()[f'Line_{int(count)}'],'Line_'+str(count))
    all_lines.append(globals()[f'Line_{int(count)}'])

# Making the Fuse of all the lines
Single_flute_line = geompy.MakeFuseList(all_lines,True,True)
geompy.addToStudy(Single_flute_line,'Single_flute_line')


# creating the points for the top lines
for count,value in enumerate (new_x_cordinates):
    globals()[f't_Vertex_{int(count+1)}'] = geompy.MakeVertex(new_x_cordinates[count],new_y_cordinates[count],0.0)
    #geompy.addToStudy(globals()[f't_Vertex_{int(count+1)}'],'t_Vertex_'+str(count+1))

# creating the lines and joining them
all_lines_top = []
for count in range(1,len(new_x_cordinates)):
    globals()[f't_Line_{int(count)}'] = geompy.MakeLineTwoPnt(globals()[f't_Vertex_{int(count)}'],globals()[f't_Vertex_{int(count+1)}'])
    #geompy.addToStudy(globals()[f't_Line_{int(count)}'],'t_Line_'+str(count))
    all_lines_top.append(globals()[f't_Line_{int(count)}'])

# Making the Fuse of all the lines
Single_flute_line_top = geompy.MakeFuseList(all_lines_top,True,True)
geompy.addToStudy(Single_flute_line_top,'Single_flute_line_top')



def calculate_angle(Ax, Ay, Bx, By, Cx, Cy):
    """
    Calculate the angle ABC (in degrees) formed by points A, B, and C.
    
    Parameters:
    Ax, Ay -- Coordinates of point A.
    Bx, By -- Coordinates of point B.
    Cx, Cy -- Coordinates of point C.
    
    Returns:
    angle -- The angle ABC in degrees.
    """
    # Convert points to numpy arrays for vector calculations
    A = np.array([Ax, Ay])
    B = np.array([Bx, By])
    C = np.array([Cx, Cy])
    
    # Vectors BA and BC
    BA = A - B
    BC = C - B
    
    # Calculate the dot product and magnitudes of BA and BC
    dot_product = np.dot(BA, BC)
    magnitude_BA = np.linalg.norm(BA)
    magnitude_BC = np.linalg.norm(BC)
    
    # Calculate the cosine of the angle
    cos_angle = dot_product / (magnitude_BA * magnitude_BC)
    
    # Ensure the cosine value is within the valid range to avoid numerical errors
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    
    # Calculate the angle in radians and then convert to degrees
    angle_radians = np.arccos(cos_angle)
    angle_degrees = np.degrees(angle_radians)
    
    return angle_degrees