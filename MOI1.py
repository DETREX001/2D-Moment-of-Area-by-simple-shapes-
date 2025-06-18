# Libraries import
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np


# File Reading
input_file = open(r"A:\code\Pro-MOI\Input.txt","r")
Read_lines = input_file.read()
input_file.close()

#Data Processsing
Read_lines_list = Read_lines.splitlines()
Shapes_Dict = dict()

for i in Read_lines_list:
    if i[0] not in Shapes_Dict.keys():
        Shapes_Dict[i[0]] = [eval(i[1:])]
    
    elif i[0] in Shapes_Dict.keys():
        Shapes_Dict[i[0]] += [eval(i[1:])]

    else: pass
        
# Screen Des
Screen, ax = plt.subplots()
ax.axhline(y=0, color='black', linestyle='--', linewidth=1)
ax.axvline(x=0, color='black', linestyle='--', linewidth=1)


# Graph draw functions
def Rectangular_Calculations(T):
    Centroid = T[1]
    length = T[0][0]
    width = T[0][1]
    Rotation = T[2]

    angle_rad = np.radians(Rotation)
    Corner_coords = [Centroid[0]-width*0.5*np.cos(angle_rad)+length*0.5*np.sin(angle_rad),
                     Centroid[1]- length*0.5*np.cos(angle_rad)+width*0.5*np.sin(angle_rad)]
    return Corner_coords

def Rectangle(T): #T: ([l,b],[Centroid],rotation,state)
    global Screen, ax
    length =  T[0][0]
    width = T[0][1]
    corner_coords = Rectangular_Calculations(T)
    Rotation = T[2]
    state = T[3]
    if state == 1:
        rectangle = patches.Rectangle(tuple(corner_coords),width,length,angle=Rotation, alpha = 0.5, edgecolor= "black", facecolor = "blue")
    elif state == 0:
        rectangle = patches.Rectangle(tuple(corner_coords),width,length,angle=Rotation, edgecolor= "black", facecolor = "white")
    else:pass
        

    ax.add_patch(rectangle)

def Triangular_calculation(T):  # T : ([side Lengths:a,b,c],[Coordinates of first vertices],Angle formed between first and 2nd vertices wrt to x axis in anti-clockwise direction, state)
    a = T[0][0]
    b = T[0][1]
    c = T[0][2]
    origin = tuple(T[1])
    rotation = -T[2]

    x1,y1 = origin
    x2, y2 = x1 + a, y1

    cos_r = (b**2 + a**2 - c**2)/(2*a*b)
    sin_r = (1-cos_r**2)**0.5
    x3 = x1 + b* cos_r
    y3 = y1 + b* sin_r

    vertices_original = np.array([[x1,y1],[x2,y2],[x3,y3]])

    angle_rad = np.radians(rotation)
    rotation_mat = np.array([
        [np.cos(angle_rad), -np.sin(angle_rad)],
        [np.sin(angle_rad), np.cos(angle_rad)]
    ])

    vertices_rotated = np.dot(vertices_original, rotation_mat)
    diffrence = vertices_rotated[0] - origin
    final_vertices = vertices_rotated - diffrence
    final_vertices_list = final_vertices.tolist()
    return final_vertices_list

def Triangle(T):
    global Screen, ax
    a,b,c = Triangular_calculation(T)
    state = T[3]
    if state == 1:
        triangle = patches.Polygon([a,b,c], closed= True, edgecolor = "black", facecolor = "blue", alpha = 0.5)
    elif state == 0:
        triangle = patches.Polygon([a,b,c], closed= True, edgecolor = "black", facecolor = "white")
    else:pass
    
    ax.add_patch(triangle)

def Circle(T): #T:([Radius],[center],rotation,sector angle,state)
    global Screen, ax
    radius = T[0][0]
    center = tuple(T[1])
    rotation = T[2]
    sector_angle = T[3]
    state = T[4]

    if state == 1:
        sector = patches.Wedge(center=center, r= radius, theta1=rotation, theta2= rotation+sector_angle, edgecolor = "black", facecolor = "blue", alpha = 0.5)
    elif state == 0:
        sector = patches.Wedge(center=center, r= radius, theta1=rotation, theta2= rotation+sector_angle, edgecolor = "black", facecolor = "white")
    else: pass
    ax.add_patch(sector)

# Moment of area calculations
def Moment_of_inertia_Rectangle(T):
    length =  T[0][0]
    width = T[0][1]
    Center_coords = T[1]
    Rotation = T[2]
    theta = np.radians(Rotation)
    area = length*width
    state = T[3]
    
    Inertia_Matrix = np.array([
        [ width*(length**3)/12 , 0 ],
        [0, length*(width**3)/12]
    ])
    
    Rotation_Matrix = np.array([
        [np.cos(theta), np.sin(theta)],
        [-np.sin(theta), np.cos(theta)]
    ])

    Roated_inertia_Matrix1 = np.matmul(Rotation_Matrix,Inertia_Matrix)
    Roated_inertia_Matrix = np.matmul(Roated_inertia_Matrix1,np.transpose(Rotation_Matrix))

    Ixx = float(Roated_inertia_Matrix[0][0]) + area*(Center_coords[1]**2)
    Ixy = float(-Roated_inertia_Matrix[0][1])  + area*(Center_coords[0]*Center_coords[1])
    Iyy = float(Roated_inertia_Matrix[1][1])  + area*(Center_coords[0]**2)

    if state == 1:
        return [Ixx,Iyy,Ixy]
    elif state == 0:
        return [-Ixx,-Iyy,-Ixy]
    else:
        pass

def Moment_of_inertia_Triangle(T):
    a,b,c = T[0][0],T[0][1],T[0][2]

    Cosc = (a**2+b**2-c**2)/(2*a*b)
    Sinc = (1-Cosc**2)**0.5
    h = b*Sinc
    d = b*Cosc
    area = h*a/2
    vertices = Triangular_calculation(T)
    Cx,Cy = 0,0
    for i in vertices:
        Cx += i[0]/3
        Cy += i[1]/3
    Centroid = [Cx,Cy]
    Rotation = T[2]
    theta = np.radians(Rotation)
    state = T[3]

    Inertia_Matrix = np.array([
        [a*(h**3)/36, -a*(h**2)*(a-2*d)/72],
        [ -a*(h**2)*(a-2*d)/72,(h*(a**3)-(a**2)*h*d+(d**2)*h*a)/36]
    ])

    Rotation_Matrix = np.array([
        [np.cos(theta), np.sin(theta)],
        [-np.sin(theta), np.cos(theta)]
    ])

    Roated_inertia_Matrix1 = np.matmul(Rotation_Matrix,Inertia_Matrix)
    Roated_inertia_Matrix = np.matmul(Roated_inertia_Matrix1,np.transpose(Rotation_Matrix))

    Ixx = float(Roated_inertia_Matrix[0][0]) + area*(Centroid[1]**2)
    Ixy = float(-Roated_inertia_Matrix[0][1])  + area*(Centroid[0]*Centroid[1])
    Iyy = float(Roated_inertia_Matrix[1][1])  + area*(Centroid[0]**2)

    if state == 1:
        return [Ixx,Iyy,Ixy]
    elif state == 0:
        return [-Ixx,-Iyy,-Ixy]
    else:
        pass

def Moment_of_inertia_Circle(T):
    Radius = T[0][0]
    Initial_Rotation = T[2]
    Sector_angle = T[3]
    Rad_Sector = np.radians(Sector_angle)
    Rad_Rotation = np.radians(Initial_Rotation)
    Center_coords = [(Radius/(3*Rad_Sector))*(np.cos(Rad_Rotation)-np.cos(Rad_Sector+Rad_Rotation)),(Radius/(3*Rad_Sector))*(np.sin(Rad_Sector+Rad_Rotation)-np.sin(Rad_Rotation))]
    state =  T[4]
    area = Sector_angle*Radius*Radius

    Inertia_Matrix = np.array([
        [0.25*(Radius**4)*(np.sin(2*(Rad_Rotation+Rad_Sector))*0.25 - 0.25*np.sin(2*Rad_Rotation) + 0.5*Rad_Sector), -0.125*(Radius**4)*((np.sin(Rad_Sector+Rad_Rotation))**2 - (np.sin(Rad_Sector))**2)],
        [-0.125*(Radius**4)*((np.sin(Rad_Sector+Rad_Rotation))**2 - (np.sin(Rad_Sector))**2),
         0.25*(Radius**4)*(0.25*np.sin(2*(Rad_Rotation+Rad_Sector))- 0.25*np.sin(2*Rad_Sector)- 0.5*Rad_Sector)]
    ])

    Ixx = float(Inertia_Matrix[0][0]) + area*(Center_coords[1]**2)
    Ixy = float(-Inertia_Matrix[0][1]) + area*(Center_coords[0]*Center_coords[1])
    Iyy = float(Inertia_Matrix[1][1]) + area(Center_coords[0]**2)

    if state == 1:
        return [Ixx,Iyy,Ixy]
    elif state == 0:
        return [-Ixx,-Iyy,-Ixy]
    else:
        pass

#Main Code
Polygons_list = ["R","C","T"]
fIxx,fIyy,fIxy = 0,0,0
for i in Polygons_list:
    if i in Shapes_Dict.keys():
        processing_data = Shapes_Dict.get(i)
        for T in processing_data:
            if i == "R":
                Rectangle(T)
                val = Moment_of_inertia_Rectangle(T)
                fIxx += val[0]
                fIyy += val[1]
                fIxy += val[2]
            elif i == "C":
                Circle(T)
                val = Moment_of_inertia_Circle(T)
                fIxx += val[0]
                fIyy += val[1]
                fIxy += val[2]
            elif i == "T":
                Triangle(T)
                val = Moment_of_inertia_Triangle(T)
                fIxx += val[0]
                fIyy += val[1]
                fIxy += val[2]
            else: pass

print([fIxx,fIyy,fIxy])
ax.set_xlim(-50,50)
ax.set_ylim(-50,50)
plt.show()