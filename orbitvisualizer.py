import numpy as np
import matplotlib.pyplot as plt
import datetime
import math

plt.style.use('dark_background')
ax = plt.figure().add_subplot(projection='3d')
cos = np.cos
sin = np.sin
pi = np.pi
step = 3600
theta = np.linspace(0,2*pi, step)

# Possible improvements?
# Allow for input of any datetime

class orbits:

    def __init__(planet, a, e, inc, node, arg, col, name):

        planet.a = a
        planet.e = e
        planet.inc = inc
        planet.node = node
        planet.arg = arg
        planet.col = col
        planet.name = name



    def orbvis(planet, plot=True):

        inc = np.radians(planet.inc)
        node = np.radians(planet.node)
        arg = np.radians(planet.arg)

        r = []
        x1 = []
        y1 = []

        if planet.e < 1:
            for j in range(len(theta)):
                r.append((planet.a*(1-planet.e**2))/(1+planet.e*cos(theta[j])))
                x1.append(r[j] * cos(theta[j]))
                y1.append(r[j] * sin(theta[j]))
        elif planet.e > 1:
            angle = np.linspace(-(pi/2), (pi/2), step)
            for j in range(len(angle)):
                r.append((planet.a*(1-planet.e**2))/(1+planet.e*cos(angle[j])))
                x1.append(r[j] * cos(90 + angle[j]))
                y1.append(r[j] * sin(90 + angle[j]))



        # implement inclination
        z1 = []
        for j in range(len(r)):
            z1.append(y1[j]*np.tan(inc))


        # implement longitude of Ascending Node
        x2 = []
        y2 = []
        for j in range(len(x1)):
            x2.append((x1[j]*(cos(node)) - y1[j]*(sin(node))))
            y2.append((x1[j]*(sin(node)) + y1[j]*(cos(node))))


        # implement argument of periapsis

        p1 = np.array([0,0,0])
        p2 = np.array([x2[50],y2[50],z1[50]])
        p3 = np.array([x2[100],y2[100],z1[100]])
        v1 = p2 - p1
        v2 = p3 - p1
        uv = np.cross(v1,v2)
        um = np.linalg.norm(uv)
        u = uv/um

        R = [
            [(u[0]**2)*(1-cos(arg)) + cos(arg), u[0]*u[1]*(1-cos(arg))-u[2]*sin(arg), u[0]*u[2]*(1-cos(arg))+u[1]*sin(arg)],
            [u[0]*u[1]*(1-cos(arg)) + u[2]*sin(arg), (u[1]**2)*(1-cos(arg))+cos(arg), u[1]*u[2]*(1-cos(arg))-u[0]*sin(arg)],
            [u[0]*u[2]*(1-cos(arg)) - u[1]*sin(arg), u[1]*u[2]*(1-cos(arg))+u[0]*sin(arg), (u[2]**2)*(1-cos(arg))+cos(arg)]
        ]
        planet.xs = []
        planet.ys = []
        planet.zs = []
        for j in range(len(r)):
            planet.xs.append((x2[j]*R[0][0]) + (y2[j]*R[0][1]) + (z1[j]*R[0][2]))
            planet.ys.append((x2[j]*R[1][0]) + (y2[j]*R[1][1]) + (z1[j]*R[1][2]))
            planet.zs.append((x2[j]*R[2][0]) + (y2[j]*R[2][1]) + (z1[j]*R[2][2]))


        # Vector pointing toward periapsis

        #per_i = np.argmin(r)
        #ax.plot((0,planet.xs[per_i]),(0,planet.ys[per_i]),(0,planet.zs[per_i]), color='pink')
        if plot == True:
          ax.plot(planet.xs,planet.ys,planet.zs,color=planet.col, label=planet.name)





    def realtimeplot(planet, Tyears, JD, ToP): # orbital period in years, Julian Date of desired time, Time of (preferrably most recent) Periapsis
        diff = abs(JD - ToP)
        Tdays = Tyears * 365.25
        percent = diff/Tdays

        if percent > 1:
            whole = int(percent)
            percent = percent - whole
            pos = percent * 3600
            pos = int(pos)
            ax.plot(planet.xs[pos],planet.ys[pos],planet.zs[pos],color=planet.col,marker='.')
        elif (percent == 1) or (percent == 0):
            ax.plot(planet.xs[0],planet.ys[0],planet.zs[0],color=planet.col,marker='.')
        elif percent < 1:
            pos = percent * 3600
            pos = int(pos)
            ax.plot(planet.xs[pos],planet.ys[pos],planet.zs[pos],color=planet.col,marker='.')
        else:
            pos = percent * 3600
            pos = int(pos)
            ax.plot(planet.xs[pos],planet.ys[pos],planet.zs[pos],color=planet.col,marker='.')
      
    def realtimevector(planet, Tyears, JD, ToP): # orbital period in years, Julian Date of desired time, Time of (preferrably most recent) Periapsis
        planet.orbvis(plot=False)
        diff = abs(JD - ToP)
        Tdays = Tyears * 365.25
        percent = diff/Tdays

        if percent > 1:
            whole = int(percent)
            percent = percent - whole
            pos = percent * 3600
            pos = int(pos)
            vector = [planet.xs[pos],planet.ys[pos],planet.zs[pos]]
            return vector
        elif (percent == 1) or (percent == 0):
            vector = [planet.xs[pos],planet.ys[pos],planet.zs[pos]]
            return vector
        elif percent < 1:
            pos = percent * 3600
            pos = int(pos)
            vector = [planet.xs[pos],planet.ys[pos],planet.zs[pos]]
            return vector
        else:
            pos = percent * 3600
            pos = int(pos)
            vector = [planet.xs[pos],planet.ys[pos],planet.zs[pos]]
            return vector

def jdn(dto): #Datetime object input

    year = dto.year
    month = dto.month
    day = dto.day

    not_march = month < 3
    if not_march:
        year -= 1
        month += 12

    fr_y = math.floor(year / 100.0000)
    reform = 2 - fr_y + math.floor(fr_y / 4.0000)
    jjs = day + (
        math.floor(365.25 * (year + 4716.0000)) + math.floor(30.6001 * (month + 1)) + reform - 1524.0000)

    day_fraction = dto.hour / 24.0 + dto.minute / 1440.0 + dto.second / 86400.0
    return jjs + day_fraction - 0.5

# ===========================================
# Add orbits here please
# ===========================================

# Earth's parameters as an example
Earth_a = 1.000001 # AU
Earth_e = 0.0167
Earth_inc = 7.155 # degrees
Earth_node = 348.74 # degrees
Earth_arg = 114.207 # degrees


now = jdn(datetime.datetime.now(datetime.timezone.utc))
print(now)
earth = orbits(Earth_a, Earth_e, Earth_inc, Earth_node, Earth_arg, 'blue', 'Earth')
earth.orbvis()
earth.realtimevector(1, now, 2461044.5)
mars = orbits(1.523680, 0.0934, 5.65, 49.578, 286.5, 'chocolate', 'Mars')
mars.orbvis()
venus = orbits(0.723332, 0.00677, 3.86, 76.68, 54.884, 'goldenrod', 'Venus')
venus.orbvis()
atlas = orbits(-0.26, 6.13, 175.11, 322, 128, 'red', 'ATLAS')
atlas.orbvis()







# =============================================
# Fixed Plot Stuff
# =============================================
#'''
X_arrow = [0,0.2]
Y_arrow = [0,0.2]
Z_arrow = [0,0.2]

ax.plot(X_arrow,(0,0),(0,0),color='red')
ax.plot((0,0),Y_arrow,(0,0),color='magenta')
ax.plot((0,0),(0,0),Z_arrow,color='magenta')
#'''
ax.set_xlabel('X Distance (AU)')
ax.set_ylabel('Y Distance (AU)')
ax.set_zlabel('Z Distance (AU)')
ax.plot(0,0,0,color='gold',marker='o',label='Sun')
ax.set_xlim(-2, 2)
ax.set_ylim(-2,2)
ax.set_zlim(-1, 1)
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.xaxis.pane.set_edgecolor('w')
ax.yaxis.pane.set_edgecolor('w')
ax.zaxis.pane.set_edgecolor('w')
ax.grid(visible=False)
ax.legend()
plt.show()
