import math
import matplotlib.pyplot as plt
import numpy as np

### Constantes

g = 9.806         # gravitation [m/s**2]

### Paramètres du système

l = 1.15         # longueur du pendule [m]
m = 0.1          # masse du prehenseur [kg]
b = 0.0002         # coefficient de frottement [kg*m^2/s]


### Paramètres de la simulation

step = 0.001                # pas (dt) [s]
end = 20.0                  # durée [s]
theta_0 = math.pi / 6       # position angulaire initiale [rad]
theta_dot_0 = 0             # vitesse_angulaire initiale [rad/s]
x_c_0 = -5.0                # position du cart initiale [m]
v_c_0 = 0.0                 # vitesse du cart initiale initiale [m/s]

### Fonction de déplacement du Cart

def get_a_c(t):
    w = 1 # pulsation
    return const(0)				# constant acceleration
    #return sinusoidal(t, w)		# sinusoidal
    #return square_pattern(t,w)	# 1 or -1
    return triangle(t,w/math.pi)		# triangle pattern
    return max(-t + 10, 0)      # gradually goes from 10 to 0 [m/s^2]




def const(x):
    return x

def sinusoidal(x, w):
    return math.sin(w * x)

def square_pattern(x,w):
    s = math.sin(w*x)
    if s != 0:
        return abs(s)/s
    return 0

def triangle(x,w):
    y = w*x
    return abs(y-int(y)-0.5) * 2 - 0.5

t = np.arange(0, end, step)
theta = np.empty_like(t)          
theta_dot = np.empty_like(t)
theta_dot_dot = np.empty_like(t)

x_c = np.empty_like(t)
v_c = np.empty_like(t)
a_c = np.empty_like(t)
    
def simulation():
    """
    pre: -
    post: exécute une simulation jusqu'à t=end par pas de dt=step.
          Remplit les listes theta, theta_dot, theta_dot_dot
          avec les positions, vitesses et accélérations angulaire du pendule.
    """
    # conditions initiales
    theta[0] = theta_0
    theta_dot[0] = theta_dot_0
    
    x_c[0] = x_c_0
    v_c[0] = v_c_0
    
    for i in range(len(t)-1):

        dt = step
        
        # calcule de la position du cart
        a_c[i] = get_a_c(t[i])
        v_c[i+1] = v_c[i] + a_c[i] * dt
        x_c[i+1] = x_c[i] + v_c[i] * dt
        
        # calcule de l'acceleration avec les conditions actuelles
        theta_dot_dot[i] = -g/l*math.sin(theta[i]) - a_c[i]/l*math.cos(theta[i]) - b/(m*l**2)*theta_dot[i]
        # calcul accélération, vitesse, position
        theta_dot[i+1] = theta_dot[i] + theta_dot_dot[i] * dt
        theta[i+1] = theta[i] + theta_dot[i] * dt


def graphiques():
    
    plt.figure(1)
    plt.subplot(3,1,1)
    plt.plot(t,theta, label="angular position")
    plt.legend()
    plt.subplot(3,1,2)
    plt.plot(t,theta_dot, label="angular speed")
    plt.legend()
    plt.subplot(3,1,3)
    plt.plot(t,theta_dot_dot, label="angular acceleration")
    plt.legend()
    plt.show()

def get_kinetic_energy(i):
    """ Calculates kinetic energy
    pre: i: int, current iteration
    post: float, kinetic energy
    """
    
    vel_c = np.array([v_c[i], 0])       # vitesse du chariot vecteur en 2 dimensions (x,y) -> (v_c, 0)
    
    speed = l * theta_dot[i]            # vitesse du pendule relative au chariot (scalaire)
    vel_rel_p = np.array([speed * math.sin(theta[i]), speed * math.cos(theta[i])])  # vitesse vectorielle
    
    vel_abs = vel_c + vel_rel_p         # vitesse absolue du pendule
    
    speed_abs = math.sqrt(vel_abs[0]**2 + vel_abs[1]**2)
    
    return 0.5 * m * speed_abs**2
    
def graphiques_energie():
    e_pot = np.empty_like(t)
    e_kin = np.empty_like(t)
    e_tot = np.empty_like(t)
    for i in range(len(t)):
        e_pot[i] = m * g * l * (1-math.cos(theta[i]))   # mgh
        e_kin[i] = get_kinetic_energy(i) # speed of cart + speedc
        e_tot[i] = e_pot[i] + e_kin[i]

    print(e_kin[len(e_kin)-1])
    plt.figure(2)
    plt.plot(t, e_pot , label="potential")
    plt.plot(t, e_kin, label="kinetic")
    plt.plot(t, e_tot, label="total")
    plt.legend()
    plt.show()

### programme principal

simulation()
graphiques()
graphiques_energie()

