import math
import matplotlib.pyplot as plt
import numpy as np

### Constantes
g = 9.81                        # gravitation [m/s**2]

### Paramètres du système
m = 0.1                         # masse du bloc [kg]
l = 0.15                        # longeur du cable [m]
b = 0.001                      #coefficient de frottement [kg*m**2/s]

### Paramètres de la simulation
step = 0.001                    # pas (dt) [s]
end = 20.0                      # durée [s]
theta0 = math.pi/6              #angle initial du pendule [rad]
omega0 = 0                      #vitesse angulaire initiale du pendule [rad/s]
v_c_0 = 0                       #vitesse initiale du chariot [m/s]

t = np.arange(0, end, step)
actual_t = 0
omega = np.empty_like(t)
theta = np.empty_like(t)
theta_deg = np.empty_like(t)
v_c = np.empty_like(t)
a_c = np.empty_like(t)

### Fonction de l'acceleration du chariot
def acceleration_chariot(t):
    w = 1                        #pulsation
    #return math.sin(w * t)      #acceleration du chariot sinusoidale
    #return triangle(t,w)        #acceleration du chariot triangulaire
    return 0.0               #acceleration du chariot constante

def triangle(x,w):
    y = w*x
    return abs(y-int(y)-0.5) * 2 - 0.5

def simulation():
    dt=step
    ### Conditions initiales
    omega[0] = omega0
    theta[0] = theta0
    v_c[0] = v_c_0

    for i in range(len(omega)-1):
        # Calcul de l'acceleration et de la vitesse du chariot
        a_c[i] = acceleration_chariot(t[i])
        v_c[i+1] = v_c[i] + a_c[i] * dt
        # Calcul de la vitesse et de l'angle du pendule
        omega[i+1]= omega[i] - ((g/l) *math.sin(theta[i]) + (b/(m *l**2)) * omega[i] + (a_c[i]/l)*math.cos(theta[i]))*dt
        theta[i+1] = theta[i] + omega[i+1]*dt

    plt.figure(1)
    plt.subplot(2,1,1)
    plt.plot(t, theta, label="theta")
    plt.legend()
    plt.subplot(2,1,2)
    plt.plot(t, omega, label="omega")
    plt.legend()
    plt.show()

def energie_cinetique(i):
    vel_c = np.array([v_c[i], 0])
    
    v_l = l * omega[i]                                                              #Vitesse lineaire du pendule relative au chariot
    vel_rel_p = np.array([v_l * math.sin(theta[i]), v_l * math.cos(theta[i])])      #Vitesse vectorielle
    vel_abs = vel_c + vel_rel_p                                                     #Vitesse vectorielle absolue du pendule
    speed_abs = math.sqrt(vel_abs[0]**2 + vel_abs[1]**2)                            #Vitesse absolue du pendule (norme)
    
    return 0.5 * m * speed_abs**2

def graphiques_energie():
    e_pot = np.empty_like(t)
    e_cin = np.empty_like(t)
    e_tot = np.empty_like(t)
    for i in range(len(t)):
        e_pot[i] = m * g * l * (1-math.cos(theta[i]))
        e_cin[i] = energie_cinetique(i)
        e_tot[i] = e_pot[i] + e_cin[i]

    plt.figure(2)
    plt.plot(t, e_pot, label="potentiel")
    plt.plot(t, e_cin, label="cinétique")
    plt.plot(t, e_tot, label="total")
    plt.legend()
    plt.show()

### Programme Principal
simulation()
graphiques_energie()