import numpy as np
from atmosfera_estandar import atmosfera_estandar
from scipy.optimize import fmin, fsolve
from scipy import interpolate
from enum import IntEnum
from airfoil_characteristics import Airfoil

GRAVITY = 9.81
FEET2METER = 0.3048
METER2FEET = 1/FEET2METER
POUND2KILOGRAM = 0.453592
KG2POUND = 1 / POUND2KILOGRAM
HP2WATT = 735.499
WATT2HP = 1/HP2WATT

class Airfoil_old:
    def __init__(self, name='unidentified', a2d=None, alpha0=None, alpha_cl_max=None):
        self.name = name
        self.a2d = a2d
        self.alpha_cl_max = alpha_cl_max
        self.alpha0 = alpha0

    def load_airfoil_data(self, airfoil_name):
        pass


#TODO: def swirl_thrust_loss():


class EngineType(IntEnum):
    ALLISON250 = 1
    PRATT_WHITNEY = 2
    LYCOMING = 3

class Helicopter:

    class Rotor:
        class Blade:
            def __init__(self, airfoil=None, length=None, twist=None, chord=None, blade_tip_lost_factor=None, cd0=None):
                self.airfoil = airfoil if airfoil else Airfoil()
                self.chord = chord
                self.length = length
                self.twist = twist
                self.Fp = blade_tip_lost_factor
                self.cd0 = cd0
                self.area = self.chord * self.length

            def local_twist(self, r):
                return self.twist * r/self.length

        def __init__(self, blade=None, number_of_blades=None, inertia=None, angular_velocity=None):
            self.blade = blade if blade else self.Blade()
            self.number_of_blades = number_of_blades
            self.radius = self.blade.length
            self.diameter = self.radius*2
            self.inertia = inertia
            self.angular_velocity = angular_velocity
            self.solidity = self.number_of_blades * self.blade.chord / (self.radius * np.pi)
            self.area = np.pi * self.radius**2
            self.blade_area = self.number_of_blades * self.blade.area
            if self.angular_velocity:
                self.tip_speed = self.angular_velocity * self.radius
            else:
                self.tip_speed = None

        def a2d_corrected(self, r, height):
            _, _, _, _, _, _, vson = atmosfera_estandar('altura', height)
            local_speed = self.angular_velocity * r * self.radius
            local_mach = local_speed / vson
            return self.blade.airfoil.a2d / np.sqrt(1-local_mach**2)

        def induce_angle(self, r, theta, height):
            local_angle = self.blade.local_twist(r) - self.blade.airfoil.alpha0 + theta
            return self.a2d_corrected(r,height)*self.solidity / (16*r) *\
                   (-1 * np.sqrt(1 + 32*r*local_angle)/self.a2d_corrected(r,height)/self.solidity)

        def local_angle(self, r, theta, height):
            return self.blade.local_twist(r) -\
                   np.arctan(self.induce_angle(r, theta, height)) + self.blade.airfoil.alpha0

        def wake_rotation_factor(self, input):
            "Look up for function definition in prouty"
            input_points = np.arange(start=0, stop=0.5, step=0.005)
            output_points = np.array([0, 0.0165, 0.027, 0.0395, 0.049, 0.060,
                                      0.070, 0.080, 0.090, 0.100, 0.110])
            spline = interpolate.splrep(input_points, output_points, s=0)
            return interpolate.splev(input, spline, der=0)

        def tip_vortex_interference(self, input):
            "Look up for function definition in prouty"
            input_points = np.arange(start=0, stop=1.2, step=0.1)
            output_points = np.array([0.0, 0.945, 0.970, 0.995, 1.015,
                                      1.030, 1.040, 1.050, 1.060, 1.065,
                                      1.070, 1.073, 1.077, 1.080])
            spline = interpolate.splrep(input_points, output_points, s=0)
            return interpolate.splev(input, spline, der=0)

    class Engine:
        def __init__(self, engine_type, potencia_maxima, ):
            if isinstance(engine_type, EngineType):
                self.engine_type = engine_type
                if engine_type == EngineType.ALLISON250:
                    self.coeficiente_potencia_en_altura = 3.4e-5*FEET2METER # [%/m]
                elif engine_type == EngineType.PRATT_WHITNEY:
                    self.coeficiente_potencia_en_altura = 3.6e-5*FEET2METER # [%/m]
                elif engine_type == EngineType.LYCOMING:
                    self.coeficiente_potencia_en_altura = 3.7e-5*FEET2METER # [%/m]
                else:
                    raise RuntimeError("Data for {} engine not available".format(engine_type.name))
            else:
                raise TypeError("{} is not an avaiable engine type.")
            self.potencia_maxima = potencia_maxima

        def potencia_disponible(self, altura):
            return self.potencia_maxima * (1 - self.coeficiente_potencia_en_altura * altura)



    class FlightCondition:
        def __init__(self, velocity, height, deltaT=0):
            self.velocity = velocity
            self.height = height
            atmospheric_data = atmosfera_estandar('altura', self.height, deltaT=deltaT)
            self.density = atmospheric_data[4]
            self.temperature = atmospheric_data[3]
            self.speed_of_sound = atmospheric_data[6]

    def __init__(self, rotor=None, tail_rotor=None, engine=None, weight=None,):
        self.rotor = rotor if rotor else self.Rotor()
        self.tail_rotor = tail_rotor if tail_rotor else self.Rotor()
        self.weight = weight
        self.disk_loading = self.weight / self.rotor.area
        self.f = 0.8 * (self.weight/1000.0)**(2.0/3.0) # equivalent plate's surface
        self.flight_condition = self.FlightCondition(0.0, 0.0)
        self.engine = engine

    def coeficiente_de_traccion_pf(self, density):
        return (self.weight * GRAVITY) / (density * self.rotor.area * self.rotor.tip_speed**2)

    def coeficiente_de_traccion_pf_max(self):
        # TODO: este valor en realidad depende de la entrada en perdida de las palas del rotor.
        # a falta herramientas para calcularlo, asumimos el valor de 1.4, tomado de los graficos
        # de Ct/solidity vs Cq/solidity, para distintos valores de alabeo, solidez, y mach de
        # de divergencia, que se encuentran al final del capitulo 1 de "Helicopter performance,
        # stability, and control".
        return 0.165

    def coeficiente_de_potencia_inducida_pf(self, density):
        return np.sqrt(2.0)/2.0 * self.coeficiente_de_traccion_pf(density)**(3.0/2.0)

    def velocidad_inducida_0(self, density):
        return np.sqrt(self.weight * GRAVITY / (2 * density * self.rotor.area)) * 1/self.rotor.blade.Fp

    def lambda_inducido_0(self, density):
        return self.velocidad_inducida_0(density) / self.rotor.angular_velocity / self.rotor.radius

    def velocidad_inducida(self, velocity, density):
        vi0 = self.velocidad_inducida_0(density)
        return vi0 * np.sqrt(np.sqrt(1 + 0.25 * (velocity / vi0) ** 4) -
                             0.5 * (velocity / vi0) ** 2)
    def lambda_inducido(self, velocity, density):
        return self.velocidad_inducida(velocity, density) / self.rotor.tip_speed

    def coeficiente_de_traccion(self, velocity, density):
        vel = velocity if velocity else self.flight_condition.velocity
        den = density if density else self.flight_condition.density
        mu = vel/self.rotor.tip_speed
        return self.coeficiente_de_traccion_pf(den) * np.sqrt(1 +
                    ((self.f * mu ** 2 / 2) / (self.rotor.blade.area
                     * self.coeficiente_de_traccion_pf(den)/self.rotor.solidity ) )**2)

    def potencia_inducida(self, velocity, density):
        vi = self.velocidad_inducida(velocity, density)
        return self.weight * GRAVITY * vi

    def potencia_parasita(self, velocity, density):
        return 0.125 * density * self.rotor.number_of_blades * self.rotor.blade.chord * self.rotor.blade.cd0 * \
               self.rotor.radius * (self.rotor.angular_velocity * self.rotor.radius)**3 *\
               (1 + 4.6 * (velocity / (self.rotor.angular_velocity * self.rotor.radius))**2)

    def potencia_fuselaje(self, velocity, density):
        vi0 = self.velocidad_inducida_0(density)
        return 0.5 * density * self.f * vi0**3 * (np.sqrt(1+0.25*(velocity/vi0)**4)
                                                  + 0.5 * (velocity/vi0)**2)**1.5

    def potencia_necesaria_base(self, velocity, density):
        pi = self.potencia_inducida(velocity, density)
        p0 = self.potencia_parasita(velocity, density)
        pfus = self.potencia_fuselaje(velocity, density)
        return pi + p0 + pfus

    def potencia_motor_anticupla(self, velocity, density):
        return 0.1*self.potencia_necesaria_base(velocity, density)

    def potencia_necesaria(self, velocity, density):
        pi = self.potencia_inducida(velocity, density)
        p0 = self.potencia_parasita(velocity, density)
        pfus = self.potencia_fuselaje(velocity, density)
        pacc = self.potencia_motor_anticupla(velocity, density)
        pstall = self.potencia_stall(velocity, density)
        pcomp = self.potencia_compressibilidad(velocity, density)
        return pi + p0 + pfus + pacc + pstall + pcomp

    def coeficiente_de_potencia_inducida(self, velocity, density):
        vel = velocity if velocity else self.flight_condition.velocity
        den = density if density else self.flight_condition.density
        pot_ind = self.potencia_inducida(vel, den)
        return pot_ind / (den * self.rotor.area * self.rotor.tip_speed**3)

    def coeficiente_de_potencia_reducida_inducida(self, velocity, density):
        vel = velocity if velocity else self.flight_condition.velocity
        den = density if density else self.flight_condition.density
        return self.coeficiente_de_potencia_inducida(vel, den) / \
            self.coeficiente_de_potencia_inducida_pf(den)

    def coeficiente_de_potencia_parasita(self, velocity, density):
        vel = velocity if velocity else self.flight_condition.velocity
        den = density if density else self.flight_condition.density
        pot_0 = self.potencia_parasita(vel, den)
        return pot_0 / (den * self.rotor.area * self.rotor.tip_speed ** 3)

    def coeficiente_de_potencia_reducida_parasita(self, velocity=None, density=None):
        vel = velocity if velocity else self.flight_condition.velocity
        den = density if density else self.flight_condition.density
        return self.coeficiente_de_potencia_parasita(vel, den) / \
            self.coeficiente_de_potencia_inducida_pf(den)

    def coeficieciente_de_potencia_fuselaje(self, velocity, density):
        vel = velocity if velocity else self.flight_condition.velocity
        den = density if density else self.flight_condition.density
        pot_0 = self.potencia_parasita(vel, den)
        return pot_0 / (den * self.rotor.area * self.rotor.tip_speed ** 3)

    def coeficiente_de_potencia_reducida_fuselaje(self, velocity=None, density=None):
        vel = velocity if velocity else self.flight_condition.velocity
        den = density if density else self.flight_condition.density
        return self.coeficiente_de_potencia_parasita(vel, den) / \
               self.coeficiente_de_potencia_inducida_pf(den)

    def coeficieciente_de_potencia_necesaria_base(self, velocity, density):
        vel = velocity if velocity else self.flight_condition.velocity
        den = density if density else self.flight_condition.density
        pot_nec = self.potencia_necesaria_base(vel, den)
        return pot_nec / (den * self.rotor.area * self.rotor.tip_speed ** 3)

    def coeficieciente_de_potencia_necesaria(self, velocity, density):
        vel = velocity if velocity else self.flight_condition.velocity
        den = density if density else self.flight_condition.density
        pot_nec = self.potencia_necesaria(vel, den)
        return pot_nec / (den * self.rotor.area * self.rotor.tip_speed ** 3)

    def coeficiente_de_potencia_reducida_necesaria_base(self, velocity=None, density=None):
        vel = velocity if velocity else self.flight_condition.velocity
        den = density if density else self.flight_condition.density
        return self.coeficieciente_de_potencia_necesaria_base(vel, den) / \
               self.coeficiente_de_potencia_inducida_pf(den)

    def coeficiente_de_potencia_reducida_necesaria(self, velocity=None, density=None):
        vel = velocity if velocity else self.flight_condition.velocity
        den = density if density else self.flight_condition.density
        return self.coeficieciente_de_potencia_necesaria(vel, den) / \
               self.coeficiente_de_potencia_inducida_pf(den)

    def velocidad_de_potencia_minima(self, density, vel_guess=10):
        func = lambda vel: self.potencia_necesaria(vel, density)
        rta = fmin(func, vel_guess, disp=False)
        return rta[0]

    def velocidad_de_maximo_alcance(self, density, vel_guess=10):
        den = density if density else self.flight_condition.density
        func = lambda vel: self.potencia_necesaria(vel, den) / vel
        xopt = fmin(func, vel_guess, disp=False)
        return xopt

    def velocidad_autorotacion(self, velocity, density):
        return self.potencia_necesaria(velocity, density) / (self.weight * GRAVITY)

    def tiempo_de_punto_fijo_equivalente(self, altura):
        """ Extraido de "Helicopter Performance, Stability and Control",
         Capitulo 4, pag. 363, Prouty. A su vez, el libro hace referencia a
         "High Energy Rotor System", Wood, AHS 32nd Forum, 1976"""
        _, _, _, _, density, _, _ = atmosfera_estandar('altura', altura)
        return self.rotor.inertia * self.rotor.angular_velocity**2 * \
        (1 - self.coeficiente_de_traccion_pf(density=density) /
         (0.8*self.coeficiente_de_traccion_pf_max()))\
        / (1.1e3*self.engine.potencia_disponible(altura) * WATT2HP)

    def indice_de_autorotacion(self, density):
        """ Indice de autorotacion en ft**3/lb.
            Extraido de "Helicopter Performance, Stability and Control",
            Capitulo 4, pag. 363, Prouty. A su vez, el libro hace referencia a
            "A Simple Autorotative Flare Index", Fradenburgh, JAHS 29-3, 1984 """
        _, _, _, _, density0, _, _ = atmosfera_estandar('altura', 0)
        return (self.rotor.inertia * self.rotor.angular_velocity**2 / self.weight) \
                * (density / density0 / self.disk_loading) * ( METER2FEET**3 / POUND2KILOGRAM)


    def velocidad_de_autorotacion_minima(self, density, vel_guess=10):
        func = lambda vel: self.velocidad_autorotacion(vel, density)
        xopt, fopt = fmin(func, vel_guess)
        return (xopt, fopt)

    def potencia_excedente(self, velocity, altura):
        vel = velocity if velocity else self.flight_condition.velocity
        h = altura if altura else self.flight_condition.density
        _, _, _, _, den, _, _ = atmosfera_estandar('altura', h)
        pot_nec = self.potencia_necesaria(vel, den)
        pot_disp = self.engine.potencia_disponible(h)
        return pot_disp - pot_nec

    def traccion_necesaria(self, velocity, density):
        resistencia = 0.5*velocity**2*density*self.f
        return np.sqrt((self.weight * GRAVITY)**2+resistencia**2)

    def velocidad_de_ascenso(self, velocity, altura):
        vel = velocity if velocity else self.flight_condition.velocity
        h = altura if altura else self.flight_condition.density
        _, _, _, _, den, _, _ = atmosfera_estandar('altura', h)
        deltaP = self.potencia_excedente(vel, h)
        traccion = self.traccion_necesaria(vel, den)
        vel_ind_0 = self.velocidad_inducida_0(den)
        return deltaP/traccion * (2*vel_ind_0 + deltaP/traccion)\
               / (vel_ind_0 + deltaP/traccion)

    def techo_practico(self, velocity, velocidad_umbral=0.5, h_guess=15000):
        vel = velocity if velocity else self.flight_condition.velocity
        func = lambda altura: self.velocidad_de_ascenso(vel, altura) - velocidad_umbral
        techo_practico = fsolve(func, x0=h_guess, xtol=1e-12)
        return techo_practico

    def alpha_PTP(self, velocity=None, density=None):
        den = density if density else self.flight_condition.density
        vel = velocity if velocity else self.flight_condition.velocity
        mu = vel/self.rotor.tip_speed
        return - self.f * mu**2 / (2*self.rotor.blade.area *
            (self.coeficiente_de_traccion_pf(den)/self.rotor.solidity))

    def lambda_pasante(self, velocity=None, density=None):
        den = density if density else self.flight_condition.density
        vel = velocity if velocity else self.flight_condition.velocity
        mu = vel / self.rotor.tip_speed
        return mu * self.alpha_PTP(vel, den) - self.lambda_inducido(vel, den)

    def theta_0(self, velocity=None, density=None):
        den = density if density else self.flight_condition.density
        vel = velocity if velocity else self.flight_condition.velocity
        mu = vel/self.rotor.tip_speed
        return (4/self.rotor.blade.airfoil.a2d * (1 + 3/2*mu**2 - 5/24*mu**4)
                * self.coeficiente_de_traccion(velocity, den)/ self.rotor.solidity
                - (1 - 3/2*mu**2 + 1/4*mu**4)*self.rotor.blade.twist/2
                - (1 + 13/24*mu**4) * self.lambda_pasante(velocity,den))\
            / (3/2 - 3/2*mu**2 - 8/9/np.pi*mu**3 + 25/36*mu**4)

        # return ((12 / self.rotor.blade.airfoil.a2d * (self.rotor.blade.Fp ** 2 / 3 + mu ** 2 / 2)\
        #     * self.coeficiente_de_traccion(vel, density) / self.rotor.solidity
        #     - self.rotor.blade.Fp ** 2 / 2 * (self.rotor.blade.Fp ** 4
        #     + 3 / 2 * mu ** 2 * (mu ** 2 - self.rotor.blade.Fp ** 2)) * self.rotor.blade.twist
        #     - (self.rotor.blade.Fp ** 4 - mu ** 2 / 2) * self.lambda_pasante(vel, den))
        #     / (self.rotor.blade.Fp * (2 / 3 * self.rotor.blade.Fp ** 4
        #     - 3 / 2 * mu ** 2 * self.rotor.blade.Fp ** 2 + 3 / 2 * mu ** 4)))

    def theta_34(self, velocity=None, density=None):
        den = density if density else self.flight_condition.density
        vel = velocity if velocity else self.flight_condition.velocity
        return self.theta_0(vel, den) + 0.75*self.rotor.blade.twist

    def aleteo_long(self, velocity=None, density=None):
        den = density if density else self.flight_condition.density
        vel = velocity if velocity else self.flight_condition.velocity
        mu = vel/self.rotor.tip_speed
        return (mu * (8 / 3 * self.rotor.blade.Fp * self.theta_0(vel, den) +
                      2 * self.rotor.blade.Fp ** 2 * self.rotor.blade.twist +
                      2 * self.lambda_pasante(vel, den)) /
                (self.rotor.blade.Fp ** 2 + 3 / 2 * mu ** 2))

    def alpha_tip_90(self, radius, velocity=None, density=None):
        """" No se usa por el momento, pero se podria usar para calcular
        el radio de entrada en perdida"""
        den = density if density else self.flight_condition.density
        vel = velocity if velocity else self.flight_condition.velocity
        mu = vel/self.rotor.tip_speed
        return (1 / (radius + mu) * (self.theta_0(vel, den) +
                                     self.rotor.blade.twist*radius -
                                    self.aleteo_long(vel, den) -
                                    self.lambda_inducido(vel, den) +
                                    mu * (self.alpha_PTP(vel, den) +
                                      self.theta_0(vel, den) +
                                      self.rotor.blade.twist * radius -
                                      self.aleteo_long(vel, den))))

    def B(self, velocity, density):
        mu = velocity / self.rotor.angular_velocity / self.rotor.radius
        return self.theta_0(velocity, density) - self.rotor.blade.airfoil.alpha_cl_max \
               + self.aleteo_long(velocity, density) - mu * self.rotor.blade.twist

    def C(self, velocity, density):
        mu = velocity / self.rotor.tip_speed
        return mu*(self.alpha_PTP(velocity, density)+self.rotor.blade.airfoil.alpha_cl_max)\
                - self.lambda_inducido(velocity, density) - mu*self.theta_0(velocity, density)\
                - mu*self.aleteo_long(velocity, density)

    def stall_relative_radius(self, velocity, density):
        discriminant = self.B(velocity, density)**2 -\
                        4*self.rotor.blade.twist*self.C(velocity, density)
        if discriminant > 0:
            return (np.sqrt(discriminant) - self.B(velocity, density))\
                   / (2*self.rotor.blade.twist)
        else:
            return 0

    def kp(self, velocity, density):
        discriminante = self.B(velocity, density) / 2.0 / self.rotor.blade.twist
        if -discriminante <= 1:
            return -(discriminante + self.stall_relative_radius(velocity, density))\
                     / (1 - self.stall_relative_radius(velocity, density))
        else:
            return 1

    def delta_cp_stall(self, velocity, density):
        mu = velocity / self.rotor.angular_velocity / self.rotor.radius
        dcp =  self.kp(velocity, density) * self.rotor.solidity / 24.0 / np.pi \
               * (1-mu)**2 \
               * (1 - self.stall_relative_radius(velocity, density))\
               * np.sqrt(1 - self.stall_relative_radius(velocity, density)**2)
        if dcp > 0:
            return dcp
        else:
            return 0

    def potencia_stall(self, velocity, density):
        return self.delta_cp_stall(velocity, density)\
               * (density * self.rotor.area * self.rotor.tip_speed ** 3)

    def mach_tip(self, altura):
        _, _, _, _, _, _, vson = atmosfera_estandar('altura', altura)
        return self.rotor.tip_speed / vson

    def mach_critico_simple(self, alpha):
        """" Ecuacion lineal simple (falta fuente)"""
        return 0.7264957 - 2.44 * alpha

    def delta_mach_divergence(self, velocity, altura):
        _, _, _, _, density, _, _ = atmosfera_estandar('altura', altura)
        mu = velocity / self.rotor.tip_speed
        mach_tip = self.mach_tip(altura)
        alpha = self.alpha_tip_90(radius=1, velocity=velocity, density=density)
        return mach_tip * (1 + mu) - self.mach_critico_simple(alpha) - 0.06

    def delta_cp_mach(self, velocity, density):
        dcp = self.rotor.solidity * (0.012 * self.delta_mach_divergence(velocity, density)
                                + 0.1 * self.delta_mach_divergence(velocity, density)**3)
        if dcp > 0:
            return dcp
        else:
            return 0

    def potencia_compressibilidad(self, velocity, density):
        return self.delta_cp_mach(velocity, density) \
               * (density * self.rotor.area * self.rotor.tip_speed ** 3)
