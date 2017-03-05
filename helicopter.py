import numpy as np
from atmosfera_estandar import atmosfera_estandar
from scipy.optimize import fmin, fsolve
from enum import IntEnum

GRAVITY = 9.81
FEET2METER = 0.3048
METER2FEET = 1/FEET2METER

class Airfoil:
    def __init__(self, name, a2d=None, alpha_perdida=None):
        self.name = name
        self.a2d = a2d
        self.alpha_perdida = alpha_perdida

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
            self.tip_speed = self.angular_velocity * self.radius

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

    def __init__(self, rotor=None, engine=None, weight=None,):
        self.rotor = rotor if rotor else self.Rotor()
        self.weight = weight
        self.f = 0.8 * (self.weight/1000.0)**(2.0/3.0) # equivalent plate's surface
        self.flight_condition = self.FlightCondition(0.0, 0.0)
        self.engine = engine

    def coeficiente_de_traccion_pf(self, density):
        return (self.weight * GRAVITY) / (density * self.rotor.area * self.rotor.tip_speed**2)

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
                    (self.f * mu ** 2 / 2 / self.rotor.blade.area
                     / self.coeficiente_de_traccion_pf(den)/self.rotor.solidity))

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

    def potencia_necesaria(self, velocity, density, coeff_rotor_anticupla=1.1):
        pi = self.potencia_inducida(velocity, density)
        p0 = self.potencia_parasita(velocity, density)
        pfus = self.potencia_fuselaje(velocity, density)
        return coeff_rotor_anticupla*(pi + p0 + pfus)

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

    def coeficieciente_de_potencia_necesaria(self, velocity, density):
        vel = velocity if velocity else self.flight_condition.velocity
        den = density if density else self.flight_condition.density
        pot_nec = self.potencia_necesaria(vel, den)
        return pot_nec / (den * self.rotor.area * self.rotor.tip_speed ** 3)

    def coeficiente_de_potencia_reducida_necesaria(self, velocity=None, density=None):
        vel = velocity if velocity else self.flight_condition.velocity
        den = density if density else self.flight_condition.density
        return self.coeficieciente_de_potencia_necesaria(vel, den) / \
               self.coeficiente_de_potencia_inducida_pf(den)

    def velocidad_de_potencia_minima(self, density, vel_guess=10):
        func = lambda vel: self.potencia_necesaria(vel, density)
        rta = fmin(func, vel_guess)
        return rta[0]

    def velocidad_de_maximo_alcance(self, density, vel_guess=10):
        den = density if density else self.flight_condition.density
        func = lambda vel: self.potencia_necesaria(vel, den) / vel
        xopt, fopt = fmin(func, vel_guess)
        return (xopt, fopt)

    def velocidad_autorotacion(self, velocity, density):
        return self.potencia_necesaria(velocity, density) / (self.weight * GRAVITY)


    def velocidad_de_autorotacion_minima(self, density, vel_guess=10):
        func = lambda vel: self.velocidad_autorotacion(vel, density)
        xopt, fopt = fmin(func, vel_guess)
        return (xopt, fopt)

    def potencia_excedente(self, velocity, altura):
        vel = velocity if velocity else self.flight_condition.velocity
        h = altura if altura else self.flight_condition.density
        den = atmosfera_estandar('altura',h)
        pot_nec = self.potencia_necesaria(vel,den)
        pot_disp = self.engine.potencia_disponible(h)
        return pot_disp - pot_nec

    def traccion_necesaria(self, velocity, density):
        resistencia = 0.5*velocity**2*density*self.f
        return np.sqrt(self.weight**2+resistencia**2)

    def velocidad_de_ascenso(self, velocity, altura):
        vel = velocity if velocity else self.flight_condition.velocity
        h = altura if altura else self.flight_condition.density
        den = atmosfera_estandar('altura',h)
        deltaP = self.potencia_excedente(vel, den)
        traccion = self.traccion_necesaria(vel, den)
        vel_ind_0 = self.velocidad_inducida_0(den)
        return deltaP/traccion * (2*vel_ind_0 + deltaP/traccion)\
               /(vel_ind_0 + deltaP/traccion)

    def techo_practico(self, velocity, velocidad_umbral=0.5):
        vel = velocity if velocity else self.flight_condition.velocity
        func = lambda altura: self.velocidad_de_ascenso(vel, altura) - velocidad_umbral
        techo_practico, _ = fsolve(func, 5000, xtol=1e-12)
        return techo_practico

    def alpha_PTP(self, velocity=None, density=None):
        den = density if density else self.flight_condition.density
        vel = velocity if velocity else self.flight_condition.velocity
        mu = vel/self.rotor.tip_speed
        return self.f*mu**2 / (2*self.rotor.blade.area*self.coeficiente_de_traccion(vel, den))

    def lambda_pasante(self, velocity=None, density=None):
        den = density if density else self.flight_condition.density
        vel = velocity if velocity else self.flight_condition.velocity
        mu = vel / self.rotor.tip_speed
        return mu * self.alpha_PTP(vel, den) - self.lambda_inducido(vel, den)

    def theta_0(self, velocity=None, density=None):
        den = density if density else self.flight_condition.density
        vel = velocity if velocity else self.flight_condition.velocity
        mu = vel/self.rotor.tip_speed
        return ((12 / self.rotor.blade.airfoil.a2d * (self.rotor.blade.Fp ** 2 / 3 + mu ** 2 / 2)\
            * self.coeficiente_de_traccion(vel, density) / self.rotor.solidity
            - self.rotor.blade.Fp ** 2 / 2 * (self.rotor.blade.Fp ** 4
            + 3 / 2 * mu ** 2 * (mu ** 2 - self.rotor.blade.Fp ** 2)) * self.rotor.blade.twist
            - (self.rotor.blade.Fp ** 4 - mu ** 2 / 2) * self.lambda_pasante(vel, den))
            / (self.rotor.blade.Fp * (2 / 3 * self.rotor.blade.Fp ** 4
            - 3 / 2 * mu ** 2 * self.rotor.blade.Fp ** 2 + 3 / 2 * mu ** 4)))

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

    def alpha_tip_90(self, velocity=None, density=None):
        den = density if density else self.flight_condition.density
        vel = velocity if velocity else self.flight_condition.velocity
        mu = vel/self.rotor.tip_speed
        return (1 / (1 + mu) * (self.theta_0(vel, den) + self.rotor.blade.twist -
                                self.aleteo_long(vel, den) -
                                self.lambda_inducido(vel, den) +
                                mu * (self.alpha_PTP(vel, den) +
                                      self.theta_0(vel, den) +
                                      self.rotor.blade.twist -
                                      self.aleteo_long(vel, den))))



