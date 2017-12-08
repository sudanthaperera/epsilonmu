import math
import numpy as np


class MaterialType:
    def __init__(self, eps_r=1, mu_r=1, sigma_e=0, sigma_m=0):
        self.eps_r = eps_r
        self.mu_r = mu_r
        self.sigma_e = sigma_e
        self.sigma_m = sigma_m


class Brick:
    def __init__(self, min_x, min_y, min_z, max_x, max_y, max_z):
        self.min_x = min_x
        self.min_y = min_y
        self.min_z = min_z
        self.max_x = max_x
        self.max_y = max_y
        self.max_z = max_z


class VoltageSource:
    def __init__(self, min_x, min_y, min_z, max_x, max_y, max_z, direction, resistance=50, magnitude=1,
                 waveform_type='gaussian'):
        self.min_x = min_x
        self.min_y = min_y
        self.min_z = min_z
        self.max_x = max_x
        self.max_y = max_y
        self.max_z = max_z
        self.direction = direction
        self.resistance = resistance
        self.magnitude = magnitude
        self.waveform_type = waveform_type

    def resistance(self, resistance=50):
        self.resistance = resistance

    def magnitude(self, magnitude=1):
        self.magnitude = magnitude

    def waveform(self, waveform_type='gaussian'):
        self.waveform_type = waveform_type

    def get_is(self):
        return self._is

    def get_js(self):
        return self._js

    def get_ks(self):
        return self._ks

    def get_ie(self):
        return self._ie

    def get_je(self):
        return self._je

    def get_ke(self):
        return self._ke

    def initWaveform(self,fdtd_domain,yee_cell,waveform):
        self._is = round((self.min_x - fdtd_domain.min_x) / yee_cell['dx'])
        self._js = round((self.min_y - fdtd_domain.min_y) / yee_cell['dy'])
        self._ks = round((self.min_z - fdtd_domain.min_z) / yee_cell['dz'])
        self._ie = round((self.max_x - fdtd_domain.min_x) / yee_cell['dx'])
        self._je = round((self.max_y - fdtd_domain.min_y) / yee_cell['dy'])
        self._ke = round((self.max_z - fdtd_domain.min_z) / yee_cell['dz'])

        if self.direction[0] == 'x':
            self.n_fields = self._ie - self._is
            self.r_magnitude_factor = (1 + self._je - self._js) * (1 + self._ke - self._ks) / (self._ie - self._is)
        elif self.direction[0] == 'y':
            self.n_fields = self._je - self._js
            self.r_magnitude_factor = (1 + self._ie - self._is) * (1 + self._ke - self._ks) / (self._je - self._js)
        elif self.direction[0] == 'z':
            self.n_fields = self._ke - self._ks
            self.r_magnitude_factor = (1 + self._ie - self._is) * (1 + self._je - self._js) / (self._ke - self._ks)

        if self.direction[1] == 'n':
            self.v_magnitude_factor = -1 * self.magnitude / self.n_fields
        elif self.direction[1] == 'p':
            self.v_magnitude_factor = 1 * self.magnitude / self.n_fields

        self.resistance_per_component = self.r_magnitude_factor * self.resistance
        self.voltage_per_e_field = waveform*self.v_magnitude_factor
        self.waveform = self.v_magnitude_factor * self.n_fields * waveform



class CurrentSource:
    def __init__(self, min_x, min_y, min_z, max_x, max_y, max_z, direction, resistance=50, magnitude=1,
                 waveform_type='gaussian'):
        self.min_x = min_x
        self.min_y = min_y
        self.min_z = min_z
        self.max_x = max_x
        self.max_y = max_y
        self.max_z = max_z
        self.direction = direction
        self.resistance = resistance
        self.magnitude = magnitude
        self.waveform_type = waveform_type

    def resistance(self, resistance=50):
        self.resistance = resistance

    def magnitude(self, magnitude=1):
        self.magnitude = magnitude

    def waveform(self, waveform_type='gaussian'):
        self.waveform_type = waveform_type

    def get_is(self):
        return self._is

    def get_js(self):
        return self._js

    def get_ks(self):
        return self._ks

    def get_ie(self):
        return self._ie

    def get_je(self):
        return self._je

    def get_ke(self):
        return self._ke

    def initWaveform(self,fdtd_domain,yee_cell,waveform):
        self._is = round((self.min_x - fdtd_domain.min_x) / yee_cell['dx'])
        self._js = round((self.min_y - fdtd_domain.min_y) / yee_cell['dy'])
        self._ks = round((self.min_z - fdtd_domain.min_z) / yee_cell['dz'])
        self._ie = round((self.max_x - fdtd_domain.min_x) / yee_cell['dx'])
        self._je = round((self.max_y - fdtd_domain.min_y) / yee_cell['dy'])
        self._ke = round((self.max_z - fdtd_domain.min_z) / yee_cell['dz'])

        if self.direction[0] == 'x':
            self.n_fields = (1 + self._je - self._js) * (1 + self._ke - self._ks)
            self.r_magnitude_factor = (1 + self._je - self._js) * (1 + self._ke - self._ks) / (self._ie - self._is)
        elif self.direction[0] == 'y':
            self.n_fields = (1 + self._ie - self._is) * (1 + self._ke - self._ks)
            self.r_magnitude_factor = (1 + self._ie - self._is) * (1 + self._ke - self._ks) / (self._je - self._js)
        elif self.direction[0] == 'z':
            self.n_fields = (1 + self._ie - self._is) * (1 + self._je - self._js)
            self.r_magnitude_factor = (1 + self._ie - self._is) * (1 + self._je - self._js) / (self._ke - self._ks)

        if self.direction[1] == 'n':
            self.i_magnitude_factor = -1 * self.magnitude / self.n_fields
        elif self.direction[1] == 'p':
            self.i_magnitude_factor = 1 * self.magnitude / self.n_fields

        self.resistance_per_component = self.r_magnitude_factor * self.resistance
        self.current_per_e_field = waveform * self.i_magnitude_factor
        self.waveform = self.i_magnitude_factor * self.n_fields * waveform

class Resistor:
    def __init__(self, min_x, min_y, min_z, max_x, max_y, max_z, direction, resistance):
        self.min_x = min_x
        self.min_y = min_y
        self.min_z = min_z
        self.max_x = max_x
        self.max_y = max_y
        self.max_z = max_z
        self.direction = direction
        self.resistance = resistance

    def get_is(self):
        return self._is

    def get_js(self):
        return self._js

    def get_ks(self):
        return self._ks

    def get_ie(self):
        return self._ie

    def get_je(self):
        return self._je

    def get_ke(self):
        return self._ke

    def initResistor(self,fdtd_domain,yee_cell):
        self._is = round((self.min_x - fdtd_domain.min_x) / yee_cell['dx'])
        self._js = round((self.min_y - fdtd_domain.min_y) / yee_cell['dy'])
        self._ks = round((self.min_z - fdtd_domain.min_z) / yee_cell['dz'])
        self._ie = round((self.max_x - fdtd_domain.min_x) / yee_cell['dx'])
        self._je = round((self.max_y - fdtd_domain.min_y) / yee_cell['dy'])
        self._ke = round((self.max_z - fdtd_domain.min_z) / yee_cell['dz'])

        if self.direction[0] == 'x':
            self.r_magnitude_factor = (1 + self._je - self._js) * (1 + self._ke - self._ks) / (self._ie - self._is)
        elif self.direction[0] == 'y':
            self.r_magnitude_factor = (1 + self._ie - self._is) * (1 + self._ke - self._ks) / (self._je - self._js)
        elif self.direction[0] == 'z':
            self.r_magnitude_factor = (1 + self._ie - self._is) * (1 + self._je - self._js) / (self._ke - self._ks)

        self.resistance_per_component = self.r_magnitude_factor * self.resistance

class Inductor:
    def __init__(self, min_x, min_y, min_z, max_x, max_y, max_z, direction, inductance):
        self.min_x = min_x
        self.min_y = min_y
        self.min_z = min_z
        self.max_x = max_x
        self.max_y = max_y
        self.max_z = max_z
        self.direction = direction
        self.inductance = inductance

    def get_is(self):
        return self._is

    def get_js(self):
        return self._js

    def get_ks(self):
        return self._ks

    def get_ie(self):
        return self._ie

    def get_je(self):
        return self._je

    def get_ke(self):
        return self._ke

    def initInductor(self,fdtd_domain,yee_cell):
        self._is = round((self.min_x - fdtd_domain.min_x) / yee_cell['dx'])
        self._js = round((self.min_y - fdtd_domain.min_y) / yee_cell['dy'])
        self._ks = round((self.min_z - fdtd_domain.min_z) / yee_cell['dz'])
        self._ie = round((self.max_x - fdtd_domain.min_x) / yee_cell['dx'])
        self._je = round((self.max_y - fdtd_domain.min_y) / yee_cell['dy'])
        self._ke = round((self.max_z - fdtd_domain.min_z) / yee_cell['dz'])

        if self.direction[0] == 'x':
            self.l_magnitude_factor = (1 + self._je - self._js) * (1 + self._ke - self._ks) / (self._ie - self._is)
        elif self.direction[0] == 'y':
            self.l_magnitude_factor = (1 + self._ie - self._is) * (1 + self._ke - self._ks) / (self._je - self._js)
        elif self.direction[0] == 'z':
            self.l_magnitude_factor = (1 + self._ie - self._is) * (1 + self._je - self._js) / (self._ke - self._ks)

        self.inductance_per_component = self.l_magnitude_factor * self.inductance


class Capacitor:
    def __init__(self, min_x, min_y, min_z, max_x, max_y, max_z, direction, capacitance):
        self.min_x = min_x
        self.min_y = min_y
        self.min_z = min_z
        self.max_x = max_x
        self.max_y = max_y
        self.max_z = max_z
        self.direction = direction
        self.capacitance = capacitance

    def get_is(self):
        return self._is

    def get_js(self):
        return self._js

    def get_ks(self):
        return self._ks

    def get_ie(self):
        return self._ie

    def get_je(self):
        return self._je

    def get_ke(self):
        return self._ke

    def initCapacitor(self,fdtd_domain,yee_cell):
        self._is = round((self.min_x - fdtd_domain.min_x) / yee_cell['dx'])
        self._js = round((self.min_y - fdtd_domain.min_y) / yee_cell['dy'])
        self._ks = round((self.min_z - fdtd_domain.min_z) / yee_cell['dz'])
        self._ie = round((self.max_x - fdtd_domain.min_x) / yee_cell['dx'])
        self._je = round((self.max_y - fdtd_domain.min_y) / yee_cell['dy'])
        self._ke = round((self.max_z - fdtd_domain.min_z) / yee_cell['dz'])

        if self.direction[0] == 'x':
            self.c_magnitude_factor = (self._ie - self._is) / ((1 + self._je - self._js) * (1 + self._ke - self._ks))
        elif self.direction[0] == 'y':
            self.c_magnitude_factor = (self._je - self._js) / ((1 + self._ie - self._is) * (1 + self._ke - self._ks))
        elif self.direction[0] == 'z':
            self.c_magnitude_factor = (self._ke - self._ks) / ((1 + self._ie - self._is) * (1 + self._je - self._js))

        self.capacitance_per_component = self.l_magnitude_factor * self.capacitance

class Impedance:
    def __init__(self, min_x, min_y, min_z, max_x, max_y, max_z, direction):
        self.min_x = min_x
        self.min_y = min_y
        self.min_z = min_z
        self.max_x = max_x
        self.max_y = max_y
        self.max_z = max_z
        self.direction = direction


class Diode:
    def __init__(self, min_x, min_y, min_z, max_x, max_y, max_z, direction):
        self.min_x = min_x
        self.min_y = min_y
        self.min_z = min_z
        self.max_x = max_x
        self.max_y = max_y
        self.max_z = max_z
        self.direction = direction


class BJT:
    def __init__(self, min_x, min_y, min_z, max_x, max_y, max_z, direction):
        self.min_x = min_x
        self.min_y = min_y
        self.min_z = min_z
        self.max_x = max_x
        self.max_y = max_y
        self.max_z = max_z
        self.direction = direction


class FET:
    def __init__(self, min_x, min_y, min_z, max_x, max_y, max_z, direction):
        self.min_x = min_x
        self.min_y = min_y
        self.min_z = min_z
        self.max_x = max_x
        self.max_y = max_y
        self.max_z = max_z
        self.direction = direction


class Waveform:
    def __init__(self, time):
        self.time = time
        self.waveform = np.zeros(len(self.time))


class SinusoidalWaveform:
    def __init__(self, time, frequency, t_0):
        self.frequency = frequency
        self.time = time
        self.t_0 = t_0
        self.waveform = math.sin(2 * math.pi() * frequency * (time - t_0))


class UnitStepWaveform:
    def __init__(self, start_time_step, number_of_time_steps):
        self.start_time_step = start_time_step
        self.number_of_time_steps = number_of_time_steps
        self.waveform[0:self.number_of_time_steps] = 1
        self.waveform[0:self.start_index - 1] = 0


class GaussianWaveform:
    def __init__(self, time, problem_definition):
        self.time = time
        self.ncpw = problem_definition.number_of_cells_per_wavelength
        self.c = problem_definition.c
        self.yee_cell = problem_definition.yee_cell
        self.max_yee_cell = max([self.yee_cell['dx'], self.yee_cell['dy'], self.yee_cell['dz']])
        self.max_frequency = self.c / (self.ncpw * self.max_yee_cell)
        self.tau = self.ncpw * self.max_yee_cell / (2 * self.c)
        self.t_0 = 4.5 * self.tau
        self.waveform = math.exp(-pow((self.time - self.t_0) / self.tau, 2))


class DerivativeGaussianWaveform:
    def __init__(self, time, problem_definition):
        self.time = time
        self.ncpw = problem_definition.number_of_cells_per_wavelength
        self.c = problem_definition.c
        self.yee_cell = problem_definition.yee_cell
        self.max_yee_cell = max([self.yee_cell['dx'], self.yee_cell['dy'], self.yee_cell['dz']])
        self.max_frequency = self.c / (self.ncpw * self.max_yee_cell)
        self.tau = self.ncpw * self.max_yee_cell / (2 * self.c)
        self.t_0 = 4.5 * self.tau
        self.waveform = -(math.sqrt(2 * math.exp(1)) / self.tau) * (self.time - self.t_0) * math.exp(
            -pow(((self.time - self.t_0) / self.tau), 2))


class CosineModulatedGaussianWaveform:
    def __init__(self, time, center_frequency, bandwidth):
        self.time = time
        self.frequency = center_frequency
        self.bandwidth = bandwidth
        self.tau = 0.966 / bandwidth
        self.t_0 = 4.5 * self.tau
        self.waveform = math.cos(2 * math.pi * self.frequency * (self.time - self.t_0)) * math.exp(
            -pow(((self.time - self.t_0) / self.tau), 2))


class VoltageProbe:
    def __init__(self, min_x, min_y, min_z, max_x, max_y, max_z, direction):
        self.min_x = min_x
        self.min_y = min_y
        self.min_z = min_z
        self.max_x = max_x
        self.max_y = max_y
        self.max_z = max_z
        self.direction = direction


class CurrentProbe:
    def __init__(self, min_x, min_y, min_z, max_x, max_y, max_z, direction):
        self.min_x = min_x
        self.min_y = min_y
        self.min_z = min_z
        self.max_x = max_x
        self.max_y = max_y
        self.max_z = max_z
        self.direction = direction


def make_frequency_list(start, end=0, step_size=0, number_of_points=0):
    if step_size != 0 and end == 0 and number_of_points != 0:
        return [start + step_size * i for i in range(number_of_points)]
    elif step_size == 0 and end != 0 and number_of_points != 0:
        return [start + (end - start) * i / number_of_points for i in range(number_of_points)]
    elif step_size != 0 and end != 0 and number_of_points == 0:
        return [start + step_size * i for i in range(int((end - start) / step_size))]
    else:
        return []


class Port:
    def __init__(self, voltage_probe_index, current_probe_index, impedance, source_port_enable):
        self.voltageProbeIndex = voltage_probe_index
        self.currentProbeIndex = current_probe_index
        self.impedance = impedance
        self.sourcePortEnable = source_port_enable


class Definition:
    def __init__(self):
        self.number_of_time_steps = 10000
        self.courant_factor = 0.9
        self.number_of_cells_per_wavelength = 20
        self.yee_cell = {'dx': 0.5e-3,
                         'dy': 0.5e-3,
                         'dz': 0.5e-3}
        self.boundary = {'cpml_xn': True,
                         'air_buffer_number_of_cells_xn': 10,
                         'cpml_number_of_cells_xn': 8,
                         'pbc_xn': True,
                         'cpml_xp': True,
                         'air_buffer_number_of_cells_xp': 10,
                         'cpml_number_of_cells_xp': 8,
                         'pbc_xp': True,
                         'cpml_yn': True,
                         'air_buffer_number_of_cells_yn': 10,
                         'cpml_number_of_cells_yn': 8,
                         'pbc_yn': True,
                         'cpml_yp': True,
                         'air_buffer_number_of_cells_yp': 10,
                         'cpml_number_of_cells_yp': 8,
                         'pbc_yp': True,
                         'cpml_zn': True,
                         'air_buffer_number_of_cells_zn': 10,
                         'cpml_number_of_cells_zn': 8,
                         'pbc_zn': True,
                         'cpml_zp': True,
                         'air_buffer_number_of_cells_zp': 10,
                         'cpml_number_of_cells_zp': 8,
                         'pbc_zp': True,
                         'cpml_order': 3,
                         'cpml_sigma_factor': 1.3,
                         'cpml_kappa_max': 7,
                         'cpml_alpha_min': 0,
                         'cpml_alpha_max': 0.05}
        self.materials = {'air': {'index': 1, 'property': MaterialType(1, 1, 0, 0)},
                          'pec': {'index': 2, 'property': MaterialType(1, 1, 1e10, 0)},
                          'pmc': {'index': 3, 'property': MaterialType(1, 1, 0, 1e10)},
                          'RO5380': {'index': 4, 'property': MaterialType(2.2, 1, 0, 0)}}
        self.bricks = {'sim_box': {'air': Brick(-2e-3, -2e-3, -2e-3, 2e-3, 2e-3, 2e-3)},
                       'substrate': {'RO5380': Brick(-2e-3, -2e-3, -2e-3, 2e-3, 2e-3, 0)},
                       'ground': {'pec': Brick(-2e-3, -2e-3, -2e-3, 2e-3, 2e-3, -2e-3)},
                       'patch': {'pec': Brick(-2e-3, -2e-3, 0, 2e-3, 2e-3, 0)}}
        self.excitations = {'H': VoltageSource(-2e-3, -2e-3, 0, 2e-3, 2e-3, 0, 'zp'),
                            'V': VoltageSource(-2e-3, -2e-3, 0, 2e-3, 2e-3, 0, 'zp')}
        self.frequencies = {'port parameter': make_frequency_list(start=2.5e9, end=3.1e9, step_size=0.1e9),
                            'radiation': make_frequency_list(start=2.5e9, end=3.1e9, step_size=0.1e9)}
        self.voltageProbe = {'H': VoltageProbe(0, 5.25e-3, -1.575e-3, 0, 5.25e-3, 0, 'zp'),
                             'V': VoltageProbe(5.25e-3, 0, -1.575e-3, 5.25e-3, 0, 0, 'zp')}
        self.currentProbe = {'H': CurrentProbe(0, 5.25e-3, -1e-3, 0, 5.25e-3, -1e-3, 'zp'),
                             'V': CurrentProbe(5.25e-3, 0, -1e-3, 5.25e-3, 0, -1e-3, 'zp')}
        self.ports = {'H': Port(1, 1, 50, True),
                      'V': Port(1, 1, 50, True)}

    def set_number_of_time_steps(self, number_of_time_steps):
        self.number_of_time_steps = number_of_time_steps

    def set_courant_factor(self, courant_factor):
        self.courant_factor = courant_factor

    def set_number_of_cells_per_wavelength(self, number_of_cells_per_wavelength):
        self.number_of_cells_per_wavelength = number_of_cells_per_wavelength

    ### Yee Cell ###
    def set_yee_cell_x(self, x):
        self.yee_cell['dx'] = x

    def set_yee_cell_y(self, y):
        self.yee_cell['dy'] = y

    def set_yee_cell_z(self, z):
        self.yee_cell['dz'] = z

    ### Boundary Xn ###
    def set_boundary_cpml_xn(self, state):
        self.boundary['cpml_xn'] = state

    def set_air_buffer_number_of_cells_xn(self, number_of_cells):
        self.boundary['air_buffer_number_of_cells_xn'] = number_of_cells

    def set_cpml_number_of_cells_xn(self, number_of_cells):
        self.boundary['cpml_number_of_cells_xn'] = number_of_cells

    def set_pbc_xn(self, state):
        self.boundary['pbc_xn'] = state

    ### Boundary Xp ###
    def set_boundary_cpml_xp(self, state):
        self.boundary['cpml_xp'] = state

    def set_air_buffer_number_of_cells_xp(self, number_of_cells):
        self.boundary['air_buffer_number_of_cells_xp'] = number_of_cells

    def set_cpml_number_of_cells_xp(self, number_of_cells):
        self.boundary['cpml_number_of_cells_xp'] = number_of_cells

    def set_pbc_xp(self, state):
        self.boundary['pbc_xp'] = state

    ### Boundary Yn ###
    def set_boundary_cpml_yn(self, state):
        self.boundary['cpml_yn'] = state

    def set_air_buffer_number_of_cells_yn(self, number_of_cells):
        self.boundary['air_buffer_number_of_cells_yn'] = number_of_cells

    def set_cpml_number_of_cells_yn(self, number_of_cells):
        self.boundary['cpml_number_of_cells_yn'] = number_of_cells

    def set_pbc_yn(self, state):
        self.boundary['pbc_yn'] = state

    ### Boundary Yp ###
    def set_boundary_cpml_yp(self, state):
        self.boundary['cpml_yp'] = state

    def set_air_buffer_number_of_cells_yp(self, number_of_cells):
        self.boundary['air_buffer_number_of_cells_yp'] = number_of_cells

    def set_cpml_number_of_cells_yp(self, number_of_cells):
        self.boundary['cpml_number_of_cells_yp'] = number_of_cells

    def set_pbc_yp(self, state):
        self.boundary['pbc_yp'] = state

    ### Boundary Zn ###
    def set_boundary_cpml_zn(self, state):
        self.boundary['cpml_zn'] = state

    def set_air_buffer_number_of_cells_zn(self, number_of_cells):
        self.boundary['air_buffer_number_of_cells_zn'] = number_of_cells

    def set_cpml_number_of_cells_zn(self, number_of_cells):
        self.boundary['cpml_number_of_cells_zn'] = number_of_cells

    def set_pbc_zn(self, state):
        self.boundary['pbc_zn'] = state

    ### Boundary Zp ###
    def set_boundary_cpml_zp(self, state):
        self.boundary['cpml_zp'] = state

    def set_air_buffer_number_of_cells_zp(self, number_of_cells):
        self.boundary['air_buffer_number_of_cells_zp'] = number_of_cells

    def set_cpml_number_of_cells_zp(self, number_of_cells):
        self.boundary['cpml_number_of_cells_zp'] = number_of_cells

    def set_pbc_zp(self, state):
        self.boundary['pbc_zp'] = state

    ### CPML Parameter ###
    def set_cpml_order(self, order):
        self.boundary['cpml_order'] = order

    def set_cpml_sigma_factor(self, factor):
        self.boundary['cpml_sigma_factor'] = factor

    def set_cpml_kappa_max(self, kappa_max):
        self.boundary['cpml_kappa_max'] = kappa_max

    def set_cpml_alpha_min(self, alpha_min):
        self.boundary['cpml_alpha_min'] = alpha_min

    def set_cpml_alpha_max(self, alpha_max):
        self.boundary['cpml_alpha_max'] = alpha_max

    ### initializing FDTD parameters ###
    def init_fdtd_param(self):
        self.eps_0 = 8.854187817e-12
        self.mu_0 = 4 * math.pi * 1e-7
        self.c = 1 / math.sqrt(self.mu_0 * self.eps_0)
        self.dt = self.courant_factor / (self.c * math.sqrt((1 / self.yee_cell.dx ^ 2) +
                                                            (1 / self.yee_cell.dy ^ 2) +
                                                            (1 / self.yee_cell.dz ^ 2)))
        self.time = np.array([(i + 0.5) * self.dt for i in range(self.number_of_time_steps)])
