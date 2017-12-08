import numpy as np


class DomainFDTD:
    def __init__(self, problem_definition):
        self.size_x = 0.0
        self.size_y = 0.0
        self.size_z = 0.0
        self.max_x = 0.0
        self.max_y = 0.0
        self.max_z = 0.0
        self.min_x = 0.0
        self.min_y = 0.0
        self.min_z = 0.0
        self.nx = 0
        self.ny = 0
        self.nz = 0
        self.cell_center_coordinates_x = []
        self.cell_center_coordinates_y = []
        self.cell_center_coordinates_z = []
        self.bricks = problem_definition.bricks
        self.boundary = problem_definition.boundary
        self.yee_cell = problem_definition.yee_cell
        self.materials = problem_definition.materials

    def check_sim_box(self, media='air'):
        for zoneKey in self.bricks.keys():
            for mediaKey in self.bricks[zoneKey].keys():
                check = self.bricks[zoneKey][mediaKey].max_x >= self.bricks['sim_box'][media].min_x
                check = check and self.bricks[zoneKey][mediaKey].max_y >= self.bricks['sim_box'][media].min_y
                check = check and self.bricks[zoneKey][mediaKey].max_z >= self.bricks['sim_box'][media].min_z
                check = check and self.bricks[zoneKey][mediaKey].min_x <= self.bricks['sim_box'][media].max_x
                check = check and self.bricks[zoneKey][mediaKey].min_y <= self.bricks['sim_box'][media].max_y
                check = check and self.bricks[zoneKey][mediaKey].min_z <= self.bricks['sim_box'][media].max_z
                if not check:
                    break
        return check

    def calculate_cell_center_coordinates(self):
        self.cell_center_coordinates_x = np.zeros(self.nx, self.ny, self.nz);
        self.cell_center_coordinates_y = np.zeros(self.nx, self.ny, self.nz);
        self.cell_center_coordinates_z = np.zeros(self.nx, self.ny, self.nz);
        for i in range(self.nx):
            self.cell_center_coordinates_x[i][:][:] = (i - 0.5) * self.dx + self.min_x
        for i in range(self.ny):
            self.cell_center_coordinates_y[:][i][:] = (i - 0.5) * self.dy + self.min_y
        for i in range(self.nx):
            self.cell_center_coordinates_z[:][:][i] = (i - 0.5) * self.dz + self.min_z

    def air_buffer_size(self,param):
        return self.yee_cell['d' + param[0]] * self.boundary['air_buffer_number_of_cells_' + param]

    def calculate_domain_limits(self,media='air'):
        if self.check_sim_box(media):
            self.min_x = self.bricks['sim_box'][media].min_x - self.air_buffer_size('xn')
            self.min_y = self.bricks['sim_box'][media].min_y - self.air_buffer_size('yn')
            self.min_z = self.bricks['sim_box'][media].min_z - self.air_buffer_size('zn')
            self.max_x = self.bricks['sim_box'][media].max_x + self.air_buffer_size('xp')
            self.max_y = self.bricks['sim_box'][media].max_y + self.air_buffer_size('yp')
            self.max_z = self.bricks['sim_box'][media].max_z + self.air_buffer_size('zp')
            self.size_x = self.min_x - self.max_x
            self.size_y = self.min_y - self.max_y
            self.size_z = self.min_z - self.max_z
            self.nx = round(self.size_x/self.yee_cell['dx'])
            self.ny = round(self.size_y/self.yee_cell['dy'])
            self.nz = round(self.size_z/self.yee_cell['dz'])
            self.size_x = self.nx*self.yee_cell['dx']
            self.size_y = self.ny*self.yee_cell['dy']
            self.size_z = self.nz*self.yee_cell['dz']
            self.max_x = self.size_x + self.min_x
            self.max_y = self.size_y + self.min_y
            self.max_z = self.size_z + self.min_z
            return self.nx, self.ny, self.nz
        else:
            print('The simulation box should enclose the structure under simulation')
            exit()

    def create_material_space(self):
        material_3d_space = np.ones([self.nx, self.ny, self.nz], dtype=float);
        for key_box in self.bricks.keys():
            for key_mat in self.bricks[key_box].keys():
                blx = round((self.bricks[key_box][key_mat].min_x - self.min_x) / self.yee_cell['dx']) + 1
                bly = round((self.bricks[key_box][key_mat].min_y - self.min_y) / self.yee_cell['dy']) + 1
                blz = round((self.bricks[key_box][key_mat].min_z - self.min_z) / self.yee_cell['dz']) + 1
                bux = round((self.bricks[key_box][key_mat].max_x - self.min_x) / self.yee_cell['dx']) + 1
                buy = round((self.bricks[key_box][key_mat].max_y - self.min_y) / self.yee_cell['dy']) + 1
                buz = round((self.bricks[key_box][key_mat].max_z - self.min_z) / self.yee_cell['dz']) + 1
                material_3d_space[blx:bux-1,bly:buy-1,blz:buz-1] = self.materials[self.bricks[key_box][key_mat]]['index']
        return material_3d_space


class InitMaterialParamSet:
    def __init__(self,problem_definition):
        print('initializing the problem')
        self.domain = DomainFDTD(problem_definition)
        self.nx, self.ny, self.nz = self.domain.calculate_domain_limits()
        self.material_space = self.domain.create_material_space();
        self.property_vectors = {'eps_r' : np.array([]),
                                 'mu_r' : np.array([]),
                                 'sigma_e' : np.array([]),
                                 'sigma_m' : np.array([])}
        self.eps_r_x = np.ones([self.nx, self.ny + 1, self.nz + 1], dtype=float)
        self.eps_r_y = np.ones([self.nx + 1, self.ny, self.nz + 1], dtype=float)
        self.eps_r_z = np.ones([self.nx + 1, self.ny + 1, self.nz], dtype=float)
        self.mu_r_x = np.ones([self.nx + 1, self.ny, self.nz], dtype=float)
        self.mu_r_y = np.ones([self.nx, self.ny + 1, self.nz], dtype=float)
        self.mu_r_z = np.ones([self.nx, self.ny, self.nz + 1], dtype=float)
        self.sigma_e_x = np.zeros([self.nx, self.ny + 1, self.nz + 1], dtype=float)
        self.sigma_e_y = np.zeros([self.nx + 1, self.ny, self.nz + 1], dtype=float)
        self.sigma_e_z = np.zeros([self.nx + 1, self.ny + 1, self.nz], dtype=float)
        self.sigma_m_x = np.zeros([self.nx + 1, self.ny, self.nz], dtype=float)
        self.sigma_m_y = np.zeros([self.nx, self.ny + 1, self.nz], dtype=float)
        self.sigma_m_z = np.zeros([self.nx, self.ny, self.nz + 1], dtype=float)

    def init_property_vectors(self):
        for key in self.domain.materials.keys():
            if not (self.domain.materials[key]['property'].eps_r == 0):
                self.property_vectors['eps_r'][key] = self.domain.materials[key]['property'].eps_r;
            else:
                self.property_vectors['eps_r'][key] = 1
            if not (self.domain.materials[key]['property'].mu_r == 0):
                self.property_vectors['mu_r'][key] = self.domain.materials[key]['property'].mu_r;
            else:
                self.property_vectors['mu_r'][key] = 1e-20
            if not (self.domain.materials[key]['property'].sigma_e == 0):
                self.property_vectors['sigma_e'][key] = self.domain.materials[key]['property'].sigma_e;
            else:
                self.property_vectors['sigma_e'][key] = 1e-20
            if not (self.domain.materials[key]['property'].sigma_m == 0):
                self.property_vectors['sigma_m'][key] = self.domain.materials[key]['property'].sigma_m;
            else:
                self.property_vectors['sigma_m'][key] = 1e-20

    def calculate_effective_electric_property_x(self, property_key):
        ep1 = np.zeros(self.nx, self.ny, self.nz)
        ep2 = np.zeros(self.nx, self.ny, self.nz)
        ep3 = np.zeros(self.nx, self.ny, self.nz)
        ep4 = np.zeros(self.nx, self.ny, self.nz)
        ep1[0: self.nx, 1: self.ny, 1: self.nz] = self.property_vector[property_key][self.material_space[0: self.nx, 1: self.ny, 1: self.nz]]
        ep2[0: self.nx, 0: self.ny - 1, 1: self.nz] = self.property_vector[property_key][self.material_space[0: self.nx, 0: self.ny - 1, 1: self.nz]]
        ep3[0: self.nx, 1: self.ny, 0: self.nz - 1] = self.property_vector[property_key][self.material_space[0: self.nx, 1: self.ny, 0: self.nz - 1]]
        ep4[0: self.nx, 0: self.ny - 1, 0: self.nz - 1] = self.property_vector[property_key][self.material_space[0: self.nx, 0: self.ny - 1, 0: self.nz - 1]]
        if property_key == 'eps_r':
            self.eps_r_x [0:self.nx,1:self.ny,1:self.nz] = 0.25*(ep1 + ep2 + ep3 + ep4)
        elif property_key == 'sigma_e':
            self.sigma_e_x [0:self.nx,1:self.ny,1:self.nz]= 0.25 * (ep1 + ep2 + ep3 + ep4)

    def calculate_effective_electric_property_y(self, property_key):
        ep1 = np.zeros(self.nx, self.ny, self.nz)
        ep2 = np.zeros(self.nx, self.ny, self.nz)
        ep3 = np.zeros(self.nx, self.ny, self.nz)
        ep4 = np.zeros(self.nx, self.ny, self.nz)
        ep1[1: self.nx, 0: self.ny, 1: self.nz] = self.property_vector[property_key][self.material_space[1: self.nx, 0: self.ny, 1: self.nz]]
        ep2[0: self.nx-1, 0: self.ny, 1: self.nz] = self.property_vector[property_key][self.material_space[0: self.nx-1, 0: self.ny, 1: self.nz]]
        ep3[1: self.nx, 0: self.ny, 0: self.nz - 1] = self.property_vector[property_key][self.material_space[1: self.nx, 0: self.ny, 0: self.nz - 1]]
        ep4[0: self.nx-1, 0: self.ny, 0: self.nz - 1] = self.property_vector[property_key][self.material_space[0: self.nx-1, 0: self.ny, 0: self.nz - 1]]
        if property_key == 'eps_r':
            self.eps_r_y [1:self.nx,0:self.ny,1:self.nz] = 0.25*(ep1 + ep2 + ep3 + ep4)
        elif property_key == 'sigma_e':
            self.sigma_e_y [1:self.nx,0:self.ny,1:self.nz]= 0.25 * (ep1 + ep2 + ep3 + ep4)

    def calculate_effective_electric_property_z(self, property_key):
        ep1 = np.zeros(self.nx, self.ny, self.nz)
        ep2 = np.zeros(self.nx, self.ny, self.nz)
        ep3 = np.zeros(self.nx, self.ny, self.nz)
        ep4 = np.zeros(self.nx, self.ny, self.nz)
        ep1[0: self.nx, 1: self.ny, 0: self.nz] = self.property_vector[property_key][self.material_space[1: self.nx, 1: self.ny, 0: self.nz]]
        ep2[0: self.nx-1, 1: self.ny, 0: self.nz] = self.property_vector[property_key][self.material_space[0: self.nx-1, 1: self.ny, 0: self.nz]]
        ep3[1: self.nx, 0: self.ny-1, 0: self.nz] = self.property_vector[property_key][self.material_space[1: self.nx, 0: self.ny-1, 0: self.nz]]
        ep4[0: self.nx-1, 0: self.ny - 1, 0: self.nz] = self.property_vector[property_key][self.material_space[0: self.nx-1, 0: self.ny - 1, 0: self.nz]]
        if property_key == 'eps_r':
            self.eps_r_z [1:self.nx,1:self.ny,0:self.nz] = 0.25*(ep1 + ep2 + ep3 + ep4)
        elif property_key == 'sigma_e':
            self.sigma_e_z [1:self.nx,1:self.ny,0:self.nz]= 0.25 * (ep1 + ep2 + ep3 + ep4)

    def calculate_effective_magnatic_property_x(self,property_key):
        mp1 = np.zeros(self.nx, self.ny, self.nz)
        mp2 = np.zeros(self.nx, self.ny, self.nz)
        mp3 = np.zeros(self.nx, self.ny, self.nz)
        mp4 = np.zeros(self.nx, self.ny, self.nz)
        mp1[1:self.nx,0:self.ny,0:self.nz] = self.property_vector[property_key][self.material_space[1:self.nx,0:self.ny,0:self.nz]]
        mp2[0:self.nx-1,0:self.ny,0:self.nz] = self.property_vector[property_key][self.material_space[0:self.nx-1,0:self.ny,0:self.nz]]
        mp3[1:self.nx,0:self.ny,0:self.nz] = self.property_vector[property_key][self.material_space[1:self.nx,0:self.ny,0:self.nz]]
        mp4[0:self.nx-1,0:self.ny,0:self.nz] = self.property_vector[property_key][self.material_space[0:self.nx-1,0:self.ny,0:self.nz]]
        if property_key == 'mu_r':
            self.mu_r_x[1:self.nx, 0:self.ny, 0:self.nz] = 2*(mp1* mp2)/(mp3+ mp4)
        elif property_key == 'sigma_m':
            self.sigma_m_x[1:self.nx, 0:self.ny, 0:self.nz] = 2*(mp1* mp2)/(mp3+ mp4)

    def calculate_effective_magnatic_property_y(self,property_key):
        mp1 = np.zeros(self.nx, self.ny, self.nz)
        mp2 = np.zeros(self.nx, self.ny, self.nz)
        mp3 = np.zeros(self.nx, self.ny, self.nz)
        mp4 = np.zeros(self.nx, self.ny, self.nz)
        mp1[0:self.nx,1:self.ny,0:self.nz] = self.property_vector[property_key][self.material_space[0:self.nx,1:self.ny,0:self.nz]]
        mp2[0:self.nx,0:self.ny-1,0:self.nz] = self.property_vector[property_key][self.material_space[0:self.nx,0:self.ny-1,0:self.nz]]
        mp3[0:self.nx,1:self.ny,0:self.nz] = self.property_vector[property_key][self.material_space[0:self.nx,1:self.ny,0:self.nz]]
        mp4[0:self.nx,0:self.ny-1,0:self.nz] = self.property_vector[property_key][self.material_space[0:self.nx,0:self.ny-1,0:self.nz]]
        if property_key == 'mu_r':
            self.mu_r_y[0:self.nx, 1:self.ny, 0:self.nz] = 2*(mp1* mp2)/(mp3+ mp4)
        elif property_key == 'sigma_m':
            self.sigma_m_y[0:self.nx, 1:self.ny, 0:self.nz] = 2*(mp1* mp2)/(mp3+ mp4)

    def calculate_effective_magnatic_property_z(self,property_key):
        mp1 = np.zeros(self.nx, self.ny, self.nz)
        mp2 = np.zeros(self.nx, self.ny, self.nz)
        mp3 = np.zeros(self.nx, self.ny, self.nz)
        mp4 = np.zeros(self.nx, self.ny, self.nz)
        mp1[0:self.nx,0:self.ny,1:self.nz] = self.property_vector[property_key][self.material_space[0:self.nx,0:self.ny,1:self.nz]]
        mp2[0:self.nx,0:self.ny,0:self.nz-1] = self.property_vector[property_key][self.material_space[0:self.nx,0:self.ny,0:self.nz-1]]
        mp3[0:self.nx,0:self.ny,1:self.nz] = self.property_vector[property_key][self.material_space[0:self.nx,0:self.ny,1:self.nz]]
        mp4[0:self.nx,0:self.ny,0:self.nz-1] = self.property_vector[property_key][self.material_space[0:self.nx,0:self.ny,0:self.nz-1]]
        if property_key == 'mu_r':
            self.mu_r_z[0:self.nx, 0:self.ny, 1:self.nz] = 2*(mp1* mp2)/(mp3+ mp4)
        elif property_key == 'sigma_m':
            self.sigma_m_z[0:self.nx, 0:self.ny, 1:self.nz] = 2*(mp1* mp2)/(mp3+ mp4)

    def create_zero_thickness_pec_plate(self):
        bricks = self.domain.bricks
        materials = self.domain.materials
        yee_cell = self.domain.yee_cell
        for key in bricks.keys():
            sigma = materials[bricks[key].keys()]['property'].sigma_e
            blx = round((bricks[key].min_x - self.domain.min_x) / yee_cell.dx)
            bly = round((bricks[key].min_y - self.domain.min_y) / yee_cell.dy)
            blz = round((bricks[key].min_z - self.domain.min_z) / yee_cell.dz)
            bux = round((bricks[key].max_x - self.domain.min_x) / yee_cell.dx)
            buy = round((bricks[key].max_y - self.domain.min_y) / yee_cell.dy)
            buz = round((bricks[key].max_z - self.domain.min_z) / yee_cell.dz)
            if blx == bux:
                self.sigma_e_y[blx, bly:buy-1,blz:buz] = sigma
                self.sigma_e_z[blx, bly:buy,blz:buz-1] = sigma
            if bly == buy:
                self.sigma_e_z[blx: bux, bly, blz: buz - 1] = sigma
                self.sigma_e_x[blx: bux - 1, bly, blz: buz] = sigma
            if blz == buz:
                self.sigma_e_x[blx: bux - 1, bly: buy, blz] = sigma
                self.sigma_e_y[blx: bux, bly: buy - 1, blz] = sigma

    def get_material_3d_space(self):
        return self.material_space

    def get_fdtd_domain(self):
        return self.domain

    def get_eps_params(self):
        return self.eps_r_x,self.eps_r_y,self.eps_r_z

    def get_mu_params(self):
        return self.mu_r_x,self.mu_r_y,self.mu_r_z

    def get_sigma_e_params(self):
        return self.sigma_e_x,self.sigma_e_y,self.sigma_e_z

    def get_sigma_m_params(self):
        return self.sigma_m_x,self.sigma_m_y,self.sigma_m_z


class ElectromagneticFieldComponents:
    def __init__(self,nx,ny,nz):
        self.Hx = np.zeros(nx+1, ny, nz)
        self.Hy = np.zeros(nx, ny+1, nz)
        self.Hz = np.zeros(nx, ny, nz+1)
        self.Ex = np.zeros(nx, ny+1, nz+1)
        self.Ey = np.zeros(nx+1, ny, nz+1)
        self.Ez = np.zeros(nx+1, ny+1, nz)

class UpdatingCoefficients:
    def __init__(self,mps,eps_0,mu_0,dt,dx,dy,dz):
        self.eps_0 = eps_0
        self.mu_0 = mu_0
        self.dt = dt
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.Cexe = (2 * mps.eps_r_x * eps_0 - dt * mps.sigma_e_x)/(2 * mps.eps_r_x * eps_0 + dt * mps.sigma_e_x)
        self.Cexhz = (2 * dt / dy)/(2 * mps.eps_r_x * eps_0 + dt * mps.sigma_e_x)
        self.Cexhy = -(2 * dt / dz)/(2 * mps.eps_r_x * eps_0 + dt * mps.sigma_e_x)

        self.Ceye = (2 * mps.eps_r_y * eps_0 - dt * mps.sigma_e_y)/(2 * mps.eps_r_y * eps_0 + dt * mps.sigma_e_y)
        self.Ceyhx = (2 * dt / dz)/ (2 * mps.eps_r_y * eps_0 + dt * mps.sigma_e_y)
        self.Ceyhz = -(2 * dt / dx)/ (2 * mps.eps_r_y * eps_0 + dt * mps.sigma_e_y)

        self.Ceze = (2 * mps.eps_r_z * eps_0 - dt * mps.sigma_e_z)/(2 * mps.eps_r_z * eps_0 + dt * mps.sigma_e_z)
        self.Cezhy = (2 * dt / dx)/ (2 * mps.eps_r_z * eps_0 + dt * mps.sigma_e_z)
        self.Cezhx = -(2 * dt / dy)/ (2 * mps.eps_r_z * eps_0 + dt * mps.sigma_e_z)

        self.Chxh = (2 * mps.mu_r_x * mu_0 - dt * mps.sigma_m_x)/(2 * mps.mu_r_x * mu_0 + dt * mps.sigma_m_x)
        self.Chxez = -(2 * dt / dy)/ (2 * mps.mu_r_x * mu_0 + dt * mps.sigma_m_x)
        self.Chxey = (2 * dt / dz)/ (2 * mps.mu_r_x * mu_0 + dt * mps.sigma_m_x)

        self.Chyh = (2 * mps.mu_r_y * mu_0 - dt * mps.sigma_m_y)/(2 * mps.mu_r_y * mu_0 + dt * mps.sigma_m_y)
        self.Chyex = -(2 * dt / dz)/ (2 * mps.mu_r_y * mu_0 + dt * mps.sigma_m_y)
        self.Chyez = (2 * dt / dx)/ (2 * mps.mu_r_y * mu_0 + dt * mps.sigma_m_y)

        self.Chzh = (2 * mps.mu_r_z * mu_0 - dt * mps.sigma_m_z)/(2 * mps.mu_r_z * mu_0 + dt * mps.sigma_m_z)
        self.Chzey = -(2 * dt / dx)/ (2 * mps.mu_r_z * mu_0 + dt * mps.sigma_m_z)
        self.Chzex = (2 * dt / dy)/ (2 * mps.mu_r_z * mu_0 + dt * mps.sigma_m_z)

    def vs_coefficients(self,vs):
        for k in vs.keys():
            _is = vs[k].get_is()
            _js = vs[k].get_js()
            _ks = vs[k].get_ks()
            _ie = vs[k].get_ie()
            _je = vs[k].get_je()
            _ke = vs[k].get_ke()
            _Res = vs[k].resistance_per_component

            if vs[k].direction[0] == 'x':
                self.Cexe[] = (2 * self.eps_0 * self.eps_r_x[] - self.dt * self.sigma_e_x[] - (self.dt*self.dx)/(_Res*self.dy*self.dz))/(2 * self.eps_0 * self.eps_r_x[] + self.dt * self.sigma_e_x[] + (self.dt*self.dx)/(self.R*self.dy*self.dz))
                self.Cexhz[] = (2 * dt / dy)/(2 * eps_r_x(fi) * eps_0 + dt * sigma_e_x(fi) + (dt*dx)/(R*dy*dz))
                self.Cexhy[] = -(2 * dt / dz)/(2 * eps_r_x(fi) * eps_0 + dt * sigma_e_x(fi) + (dt*dx)/(R*dy*dz))
                vs[k].Cexs = -(2 * dt / (R * dy * dz))/(2 * eps_r_x(fi) * eps_0 + dt * sigma_e_x(fi) + (dt*dx)/(R*dy*dz))