from sys import platform as pltm
import os, platform, pkg_resources
from ctypes import CDLL, c_int, c_double, c_char_p, POINTER
import numpy as np
import json
import _ctypes

class ASRE_Timoshenko_model:
    """
    A class to represent an ASRE Timoshenko beam model.
    
    Attributes
    ----------
    nnode : int
        Number of node.
    meshX : np.array(dtype=np.float64)
        The x coordinates of the beam nodes. Unit: m
    meshY : np.array(dtype=np.float64)
        The y coordinates of the beam nodes. Unit: m
    meshZ : np.array(dtype=np.float64)
        The z coordinates of the beam nodes. Unit: m
    dfoot : float
        The depth of the beam. Unit: m
    bfoot : float
        The width of the beam. Unit: m
    Eb : float
        The Young's modulus of the beam. Unit: N/m^2
    EoverG : float
        The ratio between Young's modulus and the shear modulus. Unit: unitless
    ni_foot : float
        The Poisson's ratio of beam. Unit: unitless
    q_foot : float
        The uniform weight applied on the beam. Unit: N/m (Weight along unit length in the longitudinal direction)
    EsNominal : float
        The nominal elastic modulus of soil. Unit: N/m^2
    nis : float
        The Poisson's ratio of soil. Unit: unitless 
    mu_int : float
        The friction coefficient between soil and beam. Unit: unitless
    """

    def __init__(self, nnode, meshX, meshY, meshZ, dfoot, bfoot, d_NA = 0, solver = 'elasto-plastic'):
        """
        Constructs the beam with dimension properties.

        Parameters
        ----------
        nnode : int
            Number of node.
        meshX : np.array(dtype=np.float64)
            The x coordinates of the beam nodes. Unit: m
        meshY : np.array(dtype=np.float64)
            The y coordinates of the beam nodes. Unit: m
        meshZ : np.array(dtype=np.float64)
            The z coordinates of the beam nodes. Unit: m
        dfoot : float
            The depth of the beam. Unit: m
        bfoot : float
            The width of the beam. Unit: m
        d_NA : float
            The distance from beam neutural axis to ground. Unit: m
        solver : str
            The solver used in the model. Default is 'elasto-plastic'.
            Another option is 'elastic'.
        """
        self.nnode = nnode
        self.meshX = meshX
        self.meshY = meshY
        self.meshZ = meshZ
        self.dfoot = dfoot
        self.bfoot = bfoot
        self.asre_dll = self._import_dll()
        self.d_NA = d_NA
        self.solver = solver

    def _import_dll(self):
        """
        import the CDLL.
        
        Returns
        -------
        CDLL
            The CDLL.
        """
        if pltm == "linux" or pltm == "linux2":
            lib_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "ASREcpp", "bin", "macOS_arm", "libASRElibTimoBeam.so")
            if os.path.exists(lib_path):
                c_lib = CDLL(lib_path)
            else:
                c_lib = None
                message = f'ASRE is not precompiled for {pltm}, please compile the ASRE cpp library'
        elif pltm == "darwin":
            if platform.processor().startswith('i'): # Intel processor
                lib_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "ASREcpp", "bin", "macOS", "libASRElibTimoBeam.dylib")
                if os.path.exists(lib_path):
                    c_lib = CDLL(lib_path)
                else:
                    c_lib = None
                    message = f'ASRE is not precompiled for {pltm} {platform.processor()}, please compile the ASRE cpp library'
            else:
                lib_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "ASREcpp", "bin", "macOS_arm", "libASRElibTimoBeam.dylib")
                if os.path.exists(lib_path):
                    c_lib = CDLL(lib_path)
                else:
                    c_lib = None
                    message = f'ASRE is not precompiled for {pltm} {platform.processor()}, please compile the ASRE cpp library'
        elif pltm == "win32":
            lib_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "ASREcpp", "bin", "win32", "ASRElibTimoBeam.dll")
            if os.path.exists(lib_path):
                c_lib = CDLL(lib_path)
            else:
                c_lib = None
                message = f'ASRE is not precompiled for {pltm}, please compile the ASRE cpp library'
            # c_lib = CDLL(pkg_resources.resource_filename('ASREp', 'ASREcpp//bin//win32//ASRElib.dll'))
            # c_lib.printName()
        if c_lib is None:
            raise ImportError(message)
        c_lib.run.argtypes = [c_int, #nnode 
                              np.ctypeslib.ndpointer(dtype=np.float64), #meshX
                              np.ctypeslib.ndpointer(dtype=np.float64), #meshY
                              np.ctypeslib.ndpointer(dtype=np.float64), #meshZ
                              np.ctypeslib.ndpointer(dtype=np.float64), #dispV
                              np.ctypeslib.ndpointer(dtype=np.float64), #dispL
                              np.ctypeslib.ndpointer(dtype=np.float64), #dispT
                              c_double, #Eb
                              c_double, #EoverG
                              c_double, #EsNominal
                              c_double, #nis
                              c_double, #dfoot
                              c_double, #bfoot
                              c_double, #ni_foot
                              c_double, #mu_int
                              c_double, #qz_foot
                              c_double, #d_NA
                              c_char_p, #solver
                              c_char_p, #output
                                  POINTER(c_double), # result
                                  c_int # result_size
                            ]
        c_lib.run.restype = c_int

        c_lib.KBern3D_foot_TIM_dNA_interface.argtypes = [c_double, # E
                                               c_double, # d_foot
                                               c_double, # b_foot
                                               c_double, # dx
                                               c_double, # EGratio
                                               c_double, # ni_str
                                               c_double, # d_NA
                                               POINTER(c_double) # result_array
                                               ]
        c_lib.KBern3D_foot_TIM_dNA_interface.restype = c_int
        return c_lib

    def KBern3D_foot_TIM_dNA(self, elem_i = 0, d_NA = None):
        """
        Call the C++ function KBern3D_foot_TIM_dNA_interface.

        Parameters
        ----------
        elem_i : int, optional, the index of the element to compute the stiffness matrix for.
            If None, compute the stiffness matrix for the first element.

        Returns
        -------
        np.array(dtype=float)
            The element stiffness matrix.
        """
        E = self.Eb
        d_foot = self.dfoot
        b_foot = self.bfoot
        if elem_i >= self.nnode - 1:
            raise IndexError(f'elem_i {elem_i} is out of range for nnode {self.nnode}')
        dx = self.meshX[elem_i+1] - self.meshX[elem_i]
        EGratio = self.EoverG
        ni_str = self.ni_foot
        if d_NA is None:
            d_NA = self.d_NA

        result_array = (c_double * 144)(*([0]*144))
        result = self.asre_dll.KBern3D_foot_TIM_dNA_interface(E, d_foot, b_foot, 
                                                              dx, EGratio, ni_str, 
                                                              d_NA, result_array)
        if result == 0:
            element_stiffness = np.array(list(result_array)).copy()
        else:
            raise RuntimeError(f'KBern3D_foot_TIM_dNA_interface failed with error code {result}')
        return element_stiffness

    def set_beam_properties(self, Eb, EoverG, q_foot, d_NA = 0):
        """
        Set the beam material properties.

        Parameters
        ----------
        Eb : float
            The Young's modulus of the beam. Unit: N/m^2
        EoverG : float
            The ratio between Young's modulus and the shear modulus. Unit: unitless
        q_foot : float
            The uniform weight applied on the beam. Unit: N/m (Weight along unit
            length in the longitudinal direction)
        """
        self.Eb = Eb
        self.EoverG = EoverG
        self.q_foot = q_foot
        self.ni_foot = self.EoverG/2 - 1 #The Poisson's ratio of beam. Unit: unitless
        self.d_NA = d_NA

    def set_soil_properties(self, EsNominal, nis, mu_int):
        """
        Set the beam material properties.

        Parameters
        ----------
        EsNominal : float
            The nominal elastic modulus of soil. Unit: N/m^2
        nis : float
            The Poisson's ratio of soil. Unit: unitless 
        mu_int : float
            The friction coefficient between soil and beam. Unit: unitless
        """
        self.EsNominal = EsNominal
        self.nis = nis
        self.mu_int = mu_int
        
    def run_model(self, dispL, dispT, dispV, output = 'disp'):
        """
        Run the SSI model under greenfield displacements

        Parameters
        ----------
        dispL : float
            The nominal elastic modulus of soil. Unit: N/m^2
        dispT : float
            The Poisson's ratio of soil. Unit: unitless 
        dispV : float
            The friction coefficient between soil and beam. Unit: unitless
        output : str
            If 'disp' then save beam displacement to self.beam_disp
            self.beam_dispL is the beam disp in longitudinal direction
            self.beam_dispT is the beam disp in transverse direction
            self.beam_dispV is the beam disp in vertical direction
            If 'strain' then save beam principal strain to self.strain
            Each element in self.strain is one principal strain of the beam
            If 'strain+disp' then save both self.disp and self.strain
        
        Returns
        -------
        bool
            Return Ture if run success and False if unsuccess 
        """
        self.ouput = output
        if self.ouput == 'disp':
            result_size = self.nnode * 6
            self.result_array_ptr = (c_double * result_size)(*([0]*result_size))
        elif self.ouput == 'strain':
            result_size = 3 * (self.nnode-1) * 2
            self.result_array_ptr = (c_double * result_size)(*([0]*result_size))
        elif self.ouput == 'strain+disp':
            result_size = 3 * (self.nnode-1) * 2 + self.nnode * 6
            self.result_array_ptr = (c_double * result_size)(*([0]*result_size))
        elif self.ouput == 'strain+disp+force':
            result_size = 3 * (self.nnode-1) * 2 + self.nnode * 6 + 3 * (self.nnode-1) * 2
            self.result_array_ptr = (c_double * result_size)(*([0]*result_size))
        else:
            raise ValueError(f'output value {output} is not permitted in ASRE_Timoshenko_model')
        
        self.ouput = self.ouput.encode('utf-8')
        self.solver = self.solver.encode('utf-8')
        try:
            result = self.asre_dll.run(self.nnode, self.meshX, self.meshY, self.meshZ,
                                       dispV, dispL, dispT, self.Eb, self.EoverG,
                                       self.EsNominal, self.nis, self.dfoot,
                                       self.bfoot, self.ni_foot, self.mu_int,
                                       self.q_foot, self.d_NA, self.solver, self.ouput,
                                       self.result_array_ptr, result_size)
            self.result = result
            self.result_array_ptr = np.array(list(self.result_array_ptr))
            if self.ouput.decode('utf-8') == 'disp':
                # The cpp lib returns the disp at the soil-beam interface (ground level)
                self.beam_DispL = self.result_array_ptr[0::6]
                self.beam_DispT = self.result_array_ptr[1::6]
                self.beam_DispV = self.result_array_ptr[2::6]
                self.beam_RotaL = self.result_array_ptr[3::6]
                self.beam_RotaT = self.result_array_ptr[4::6]
                self.beam_RotaV = self.result_array_ptr[5::6]
            elif self.ouput.decode('utf-8') == 'strain':
                self.beam_strain_top = self.result_array_ptr[0:(self.nnode-1)*2]
                self.beam_strain_bottom = self.result_array_ptr[(self.nnode-1)*2:(self.nnode-1)*4]
                self.beam_strain_diagonal = self.result_array_ptr[(self.nnode-1)*4:(self.nnode-1)*6]
            elif self.ouput.decode('utf-8') == 'strain+disp':
                self.beam_strain_top = self.result_array_ptr[0:(self.nnode-1)*2]
                self.beam_strain_bottom = self.result_array_ptr[(self.nnode-1)*2:(self.nnode-1)*4]
                self.beam_strain_diagonal = self.result_array_ptr[(self.nnode-1)*4:(self.nnode-1)*6]
                result_array = self.result_array_ptr[(self.nnode-1)*6:]
                self.beam_DispL = result_array[0::6]
                self.beam_DispT = result_array[1::6]
                self.beam_DispV = result_array[2::6]
                self.beam_RotaL = result_array[3::6]
                self.beam_RotaT = result_array[4::6]
                self.beam_RotaV = result_array[5::6]
            elif self.ouput.decode('utf-8') == 'strain+disp+force':
                self.beam_strain_top = self.result_array_ptr[0:(self.nnode-1)*2]
                self.beam_strain_bottom = self.result_array_ptr[(self.nnode-1)*2:(self.nnode-1)*4]
                self.beam_strain_diagonal = self.result_array_ptr[(self.nnode-1)*4:(self.nnode-1)*6]
                result_array = self.result_array_ptr[(self.nnode-1)*6:(self.nnode-1)*6+self.nnode*6]
                self.beam_DispL = result_array[0::6]
                self.beam_DispT = result_array[1::6]
                self.beam_DispV = result_array[2::6]
                self.beam_RotaL = result_array[3::6]
                self.beam_RotaT = result_array[4::6]
                self.beam_RotaV = result_array[5::6]
                result_array = self.result_array_ptr[(self.nnode-1)*6+self.nnode*6:]
                self.moment = result_array[0:(self.nnode-1)*2]
                self.axialForce = result_array[(self.nnode-1)*2:(self.nnode-1)*4]
                self.shearForce = result_array[(self.nnode-1)*4:(self.nnode-1)*6]
            # Get the disp at the beam-axis
            _ = self.get_beam_axis_disp()
            self.compute_internal_forces()
            self.compute_tensile_strain()
            
            return True
        except:
            self.release_cdll_handle()
            raise RuntimeError(f'ASRE_Timoshenko_model failed to run the ASRE cpp library')
    
    def compute_internal_forces(self):

        # Number of elements
        num_elements = self.nnode - 1

        # Initialize the internal forces vectors
        F_M_deltaT_el = np.zeros(2 * num_elements)
        F_N_deltaT_el = np.zeros(2 * num_elements)
        F_S_deltaT_el = np.zeros(2 * num_elements)

        # Compute the internal forces for each element
        for j in range(num_elements):
            # use the stiffness and disp at element axis
            K_local = self.KBern3D_foot_TIM_dNA(elem_i=j, d_NA=0).reshape((12, 12))  # Get the local stiffness matrix for the current element
            # Extract the nodal displacements for the current element
            element_disp = np.zeros((12, 1))
            element_disp[0:6] = self.axis_disp[j * 6:(j + 1) * 6]
            if (j + 1) * 6 < len(self.axis_disp):
                element_disp[6:12] = self.axis_disp[(j + 1) * 6:(j + 2) * 6]

            # Compute the internal forces for the current element
            Fin_P = np.dot(K_local, element_disp)

            # Store the internal forces in the vectors
            F_M_deltaT_el[2 * j] = -Fin_P[4, 0]
            F_M_deltaT_el[2 * j + 1] = Fin_P[10, 0]
            F_N_deltaT_el[2 * j] = -Fin_P[0, 0]
            F_N_deltaT_el[2 * j + 1] = Fin_P[6, 0]
            F_S_deltaT_el[2 * j] = -Fin_P[2, 0]
            F_S_deltaT_el[2 * j + 1] = Fin_P[8, 0]

        self.moment = F_M_deltaT_el
        self.axialForce = F_N_deltaT_el
        self.shearForce = F_S_deltaT_el
        return
    
    def compute_tensile_strain(self):
        Ab = self.bfoot * self.dfoot  # Cross-sectional area of the beam
        poissons_ratio = self.ni_foot  # Poisson's ratio of the beam
        Gb = self.Eb / (2 * (1 + poissons_ratio))  # Shear modulus of the beam
        As = Ab * (10 + (poissons_ratio * 10)) / (12 + (11 * poissons_ratio))  # According to Wikipedia
        GAs = Gb * As
        shear_factor_midpoint = 1.5 # for Timoshenko beam theory, shear factor 
        EA = self.Eb * Ab  # Axial stiffness of the beam
        EI = self.Eb * (self.bfoot * self.dfoot**3) / 12  # Bending stiffness of the beam
        # N / EA
        strain_axial_normal = self.axialForce / EA
        # ( M / EI ) * d
        strain_axial_bending_top = (self.moment / EI) * (self.dfoot - self.d_NA)  # Compressive strains in the top
        strain_axial_bending_bottom = (self.moment / EI) * self.d_NA

        # Navier's formula
        tensile_strain_top = strain_axial_normal + strain_axial_bending_top
        tensile_strain_bottom = strain_axial_normal + strain_axial_bending_bottom

        true_shear_strain = (self.shearForce * shear_factor_midpoint / GAs) / 2

        # Bending's contribution e_xx * (1 - v) / 2
        strain_exx_contribution = ((strain_axial_normal * (1 - poissons_ratio)) / 2)
        # Formula 7 from slides
        tensile_strain_midpoint = strain_exx_contribution + np.sqrt(
            (((strain_axial_normal * (1 + poissons_ratio)) / 2) ** 2)
            + (true_shear_strain ** 2))

        eps_t = np.zeros_like(tensile_strain_midpoint)  # Create return vector
        for i in range(len(tensile_strain_midpoint.flatten())):
            eps_t[i] = max(tensile_strain_midpoint[i], tensile_strain_top[i], tensile_strain_bottom[i])
        self.eps_t = eps_t
        self.tensile_strain_top = tensile_strain_top
        self.tensile_strain_bottom = tensile_strain_bottom
        self.tensile_strain_midpoint = tensile_strain_midpoint
        return

    def get_beam_axis_disp(self):
        """
        Get the beam axis displacement.

        Returns
        -------
        np.array(dtype=float)
            The beam axis displacement in the unit of m.
        """
        if hasattr(self, 'beam_DispL'):
            axis_DispL = self.beam_DispL + self.beam_RotaT * self.d_NA
            numNodes = self.nnode
            total_disp = np.zeros([numNodes * 6, 1])
            total_disp[0::6] = axis_DispL.reshape(-1, 1)  # [m] Longitudinal / axial (x)
            total_disp[1::6] = self.beam_DispT.reshape(-1, 1)  # [m] Transversal (y)
            total_disp[2::6] = self.beam_DispV.reshape(-1, 1)  # [m] Vertical (z)
            total_disp[3::6] = self.beam_RotaL.reshape(-1, 1)  # [rad] Horizontal rotations
            total_disp[4::6] = self.beam_RotaT.reshape(-1, 1)  # [rad] Transversal rotations
            total_disp[5::6] = self.beam_RotaV.reshape(-1, 1)  # [ras] Vertical rotations
            self.axis_disp = total_disp
            return axis_DispL
        else:
            raise AttributeError('The beam axis displacement is not computed yet. Please run the model first.')

    def categorize_damage(self):
        """
        Categorize the damage based on the computed tensile strain.

        Returns:
        damage_category : str - The damage category based on the tensile strain based
        on Boscardin and Cording, 1989
        """

        tensile_strain = self.eps_t.max() * 100  # Convert to percent
        # print("Maximum tensile strain is %", tensile_strain)
        if tensile_strain < 0.05:
            return 0, tensile_strain
        elif tensile_strain < 0.075:
            return 1, tensile_strain
        elif tensile_strain < 0.15:
            return 2, tensile_strain
        elif tensile_strain < 0.3:
            return 3, tensile_strain
        elif tensile_strain < 0.6:
            return 4, tensile_strain
        else:
            return 5, tensile_strain
    
    def _write_input_file(self,  dispL, dispT, dispV, output = 'disp',
                          folder=None, filename=None):
        if filename is None:
            filename = 'input.json'
        if folder is None:
            folder = os.getcwd()
        input_json = {
            'nnode': self.nnode,
            'meshX': self.meshX.tolist(),
            'meshY': self.meshY.tolist(),
            'meshZ': self.meshZ.tolist(),
            'dispL': dispL.tolist(),
            'dispT': dispT.tolist(),
            'dispV': dispV.tolist(),
            'Eb': self.Eb,
            'EoverG': self.EoverG,
            'EsNominal': self.EsNominal,
            'nis': self.nis,
            'dfoot': self.dfoot,
            'bfoot': self.bfoot,
            'ni_foot': self.ni_foot,
            'd_NA': self.d_NA,
            'mu_int': self.mu_int,
            'q_foot': self.q_foot,
            'solver': self.solver,
            'output': output
        }
        with open(os.path.join(folder, filename), 'w') as f:
            json.dump(input_json, f, indent=2)
        return os.path.join(folder, filename)
    
    def release_cdll_handle(self):
        if self.asre_dll is None:
            pass
        else:
            lib_handle = self.asre_dll._handle
            if pltm == "win32": 
                del self.asre_dll
                _ctypes.FreeLibrary(lib_handle)
                self.asre_dll = None
            elif pltm == "darwin":
                dl = CDLL('libdl.dylib')
                dl.dlclose(lib_handle)

 
