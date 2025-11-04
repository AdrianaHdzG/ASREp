from sys import platform as pltm
import os, platform, pkg_resources
from ctypes import CDLL, c_int, c_double, c_char_p, POINTER
import numpy as np
import json
import _ctypes

class ASRE_3D_solid_model:
    def __init__(self, building_node_coords, building_elements,
                  interface_node_coords, interface_elements,
                 solver = 'elasto-plastic'):
        if building_elements.dtype != np.int32:
            building_elements = building_elements.astype(np.int32)
        if interface_elements.dtype != np.int32:
            interface_elements = interface_elements.astype(np.int32)
        self.building_node_coords = building_node_coords
        self.building_elements = building_elements
        self.interface_node_coords = interface_node_coords
        self.interface_elements = interface_elements
        # self.stiffness_matrix = self._create_building_stiffness(
        #     building_node_coords, building_elements,
        #     Eb, nu, rho
        # )
        self.solver = solver
        self.asre_dll = self._import_dll()
        # Create a temp dir
        # self.run_dir = run_dir
        # if not os.path.exists(self.run_dir):
        #     raise ValueError(f'run_dir {self.run_dir} does not exist')
        # try:
        #     self.temp_dir = os.path.join(self.run_dir, 'temp')
        #     os.mkdir(self.temp_dir, exist_ok=True)
        # except:
        #     raise RuntimeError(f'Failed to create temp directory in {self.run_dir}')

    def set_building_properties(self, Eb, nu, rho):
        """
        Set the material properties of the building.
        Parameters
        ----------
        Eb : float
            The elastic modulus of building. Unit: N/m^2
        nu : float
            The Poisson's ratio of building. Unit: unitless
        rho : float
            The unit weight of building. Unit: N/m^3
        """
        self.Eb = Eb
        self.nu = nu
        self.rho = rho

    # Define __getstate__ and __setstate__ methods to handle unpicklable attributes
    def __getstate__(self):
        state = self.__dict__.copy()
        # Remove the unpicklable attribute
        state['asre_dll'] = None
        return state
    def __setstate__(self, state):
        self.__dict__.update(state)
        # Reinitialize the unpicklable attribute
        self.asre_dll = self._import_dll()
    
    # def set_building_properties(self, Eb, nu, rho):
    #     original_Eb = self.Eb
    #     original_nu = self.nu
    #     self.Eb = Eb
    #     self.nu = nu
    #     self.rho = rho
    #     if (original_nu != self.nu):
    #         self._create_building_stiffness(
    #             self.building_node_coords, self.building_elements,
    #             Eb, nu, rho
    #         )
    #     else:
    #         # the stiffness matrix can be simply scaled
    #         scale_factor = self.Eb / original_Eb
    #         self.stiffness_matrix = scale_factor * self.stiffness_matrix
            

    # def _create_building_stiffness(self, building_node_coords, 
    #                                building_elements, Eb, nu, rho):
    #     # Save the building stiffness matrix and self weight to files in the temp dir
    #     try: 
    #         self.asre_dll.computeKKfootAndBodyForce(
    #             self.building_node_coords,
    #             self.building_elements,
    #             Eb,
    #             nu,
    #             rho,
    #             self.temp_dir
    #         )
    #     except:
    #         self.release_cdll_handle()
    #         raise RuntimeError('ASRE_3D_solid_model failed to create building stiffness matrix using the ASRE cpp library')
        # Load the stiffness matrix and self weight from files in the temp dir


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
                                    "ASREcpp", "bin", "macOS_arm", "libASRElib3DSolid.so")
            if os.path.exists(lib_path):
                c_lib = CDLL(lib_path)
            else:
                c_lib = None
                message = f'ASRE is not precompiled for {pltm}, please compile the ASRE cpp library'
        elif pltm == "darwin":
            if platform.processor() == 'arm':
                lib_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "ASREcpp", "bin", "macOS_arm", "libASRElib3DSolid.dylib")
                if os.path.exists(lib_path):
                    c_lib = CDLL(lib_path)
                else:
                    c_lib = None
                    message = f'ASRE is not precompiled for {pltm} {platform.processor()}, please compile the ASRE cpp library'
            else:
                lib_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "ASREcpp", "bin", "macOS", "libASRElib3DSolid.dylib")
                if os.path.exists(lib_path):
                    c_lib = CDLL(lib_path)
                else:
                    c_lib = None
                    message = f'ASRE is not precompiled for {pltm} {platform.processor()}, please compile the ASRE cpp library'
        elif pltm == "win32":
            lib_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "ASREcpp", "bin", "win32", "ASRElib3DSolid.dll")
            print("[ASRE] lib_path:", lib_path)
            print("[ASRE] exists? ", os.path.exists(lib_path))
                                    
            if os.path.exists(lib_path):
                c_lib = CDLL(lib_path)
            else:
                c_lib = None
                message = f'ASRE is not precompiled for {pltm}, please compile the ASRE cpp library'
            # c_lib = CDLL(pkg_resources.resource_filename('ASREp', 'ASREcpp//bin//win32//ASRElib.dll'))
            # c_lib.printName()
        if c_lib is None:
            raise ImportError(message)
        # Define the argument and return types of the run function
        c_lib.run.argtypes = [
            c_int, # nnode
            np.ctypeslib.ndpointer(dtype=np.float64), # meshX
            np.ctypeslib.ndpointer(dtype=np.float64), # meshY
            np.ctypeslib.ndpointer(dtype=np.float64), # meshZ
            c_int, #nelem
            np.ctypeslib.ndpointer(dtype=np.int32), # elemNode1
            np.ctypeslib.ndpointer(dtype=np.int32), # elemNode2
            np.ctypeslib.ndpointer(dtype=np.int32), # elemNode3
            np.ctypeslib.ndpointer(dtype=np.int32), # elemNode4
            np.ctypeslib.ndpointer(dtype=np.int32), # elemNode5
            np.ctypeslib.ndpointer(dtype=np.int32), # elemNode6
            np.ctypeslib.ndpointer(dtype=np.int32), # elemNode7
            np.ctypeslib.ndpointer(dtype=np.int32), # elemNode8
            c_int, # nnode_inter
            np.ctypeslib.ndpointer(dtype=np.float64), # interNodeX
            np.ctypeslib.ndpointer(dtype=np.float64), # interNodeY
            np.ctypeslib.ndpointer(dtype=np.float64), # interNodeZ
            c_int, # nelem_inter
            np.ctypeslib.ndpointer(dtype=np.int32), # interElemNode1
            np.ctypeslib.ndpointer(dtype=np.int32), # interElemNode2
            np.ctypeslib.ndpointer(dtype=np.int32), # interElemNode3
            np.ctypeslib.ndpointer(dtype=np.int32), # interElemNode4
            c_double, # buildingE
            c_double, # buildingNu
            c_double, # buildingRho
            c_double, # Gs
            c_double, # nus
            c_double, # mu_int
            np.ctypeslib.ndpointer(dtype=np.float64), # dispX
            np.ctypeslib.ndpointer(dtype=np.float64), # dispY
            np.ctypeslib.ndpointer(dtype=np.float64), # dispZ
            c_double, # lim_t_int
            c_double, # lim_c_int
            c_char_p, # solver
            c_int, # print_iteration
            POINTER(c_double), # result array
        ]
        c_lib.run.restype = c_int
        # Define the argument and return types of the calculateStrain function
        c_lib.calculateStrain.argtypes = [
            c_int, # nnode
            np.ctypeslib.ndpointer(dtype=np.float64), # dispX
            np.ctypeslib.ndpointer(dtype=np.float64), # dispY
            np.ctypeslib.ndpointer(dtype=np.float64), # dispZ
            c_int, # nelem
            np.ctypeslib.ndpointer(dtype=np.int32), # elemNode1
            np.ctypeslib.ndpointer(dtype=np.int32), # elemNode2
            np.ctypeslib.ndpointer(dtype=np.int32), # elemNode3
            np.ctypeslib.ndpointer(dtype=np.int32), # elemNode4
            np.ctypeslib.ndpointer(dtype=np.int32), # elemNode5
            np.ctypeslib.ndpointer(dtype=np.int32), # elemNode6
            np.ctypeslib.ndpointer(dtype=np.int32), # elemNode7
            np.ctypeslib.ndpointer(dtype=np.int32), # elemNode8
            np.ctypeslib.ndpointer(dtype=np.float64), # displacement
            POINTER(c_double), # principal tensile strain
            POINTER(c_double), # principal compressive strain
        ]
        c_lib.calculateStrain.restype = c_int

        return c_lib


    def set_soil_properties(self, EsNominal, nus):
        """
        Set the beam material properties.

        Parameters
        ----------
        EsNominal : float
            The nominal elastic modulus of soil. Unit: N/m^2
        nus : float
            The Poisson's ratio of soil. Unit: unitless
        """
        self.EsNominal = EsNominal
        self.nus = nus
        self.Gs = EsNominal / (2*(1+nus))
    
    def set_interface_properties(self, mu_int, lim_t_int, lim_c_int):
        """
        Set the interface properties between soil and beam.

        Parameters
        ----------
        mu_int : float
            The friction coefficient between soil and beam. Unit: unitless
        lim_t_int : float
            The shear strength limit of the interface. Unit: N/m^2
        lim_c_int : float
            The compressive strength limit of the interface. Unit: N/m^2
        """
        self.mu_int = mu_int
        self.lim_t_int = lim_t_int
        self.lim_c_int = lim_c_int

    def calculate_principal_strain(self):
        """
        Calculate the principal tensile and compressive strains of the building elements.
        Returns
        -------
        bool
            Return True if calculation is successful, False otherwise.
        """
        if not hasattr(self, 'result_array_ptr'):
            raise RuntimeError('No displacement results found. Please run the model first.')
        result_size = self.building_elements.shape[0]
        self.result_tensile_ptr = (c_double * result_size)(*([0]*result_size))
        self.result_compressive_ptr = (c_double * result_size)(*([0]*result_size))
        try:
            result = self.asre_dll.calculateStrain(
                self.building_node_coords.shape[0],
                self.building_node_coords[:,0],
                self.building_node_coords[:,1],
                self.building_node_coords[:,2],
                self.building_elements.shape[0],
                self.building_elements[:,0],
                self.building_elements[:,1],
                self.building_elements[:,2],
                self.building_elements[:,3],
                self.building_elements[:,4],
                self.building_elements[:,5],
                self.building_elements[:,6],
                self.building_elements[:,7],
                self.result_array_ptr,
                self.result_tensile_ptr,
                self.result_compressive_ptr
            )
            self.result_tensile_ptr = np.array(list(self.result_tensile_ptr))
            self.result_compressive_ptr = np.array(list(self.result_compressive_ptr))
            return True
        except:
            self.release_cdll_handle()
            raise RuntimeError('ASRE_3D_solid_model failed to calculate principal strains using the ASRE cpp library')

    def run_model(self, dispX, dispY, dispZ, print_iteration = 0):
        """
        Run the SSI model under greenfield displacements

        Parameters
        ----------
        dispX : float
            The greenfield displacement in X direction. Unit: m
        dispY : float
            The greenfield displacement in Y direction. Unit: m
        dispZ : float
            The greenfield displacement in Z direction. Unit: m
        print_iteration : int, optional
            If print_iteration == 1, print the iteration information during the model run.
            Default is 0, no print.
        Returns
        -------
        bool
            Return True if run success and False if unsuccess
            self.building_dispX is the beam disp in longitudinal direction
            self.building_dispY is the beam disp in transverse direction
            self.building_dispZ is the beam disp in vertical direction
        """
        result_size = self.building_node_coords.shape[0] * 3
        self.result_size = result_size
        self.result_array_ptr = (c_double * result_size)(*([0]*result_size))
        if isinstance(self.solver, str):
            self.solver = self.solver.encode('utf-8')
        try:
            result = self.asre_dll.run(
                self.building_node_coords.shape[0],
                self.building_node_coords[:,0],
                self.building_node_coords[:,1],
                self.building_node_coords[:,2],
                self.building_elements.shape[0],
                self.building_elements[:,0],
                self.building_elements[:,1],
                self.building_elements[:,2],
                self.building_elements[:,3],
                self.building_elements[:,4],
                self.building_elements[:,5],
                self.building_elements[:,6],
                self.building_elements[:,7],
                self.interface_node_coords.shape[0],
                self.interface_node_coords[:,0],
                self.interface_node_coords[:,1],
                self.interface_node_coords[:,2],
                self.interface_elements.shape[0],
                self.interface_elements[:,0],
                self.interface_elements[:,1],
                self.interface_elements[:,2],
                self.interface_elements[:,3],
                self.Eb,
                self.nu,
                self.rho,
                self.Gs,
                self.nus,
                self.mu_int,
                dispX,
                dispY,
                dispZ,
                self.lim_t_int,
                self.lim_c_int,
                self.solver,
                print_iteration,
                self.result_array_ptr
            )

            self.result = result
            self.result_array_ptr = np.array(list(self.result_array_ptr))
            self.building_dispX = self.result_array_ptr[0::3]
            self.building_dispY = self.result_array_ptr[1::3]
            self.building_dispZ = self.result_array_ptr[2::3]
            return True
        except:
            self.release_cdll_handle()
            raise RuntimeError('ASRE_Timoshenko_model failed to run the ASRE cpp library')
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
 
