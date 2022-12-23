from copy import deepcopy

from Algorithms.DesignToolAlgorithmV3_1D.EXTRAS.joinBlocks import jointBlock
from Algorithms.DesignToolAlgorithmV3_1D.BoundaryConditions.FormBoundaryConditionInformation import FormBoundaryConditionInformation
class Integrate():
    def __init__(self, mesh, cfl_flag, tCurrent, currentStep) -> None:
        """
        cfl_flag = [Bool, float] -> Bool = True if want cfl criteria to be used, false if fixed dt
                                 -> float = cfl_max if Bool = True, dt_max if Bool = False
        """
        self.mesh = mesh
        
        ### Want to make copies of mapping arrays because it's inefficient to remake them after removing the ghost cells at the end of the time step
        oldMapCellIDToWestInterfaceIdx = deepcopy(mesh.mapCellIDToWestInterfaceIdx)
        oldMapCellIDToEastInterfaceIdx = deepcopy(mesh.mapCellIDToEastInterfaceIdx)
        oldMapInterfaceIDToWestCellIdx = deepcopy(mesh.mapInterfaceIDToWestCellIdx)
        oldMapInterfaceIDToEastCellIdx = deepcopy(mesh.mapInterfaceIDToEastCellIdx)
        boundaryConditions = deepcopy(mesh.boundaryConditions)
        boundaryInterfaceIDs = deepcopy(mesh.boundaryInterfaceIDs)
        ### Add boundary conditions
        #for cell in mesh.cellArray:
            #print(cell.cell_ID, cell.conservedProperties["mass"], cell.conservedProperties["xMom"], cell.conservedProperties["energy"])
        self.addBoundaryConditions(BCs = self.mesh.boundaryConditions)

        #print("Printing cell properties at start of time step")
        #for cell in self.mesh.cellArray:
            #print(cell.cell_ID, cell.fs.fluid_state.p, cell.fs.fluid_state.T, cell.fs.fluid_state.rho, cell.fs.vel_x)
        #for cell in self.mesh.cellArray:
            #print(cell.cell_ID, cell.fs.fluid_state.p, cell.fs.fluid_state.T, cell.fs.fluid_state.rho, cell.fs.vel_x)
        #print(self.mesh.mapCellIDToWestInterfaceIdx)
        #print(self.mesh.mapCellIDToEastInterfaceIdx)
        #print(self.mesh.mapInterfaceIDToWestCellIdx)
        #print(self.mesh.mapInterfaceIDToEastCellIdx)

        ### Find out inviscid flux dt from CFL
        if not cfl_flag[0]: #Using fixed dt, cfl_flag has form [False, float dt_fixed]
            self.dtTotal = cfl_flag[1]
        else:
            if cfl_flag[1] == "no_ramping": 
                # Using fixed cfl
                # cfl_flag has form [True, "no_ramping", float cfl_fixed]
                self.dtTotal = 1e6 #Large number to be initialised, won't be used
                for cell in self.mesh.cellArray:
                    self.dtTotal = min(self.dtTotal, cell.maxAllowableDt(cfl = cfl_flag[2]))
                
            elif cfl_flag[1] == "cfl_ramp_from_tCurrent": 
                # Ramping cfl based on tCurrent
                # cfl_flag has form [True, "cfl_ramp_from_tCurrent", [[cfl_1, t_1], [cfl_2, t_2], ..., [cfl_n, t_n]]]
                # If we are outside of t_n, use cfl_n. Else, find first instance tCurrent < t_i, then use cfl_i
                if tCurrent > cfl_flag[2][-1][1]: 
                    cfl_max = cfl_flag[2][-1][0]
                    #print("cfl value used: ", cfl_max)
                    self.dtTotal = 1e6 #Large number to be initialised, won't be used
                    for cell in self.mesh.cellArray:
                        self.dtTotal = min(self.dtTotal, cell.maxAllowableDt(cfl = cfl_max))
                else:
                    for cfl_pair in range(len(cfl_flag[2])):
                        if tCurrent > cfl_flag[2][cfl_pair][1]:
                            pass
                        else:
                            cfl_current = cfl_flag[2][cfl_pair][0]
                            #print("cfl value used: ", cfl_current)
                            self.dtTotal = 1e6 #Large number to be initialised, won't be used
                            for cell in self.mesh.cellArray:
                                self.dtTotal = min(self.dtTotal, cell.maxAllowableDt(cfl = cfl_current))
                            break
                
            elif cfl_flag[1] == "cfl_ramp_from_step":
                # Ramping cfl based on currentStep
                # cfl_flag has form [True, "cfl_ramp_from_step", [[cfl_1, step_1], [cfl_2, step_2], ..., [cfl_n, step_n]]]
                # If we are outside of step_n, use cfl_n. Else, find first instance currentStep < step_i, then use cfl_i
                if currentStep > cfl_flag[2][-1][1]: 
                    cfl_max = cfl_flag[2][-1][0]
                    #print("step: ", currentStep, " cfl used: ", cfl_max)
                    self.dtTotal = 1e6 #Large number to be initialised, won't be used
                    for cell in self.mesh.cellArray:
                        self.dtTotal = min(self.dtTotal, cell.maxAllowableDt(cfl = cfl_max))

                else:
                    for cfl_pair in range(len(cfl_flag[2])):
                        if currentStep > cfl_flag[2][cfl_pair][1]:
                            pass
                        
                        else:
                            cfl_current = cfl_flag[2][cfl_pair][0]
                            #print("step: ", currentStep, " cfl used: ", cfl_current)
                            self.dtTotal = 1e6 #Large number to be initialised, won't be used
                            for cell in self.mesh.cellArray:
                                self.dtTotal = min(self.dtTotal, cell.maxAllowableDt(cfl = cfl_current))
                            break
                
            elif cfl_flag[1] == "dt_ramp_then_cfl_ramp_from_step":
                # Ramping dt initially, then use cfl based on current step number
                # cfl_flag has the form [True, "dt_ramp_then_cfl_ramp_from_step", [dt_init, scale], [[cfl_1, step_1], cfl_2, step_2], ..., [cfl_n, step_n]]]
                # Starting at dt = dt_init, ramp this every time step by scale until the time step is greater than what comes from
                # cfl_1. Then ramp cfl similarly to "cfl_ramp_from_step" method.
                [dt_init, scale] = cfl_flag[2]
                dt_ramped = dt_init * scale ** currentStep
                min_dt_from_cfl = 1e6
                for cell in self.mesh.cellArray:
                    min_dt_from_cfl = min(min_dt_from_cfl, cell.maxAllowableDt(cfl = cfl_flag[3][0][0]))
                if dt_ramped >= min_dt_from_cfl:
                    if currentStep > cfl_flag[3][-1][1]:
                        self.dtTotal = 1e6
                        for cell in self.mesh.cellArray:
                            self.dtTotal = min(self.dtTotal, cell.maxAllowableDt(cfl = cfl_flag[3][-1][0]))
                        print("step = ", currentStep, " final cfl criteria used, cfl = ", cfl_flag[3][-1][0])
                    else:
                        for cfl_pair in range(len(cfl_flag[3])):
                            if currentStep > cfl_flag[3][cfl_pair][1]:
                                pass
                            else:
                                cfl_current = cfl_flag[3][cfl_pair][0]
                                print("step = ", currentStep, " Intermediate cfl value used = ", cfl_current)
                                self.dtTotal = 1e6 #Large number to be initialised, won't be used
                                for cell in self.mesh.cellArray:
                                    self.dtTotal = min(self.dtTotal, cell.maxAllowableDt(cfl = cfl_current))
                                break
                else:
                    print("step = ", currentStep, " ramped dt used: ", dt_ramped, " cfl dt: ", min_dt_from_cfl)
                    self.dtTotal = dt_ramped

            elif cfl_flag[1] == "dt_ramp_then_cfl_ramp_from_tCurrent":
                # Ramping dt initially, then use cfl based on current time.
                # cfl_flag has the form [True, "dt_ramp_then_cfl_ramp_from_tCurrent", [dt_init, scale], [[cfl_1, t_1], [cfl_2, t_2], ..., [cfl_n, t_n]]]
                # Starting at dt = dt_init, ramp this every time step by scale until the time step is greater than what comes from
                # cfl_1. Then ramp cfl similarly to "cfl_ramp_from_tCurrent".
                [dt_init, scale] = cfl_flag[2]
                dt_ramped = dt_init * scale ** currentStep
                min_dt_from_cfl = 1e6
                for cell in self.mesh.cellArray:
                    min_dt_from_cfl = min(min_dt_from_cfl, cell.maxAllowableDt(cfl = cfl_flag[3][0][0]))
                if dt_ramped >= min_dt_from_cfl:
                    if tCurrent > cfl_flag[3][-1][1]:
                        self.dtTotal = 1e6
                        for cell in self.mesh.cellArray:
                            self.dtTotal = min(self.dtTotal, cell.maxAllowableDt(cfl = cfl_flag[3][-1][0]))
                    else:
                        for cfl_pair in range(len(cfl_flag[3])):
                            if tCurrent > cfl_flag[3][cfl_pair][1]:
                                pass
                            else:
                                cfl_current = cfl_flag[3][cfl_pair][0]
                                self.dtTotal = 1e6 #Large number to be initialised, won't be used
                                for cell in self.mesh.cellArray:
                                    self.dtTotal = min(self.dtTotal, cell.maxAllowableDt(cfl = cfl_current))
                                break
                else:
                    self.dtTotal = dt_ramped
            else:
                print("Invalid CFL criteria given.")
            
        #print("This is the pre-convection properties")
        #for ind, cell in enumerate(self.mesh.cellArray):
            #print("index: ", ind, "p", cell.fs.fluid_state.p, ", T", cell.fs.fluid_state.T, ", rho", cell.fs.fluid_state.rho, ", vel_x", cell.fs.vel_x, "a", cell.fs.fluid_state.a)
            #print("index: ", ind, "Mass:", cell.conservedProperties["mass"], "xMom:", cell.conservedProperties["xMom"], "Energy:", cell.conservedProperties["energy"])
        ### Perform iteration over interfaces and update conserved properties

        for interface in self.mesh.interfaceArray:
            ### Do all the interface specific methods that are defined within each interface (reconstruction, fluxes, perform conserved property updates etc etc)
            interface.completeInterfaceMethods(cellArray = self.mesh.cellArray, \
                                                interfaceArray = self.mesh.interfaceArray, \
                                                mapCellIDToWestInterfaceIdx = self.mesh.mapCellIDToWestInterfaceIdx, \
                                                mapCellIDToEastInterfaceIdx = self.mesh.mapCellIDToEastInterfaceIdx, \
                                                mapInterfaceIDToWestCellIdx = self.mesh.mapInterfaceIDToWestCellIdx, \
                                                mapInterfaceIDToEastCellIdx = self.mesh.mapInterfaceIDToEastCellIdx, \
                                                dt_inv = self.dtTotal)
        ### Update properties after doing interface methods. Don't have to do this for boundary cells but we will for now
        #print("This is the post-convection properties")
        for ind, cell in enumerate(self.mesh.cellArray):
            cell.update_properties()
            #print("index: ", ind, "p", cell.fs.fluid_state.p, ", T", cell.fs.fluid_state.T, ", rho", cell.fs.fluid_state.rho, ", vel_x", cell.fs.vel_x, "a", cell.fs.fluid_state.a)
            #print("index: ", ind, "Mass:", cell.conservedProperties["mass"], "xMom:", cell.conservedProperties["xMom"], "Energy:", cell.conservedProperties["energy"])
        

        ### Perform iteration over cells and update conserved properties
        for cell in self.mesh.cellArray:
            cell.completeCellMethods(cellArray = self.mesh.cellArray, dt_inv = self.dtTotal)
        
        for ind, cell in enumerate(self.mesh.cellArray):
            cell.update_properties()
        
        #print("Printing cell properties at end of time step")
        #for cell in self.mesh.cellArray:
            #print(cell.cell_ID, cell.fs.fluid_state.p, cell.fs.fluid_state.T, cell.fs.fluid_state.rho, cell.fs.vel_x)
        
        ### Remove boundary cells and reset mapping arrays to what they were before adding boundary conditions
        self.mesh.mapCellIDToWestInterfaceIdx = oldMapCellIDToWestInterfaceIdx
        self.mesh.mapCellIDToEastInterfaceIdx = oldMapCellIDToEastInterfaceIdx
        self.mesh.mapInterfaceIDToWestCellIdx = oldMapInterfaceIDToWestCellIdx
        self.mesh.mapInterfaceIDToEastCellIdx = oldMapInterfaceIDToEastCellIdx
        self.mesh.boundaryConditions = boundaryConditions
        self.mesh.boundaryInterfaceIDs = boundaryInterfaceIDs
        self.mesh.cellArray = self.mesh.cellArray[:-1 * self.totalGhostCells]
        self.mesh.interfaceArray = self.mesh.interfaceArray[:-1 * self.totalGhostCells]

        #print("This is the end of time step conserved quantities")
        #for ind, cell in enumerate(self.mesh.cellArray):
            #print("index: ", ind, "p", cell.fs.fluid_state.p, ", T", cell.fs.fluid_state.T, ", rho", cell.fs.fluid_state.rho, ", vel_x", cell.fs.vel_x, "a", cell.fs.fluid_state.a)
            #print("index: ", ind, "Mass:", cell.conservedProperties["mass"], "xMom:", cell.conservedProperties["xMom"], "Energy:", cell.conservedProperties["energy"])

    def addBoundaryConditions(self, BCs):
        self.totalGhostCells = 0
        
        for boundaryCondition in BCs:
            self.totalGhostCells += boundaryCondition[1][1]
            if boundaryCondition[1][0] == "ConstantFlux_BC":
                pass
            else:
                ghostCellLayer = FormBoundaryConditionInformation(mesh = self.mesh, BC = boundaryCondition).ghostCellLayer
                self.mesh = jointBlock(meshObject1 = self.mesh, meshObject2 = ghostCellLayer, \
                    block1InterfaceIDBeingReplaced = boundaryCondition[0], block2InterfaceIDBeingReplaced = 0, \
                    newInterface = deepcopy(self.mesh.interfaceArray[boundaryCondition[0]]), AddingGhostCellsBool = True)
            