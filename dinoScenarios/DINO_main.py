import sys, os, inspect

bskName = 'Basilisk'
bskPath = '../..' + '/' + bskName + '/'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')

import SimulationBaseClass
import DINO_DKE
import DINO_FSW
import DINO_multiScenarios as scene
import sim_model


class DINO_DynSim(SimulationBaseClass.SimBaseClass):
    def __init__(self):
        # Create a sim module as an empty container
        SimulationBaseClass.SimBaseClass.__init__(self)
        # Create simulation process names
        self.DynamicsProcessName = "DynamicsProcess"
        self.FSWProcessName = "FSWProcess"
        self.DynPyProcessName = "DynPyProcess"

        # Create processes
        self.dynProc = self.CreateNewProcess(self.DynamicsProcessName)
        self.fswProc = self.CreateNewProcess(self.FSWProcessName)
        self.dynPyProc = self.CreateNewPythonProcess(self.DynPyProcessName)

        #   Create SysInterfaces for each process (jesus)
        self.dyn2dynPyInterface = sim_model.SysInterface()
        self.dyn2dynPyInterface.addNewInterface(
            self.DynPyProcessName, self.DynamicsProcessName)
        self.dynPyProc.addInterfaceRef(self.dyn2dynPyInterface)

        # self.dyn2fswInterface = sim_model.SysInterface()
        # self.dyn2pyFswInterface = sim_model.SysInterface()
        # self.fsw2pyFswInterface = sim_model.SysInterface()
        # self.dyn2fswInterface.addNewInterface(self.DynamicsProcessName,self.FSWProcessName)
        # self.dyn2pyFswInterface.addNewInterface(self.DynamicsProcessName, self.FSWPyProcessName)
        # self.fsw2pyFswInterface.addNewInterface(self.FSWPyProcessName, self.FSWPyProcessName)
        # self.dynProc.addInterfaceRef(self.dyn2pyFswInterface)
        # self.dynProc.addInterfaceRef(self.dyn2fswInterface)
        # self.fswProc.addInterfaceRef(self.dyn2fswInterface)
        # self.fswPyProc.addInterfaceRef(self.dyn2pyFswInterface)
        # self.fswPyProc.addInterfaceRef(self.fsw2pyFswInterface)

        # Crate sim subclasses
        self.DynClass = DINO_DKE.DynamicsClass(self)
        self.FSWClass = DINO_FSW.FSWClass(self)

if __name__ == "__main__":
    scene.opnavCamera_dynScenario(DINO_DynSim())
    # scene.attFilter_dynScenario(DINO_DynSim())
