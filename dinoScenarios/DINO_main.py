import sys, os, inspect

bskName = 'Basilisk'
bskPath = '../..' + '/' + bskName + '/'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')

try:
    import SimulationBaseClass
except ImportError:
    from Basilisk.utilities import SimulationBaseClass
import DINO_DKE
import DINO_FSW
import DINO_multiScenarios as scene


class DINO_DynSim(SimulationBaseClass.SimBaseClass):
    def __init__(self):
        # Create a sim module as an empty container
        SimulationBaseClass.SimBaseClass.__init__(self)
        # Create simulation process names
        self.DynamicsProcessName = "DynamicsProcess"
        self.FSWProcessName = "FSWProcess"
        # Create processes
        self.dynProc = self.CreateNewProcess(self.DynamicsProcessName)
        self.fswProc = self.CreateNewProcess(self.FSWProcessName)
        # Crate sim subclasses
        self.DynClass = DINO_DKE.DynamicsClass(self)
        self.FSWClass = DINO_FSW.FSWClass(self)

if __name__ == "__main__":
    scene.attFilter_dynScenario(DINO_DynSim())
