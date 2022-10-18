#-----------------------------------------------------------------------------
#running different configuraitons of tilted V
#-----------------------------------------------------------------------------
import numpy as np
import os
from parflow import Run
from parflow.tools.fs import mkdir, cp, get_absolute_path
from parflowio.pyParflowio import PFData
from parflow.tools import io
from Make_Solid_Files import make_solid_file


overland = Run("overland_tiltedV", __file__)

#-----------------------------------------------------------------------------

overland.FileVersion = 4

overland.Process.Topology.P = 1
overland.Process.Topology.Q = 1
overland.Process.Topology.R = 1

#---------------------------------------------------------
# Computational Grid
#---------------------------------------------------------
gridnx = int(25)
gridny = int(25)
gridnz = int(2)
griddx = 10.0
griddy = 10.0
griddz = 5.0

overland.ComputationalGrid.Lower.X = 0.0
overland.ComputationalGrid.Lower.Y = 0.0
overland.ComputationalGrid.Lower.Z = 0.0

overland.ComputationalGrid.NX = gridnx
overland.ComputationalGrid.NY = gridny
overland.ComputationalGrid.NZ = gridnz

overland.ComputationalGrid.DX = griddx
overland.ComputationalGrid.DY = griddy
overland.ComputationalGrid.DZ = griddz

#---------------------------------------------------------
# The Names of the GeomInputs
#---------------------------------------------------------
overland.GeomInput.Names = "domaininput"

overland.GeomInput.domaininput.GeomName = 'domain'
overland.GeomInput.domaininput.GeomNames = 'domain'
overland.GeomInput.domaininput.InputType = 'SolidFile'
#overland.GeomInput.domaininput.FileName = 'TV25_Yriver.pfsol'
overland.Geom.domain.Patches = "reservoir bottom slope1 channel slope2 side"

#-----------------------------------------------------------------------------
# Perm
#-----------------------------------------------------------------------------
overland.Geom.Perm.Names = 'domain'
overland.Domain.GeomName = "domain"
overland.Geom.domain.Perm.Type = 'Constant'
overland.Geom.domain.Perm.Value = 0.0000001

overland.Perm.TensorType = 'TensorByGeom'

overland.Geom.Perm.TensorByGeom.Names = 'domain'

overland.Geom.domain.Perm.TensorValX = 1.0
overland.Geom.domain.Perm.TensorValY = 1.0
overland.Geom.domain.Perm.TensorValZ = 1.0

#-----------------------------------------------------------------------------
# Specific Storage
#-----------------------------------------------------------------------------
overland.SpecificStorage.Type = 'Constant'
overland.SpecificStorage.GeomNames = 'domain'
overland.Geom.domain.SpecificStorage.Value = 1.0e-4

#-----------------------------------------------------------------------------
# Phases
#-----------------------------------------------------------------------------
overland.Phase.Names = 'water'

overland.Phase.water.Density.Type = 'Constant'
overland.Phase.water.Density.Value = 1.0

overland.Phase.water.Viscosity.Type = 'Constant'
overland.Phase.water.Viscosity.Value = 1.0

#-----------------------------------------------------------------------------
# Contaminants
#-----------------------------------------------------------------------------
overland.Contaminants.Names = ''

#-----------------------------------------------------------------------------
# Retardation
#-----------------------------------------------------------------------------
overland.Geom.Retardation.GeomNames = ''

#-----------------------------------------------------------------------------
# Gravity
#-----------------------------------------------------------------------------
overland.Gravity = 1.0

#-----------------------------------------------------------------------------
# Setup timing info
#-----------------------------------------------------------------------------
overland.TimingInfo.BaseUnit = 1.0
overland.TimingInfo.StartCount = 0
overland.TimingInfo.StartTime = 0.0
overland.TimingInfo.StopTime = 10.0
# overland.TimingInfo.StopTime = 5.0
overland.TimingInfo.DumpInterval = 1.0
overland.TimeStep.Type = 'Constant'
overland.TimeStep.Value = 0.001

#-----------------------------------------------------------------------------
# Porosity
#-----------------------------------------------------------------------------

overland.Geom.Porosity.GeomNames = 'domain'
overland.Geom.domain.Porosity.Type = 'Constant'
overland.Geom.domain.Porosity.Value = 0.01

#-----------------------------------------------------------------------------
# Domain
#-----------------------------------------------------------------------------
overland.Domain.GeomName = 'domain'

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------
overland.Phase.RelPerm.Type = 'VanGenuchten'
overland.Phase.RelPerm.GeomNames = 'domain'

overland.Geom.domain.RelPerm.Alpha = 6.0
overland.Geom.domain.RelPerm.N = 2.

#---------------------------------------------------------
# Saturation
#---------------------------------------------------------
overland.Phase.Saturation.Type = 'VanGenuchten'
overland.Phase.Saturation.GeomNames = 'domain'

overland.Geom.domain.Saturation.Alpha = 6.0
overland.Geom.domain.Saturation.N = 2.
overland.Geom.domain.Saturation.SRes = 0.2
overland.Geom.domain.Saturation.SSat = 1.0

#-----------------------------------------------------------------------------
# Wells
#-----------------------------------------------------------------------------
overland.Wells.Names = ''

#-----------------------------------------------------------------------------
# Time Cycles
#-----------------------------------------------------------------------------
overland.Cycle.Names = 'constant rainrec'
overland.Cycle.constant.Names = 'alltime'
overland.Cycle.constant.alltime.Length = 1
overland.Cycle.constant.Repeat = -1

overland.Cycle.rainrec.Names = 'rain rec'
overland.Cycle.rainrec.rain.Length = 10
overland.Cycle.rainrec.rec.Length = 300
overland.Cycle.rainrec.Repeat = -1

#-----------------------------------------------------------------------------
# Boundary Conditions: Pressure
#-----------------------------------------------------------------------------
# changed this
overland.BCPressure.PatchNames = overland.Geom.domain.Patches

overland.Patch.side.BCPressure.Type = 'FluxConst'
overland.Patch.side.BCPressure.Cycle = 'constant'
overland.Patch.side.BCPressure.alltime.Value = 0.0

overland.Patch.bottom.BCPressure.Type = 'FluxConst'
overland.Patch.bottom.BCPressure.Cycle = 'constant'
overland.Patch.bottom.BCPressure.alltime.Value = 0.0

overland.Patch.slope1.BCPressure.Type = 'OverlandKinematic'
overland.Patch.slope1.BCPressure.Cycle = 'constant'
overland.Patch.slope1.BCPressure.alltime.Value = -0.01

overland.Patch.slope2.BCPressure.Type = 'OverlandKinematic'
overland.Patch.slope2.BCPressure.Cycle = 'constant'
overland.Patch.slope2.BCPressure.alltime.Value = -0.01

overland.Patch.channel.BCPressure.Type = 'OverlandKinematic'
overland.Patch.channel.BCPressure.Cycle = 'constant'
overland.Patch.channel.BCPressure.alltime.Value = -0.01

overland.Patch.reservoir.BCPressure.Type = 'OverlandKinematic'
overland.Patch.reservoir.BCPressure.Cycle = 'constant'
overland.Patch.reservoir.BCPressure.alltime.Value = -0.01

# overland.Patch.reservoir.BCPressure.Type = 'SeepageFace'
# overland.Solver.OverlandKinematic.SeepageOne = 1
# overland.Solver.OverlandKinematic.SeepageTwo = 2

#---------------------------------------------------------
# Mannings coefficient
#---------------------------------------------------------
overland.Mannings.Type = 'Constant'
overland.Mannings.GeomNames = 'domain'
overland.Mannings.Geom.domain.Value = 3.e-6

#-----------------------------------------------------------------------------
# Phase sources:
#-----------------------------------------------------------------------------
overland.PhaseSources.water.Type = 'Constant'
overland.PhaseSources.water.GeomNames = 'domain'
overland.PhaseSources.water.Geom.domain.Value = 0.0

#-----------------------------------------------------------------------------
# Exact solution specification for error calculations
#-----------------------------------------------------------------------------
overland.KnownSolution = 'NoKnownSolution'

#-----------------------------------------------------------------------------
# Slopes
#-----------------------------------------------------------------------------
overland.TopoSlopesX.Type = "PFBFile"
overland.TopoSlopesX.GeomNames = "domain"
overland.TopoSlopesX.FileName = "slopex.pfb"

overland.TopoSlopesY.Type = "PFBFile"
overland.TopoSlopesY.GeomNames = "domain"
overland.TopoSlopesY.FileName = "slopey.pfb"

#-----------------------------------------------------------------------------
# Set solver parameters
#-----------------------------------------------------------------------------
overland.Solver = 'Richards'
overland.Solver.MaxIter = 50000

overland.Solver.Nonlinear.MaxIter = 100
overland.Solver.Nonlinear.ResidualTol = 1e-9
overland.Solver.Nonlinear.EtaChoice = 'EtaConstant'
overland.Solver.Nonlinear.EtaValue = 0.01
# Changed this
# overland.Solver.Nonlinear.UseJacobian = False
overland.Solver.Nonlinear.DerivativeEpsilon = 1e-15
overland.Solver.Nonlinear.StepTol = 1e-20
overland.Solver.Nonlinear.Globalization = 'LineSearch'
overland.Solver.Linear.KrylovDimension = 50
overland.Solver.Linear.MaxRestart = 2
overland.Solver.OverlandKinematic.Epsilon = 1E-5
# Changed this
overland.Solver.TerrainFollowingGrid = True

overland.Solver.Linear.Preconditioner = 'PFMG'
overland.Solver.PrintSubsurf = True
overland.Solver.PrintSlopes = True
overland.Solver.PrintMannings = True
overland.Solver.Drop = 1E-20
overland.Solver.AbsTol = 1E-10


## run with KWE upwinding and analytical jacobian
# Changed this
overland.Solver.Nonlinear.UseJacobian = True
overland.Solver.Linear.Preconditioner.PCMatrixType = 'PFSymmetric'
#
## run with KWE upwinding and analytical jacobian and nonsymmetric preconditioner
# overland.Solver.Nonlinear.UseJacobian = True
# overland.Solver.Linear.Preconditioner.PCMatrixType = 'FullJacobian'

#overland.Solver.Nonlinear.UseJacobian = False
#overland.Solver.Linear.Preconditioner.PCMatrixType = 'PFSymmetric'

#---------------------------------------------------------
# Initial conditions: water pressure
#---------------------------------------------------------
# set water table to be at the bottom of the domain, the top layer is initially dry
INITIAL_PRESSURE_FILE = "../initial_pressure.pfb"
overland.ICPressure.Type =                                  "PFBFile"
overland.Geom.domain.ICPressure.FileName	=		              INITIAL_PRESSURE_FILE

overland.ICPressure.GeomNames   =                           "domain"
#overland.Geom.domain.ICPressure.Value =                     1191

overland.Geom.domain.ICPressure.RefGeom  =                  "domain"
overland.Geom.domain.ICPressure.RefPatch  =                 "bottom"


#-----------------------------------------------------------------------------
# defining write function to save some space in the loop
#-----------------------------------------------------------------------------
def write_pfb_to_run_dir(myarray, fout, run_dir):
    array_pfb = PFData(myarray)
    fout = get_absolute_path(os.path.join(run_dir, fout))
    array_pfb.writeFile(fout)
    array_pfb.close()

def dist_and_run(run_dir):
    overland.ComputationalGrid.NZ = 1
    overland.dist(os.path.join(run_dir, 'slopex.pfb'))
    overland.dist(os.path.join(run_dir, 'slopey.pfb'))
    overland.ComputationalGrid.NZ = gridnz

    overland.write(working_directory=run_dir)
    overland.write(file_format='yaml', working_directory=run_dir)
    overland.write(file_format='json', working_directory=run_dir)
    overland.run(working_directory=run_dir)

#-----------------------------------------------------------------------------
# setup solid files in x and y direction for the grid dimensions
#-----------------------------------------------------------------------------
# Make a y direction solid file
solid_name = "Yriver"+"_nx" + str(gridnx) + '_dz' + str(griddz)
solid_fname = solid_name + '.pfsol'
make_solid_file(nx=gridnx, ny=gridny,
                bottom_val=2, side_val=6, top_val1=3, top_val2=4, top_val3=5, reservoir_value=1,
                latsize=griddx, zdepth=griddz*gridnz, river_dir=2,
                root_name=solid_name, out_dir='./', pftools_path="/usr/local/bin")

# Make an x direction solid file
solid_name = "Xriver"+"_nx" + str(gridnx) + '_dz' + str(griddz)
solid_fname = solid_name + '.pfsol'
make_solid_file(nx=gridnx, ny=gridny,
                bottom_val=2, side_val=6, top_val1=3, top_val2=4, top_val3=5, reservoir_value=1,
                latsize=griddx, zdepth=griddz*gridnz, river_dir=1,
                root_name=solid_name, out_dir='./', pftools_path="/usr/local/bin/")


#-----------------------------------------------------------------------------
# setup the overland patch types 
#-----------------------------------------------------------------------------


#name the patch configuration for run names
patch_name = 'Overland_Kin' 
#-----------------------------------------------------------------------------
# Draining bottom -- (positive Y slope channel)
#-----------------------------------------------------------------------------
# run_dir = patch_name + '_ypos'
# mkdir(run_dir)

# #copy in the solid file
# solid_fname = "Yriver"+"_nx" + str(gridnx) + '_dz' + str(griddz) + '.pfsol'
# cp(solid_fname, run_dir)
# overland.GeomInput.domaininput.FileName = solid_fname

# #Make slope files
# slopex = np.full((gridnz, gridny, gridnx), 1.0)
# slopex[:, :, 0:int(np.floor(gridnx/2))] = np.round(-0.01,2)
# slopex[:, :, int(np.floor(gridnx/2)):] = np.round(0.01,2)
# #slopex[:, :, int(np.floor(gridnx/2))] = 0.01
# #slopex[:, :, (int(np.floor(gridnx/2))+1):] = 0.01
# write_pfb_to_run_dir(slopex, 'slopex.pfb', run_dir)
    
# slopey = np.ones((gridnz, gridny, gridnx))
# slopey = np.round(slopey*0.01,2)
# write_pfb_to_run_dir(slopey, 'slopey.pfb', run_dir)

# dist_and_run(run_dir)

# #-----------------------------------------------------------------------------
# # Draining top -- (negative Y slope channel)
# #-----------------------------------------------------------------------------
# run_dir = patch_name + '_yneg'
# mkdir(run_dir)

# #copy in the solid file
# solid_fname = "Yriver"+"_nx" + str(gridnx) + '_dz' + str(griddz) + '.pfsol'
# cp(solid_fname, run_dir)
# overland.GeomInput.domaininput.FileName = solid_fname

# #Make slope files
# slopex = np.full((gridnz, gridny, gridnx), 1.0)
# slopex[:, :, 0:int(np.floor(gridnx/2))] = np.round(-0.01,2)
# slopex[:, :, int(np.floor(gridnx/2)):] = np.round(0.01, 2)
# #slopex[:, :, int(np.floor(gridnx/2))] = np.round(0.01, 2)
# #slopex[:, :, (int(np.floor(gridnx/2))+1):] = 0.01
# write_pfb_to_run_dir(slopex, 'slopex.pfb', run_dir)

# slopey = np.ones((gridnz, gridny, gridnx))
# slopey = np.round(slopey*-0.01,2)
# write_pfb_to_run_dir(slopey, 'slopey.pfb', run_dir)

# dist_and_run(run_dir)

# #-----------------------------------------------------------------------------
# # Draining left -- (postive X slope channel)
# #-----------------------------------------------------------------------------
# run_dir = patch_name + '_xpos'
# mkdir(run_dir)

# #copy in the solid file
# solid_fname = "Xriver"+"_nx" + str(gridnx) + '_dz' + str(griddz) + '.pfsol'
# cp(solid_fname, run_dir)
# overland.GeomInput.domaininput.FileName = solid_fname

# #Make slope files
# slopey = np.full((gridnz, gridny, gridnx), 1.0)
# slopey[:, 0:int(np.floor(gridny/2)), :] = np.round(-0.01,2)
# slopey[:, int(np.floor(gridny/2)):, :] = np.round(0.01,2)
# #slopey[:, int(np.floor(gridny/2)), :] = 0.01
# #slopey[:, (int(np.floor(gridny/2))+1):, :] = 0.01
# write_pfb_to_run_dir(slopey, 'slopey.pfb', run_dir)

# slopex = np.ones((gridnz, gridny, gridnx))
# slopex = np.round(slopex*0.01,2)
# write_pfb_to_run_dir(slopex, 'slopex.pfb', run_dir)

# dist_and_run(run_dir)
#-----------------------------------------------------------------------------
# Draining right-- (negative X slope channel)
#-----------------------------------------------------------------------------
run_dir = patch_name + '_xneg'
mkdir(run_dir)

#copy in the solid file
solid_fname = "Xriver"+"_nx" + str(gridnx) + '_dz' + str(griddz) + '.pfsol'
cp(solid_fname, run_dir)
overland.GeomInput.domaininput.FileName = solid_fname

#Make slope files
slopey = np.full((gridnz, gridny, gridnx), 1.0)
slopey[:, 0:int(np.floor(gridny/2)), :] = np.round(-0.01,2)
slopey[:, int(np.floor(gridny/2)):, :] = np.round(0.01,2)
slopey[:, int(np.floor(gridny/2)):, :] = np.round(0.01,2)
write_pfb_to_run_dir(slopey, 'slopey.pfb', run_dir)

slopex = np.ones((gridnz, gridny, gridnx))
slopex = np.round(slopex*-0.01,2)
write_pfb_to_run_dir(slopex, 'slopex.pfb', run_dir)

dist_and_run(run_dir)
