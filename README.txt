Author: Luke Bartholomew
Date: September 8 2022

Development:						Completed:
- Single Phase 1D shock tube case
	- Initialisation of problem:
	- flexibility:
- Single Phase 1D combustion chamber	
- Multi-component fluid domain capability
- Boundary conditions implementations:
	OutFlow
	InFlow
	FromStag
	Wall
	FixedP
	Cyclic (maybe)
	






Usage notes:
- Boundary conditions are defined by [int(interfaceID), [str(label), int(ghost cell layers), [optional args]]]
	optional args are dependent of str(label) -> "OutFlow" requires none because it's all zero-gradient
						  -> "InFlow" requires the definition of a flowState
						  -> "Wall" requires none

- momentum flux does not include pressure term, so this must be accounted for separately

- Joining components alters the individual components, so need to take care when adding ghost cells. 
	Can get around this by using deepcopy when inputting the arguments for jointBlock

- Reconstruction

- Python path has been added to in SinglePhaseStraightPipeCell and SinglePhaseInterface. 
	When making new versions of the design tool, need to make sure these point to the new folders, not old ones.

POTENTIAL CHANGES:
- jointBlock module:
	- Join for loops together that are currently separate
	- Update mesh2 cell and interface IDs before forming 
	  joined cellArray and interfaceArray to avoid unnecessary 
	  updates of cell and interface IDs that don't change 
	  (the ones in mesh1).
- Error catching:
	- Checking if boundary interfaces are actually boundary interfaces
	- Checking if east boundaries are attached to west boundaries and vice-versa
		(Don't want to find that correct interface on one block is not joined to
		 a real boundary on the other block)\
	- Check that when joining blocks (except when adding boundary condition), 
		the chosen interface is not linked to a boundary condition.
- Because adding components makes a non-trivial change to the boundary interface indices, 
	would be good to print the new boundary indices after adding components together.
- To avoid reforming interpolation stencils every time step, 
	on the first time step assign a list of indices to where to find the cells.
- For interfaces, split reconstruction stencil into left and right stencils 
	-> Eg [L2, L1, L0] 
	Allows for strictly one sided reconstructions
	


NEED TO DO:
- Multi-component plotting
	- Either write select multiple components to a single output file or read in multiple component files.
	- Writing multiple components to one file has the risk of unmatched variables.
- Write animation modules for single and multi-component animations
	- Be able to model comparable properties between components that might not have exactly the same property names
		- Eg gas pipe "rho" and two phase pipe "rho_g"
	- For simplicity, make all files .gif
	- ALL GIF FILES HAVE TO BE AT 30FPS WHO THE FUCK KNOWS WHY 
- Switch initialisation and boundary condition generation to use Eilmer GasStates rather than specifying properties in lists.

Last Updated:
December 1 2022

