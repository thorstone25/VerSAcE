# VSX-OOD
VSX-OOD is an object oriented design wrapper for Verasonics Vantage programming.


## Getting started

You can utilize this project at 3 levels:
* High-level - export a [QUPS](https://github.com/thorstone25/qups) `UltrasoundSystem` configuration to VSX structs
* Mid-level - use the `VSXStruct` wrappers and linking stage to program with a tree structure and avoid indexing errors 
* Low-level - use the `VSXStruct` wrappers for argument validation to avoid typos, then use `struct` to convert properties.


## Setting up a MATLAB project for path management
It may be helpful to create a MATLAB project to manage your QUPS, VSX-OOD, and Vantage paths.
1. Change directory to your vsx-ood repository
2. Open a matlab editor, navigate to the "Home" tab and select New->Project->From Folder->Create (were vsx-ood path is chosen)
3. In the "Project" tab, select References->Browse->*navigate to qups repo*->Add
4. In the "Project" tab, select References->Browse->*navigate to Vantage repo (e.g., Vantage-4.8.4)*->Add
5. In your "Current Folder" MATLAB pane, double-click the "Vsxood.prj" file. You will need to open the file twice due to some MATLAB bug. Now, you should have all necessary QUPS, VSX-OOD, and Vantage paths added to your MATLAB $PATH environment variable.
6. Change directory to your Vantage directory.
7. In the "Home" tab, select Open->*navigate to vsx-ood/scripts dir*->*Choose custom script*. You should be able to run this custom script while your current working dir is still on Vantage.

