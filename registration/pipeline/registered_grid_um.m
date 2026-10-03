function px_um = registered_grid_um()
%REGISTERED_GRID_UM  Voxel size of every brain's registered volumes, in micrometres.
%   px_um = REGISTERED_GRID_UM() returns 10. run_register_to_atlas 'register'
%   and run_add_sep_channel sample each slice at this spacing when they place
%   it in the atlas (LightSuite's px_atlas, in generateRegisteredSliceVolume),
%   whatever px_atlas the brain's sliceinfo.mat holds, so every brain, adult or
%   young, extracted once or again, comes out on the same grid: the one
%   run_collect_by_group stacks and the Python route block-averages to 20 um.
%
%   LightSuite places the slices in 3D on a grid of half the registration voxel
%   (px_register), so the two spacings agree only with px_register = 20, which
%   run_register_to_atlas checks before registering.
%
%   Run by register_to_atlas and add_sep_channel.

px_um = 10;

end
