# @Dataset input
# @OpService opService
# @UIService uiService
# @OUTPUT Dataset volume_edt

# We assume the image is sufficiently preprocessed.

print( "Starting script" )

# Minor preprocessing
#blur = opService.filter.mean(input, [2, 2, 2])
#thresholded = opService.threshold.moments(blur)

thresholded = opService.convert().bit(input)

# Compute our 3D mesh
mesh = opService.geom().marchingCubes(thresholded)

print( "Mesh" )
print( mesh )

# Voxelize the mesh
voxelization = opService.geom().voxelization(mesh)

# Some coordinate that is interior to our mesh. This is hard coded, and generally needs to be manually determined per image.
interior_coord = [28, 26, 24]
# TODO: voxelization returns the surface of a mesh. Previous code used `process3d.Flood_Fill` to fill the interior.
#   there is no readily available IJ2 option for flood fill.

voxelized_volume = opService.stats().sum(voxelization)# TODO: once voxelization is filled make sure this is called on correct variable

print( 'Volume is: ' + str(voxelized_volume) )

# Distance transform
volume_edt = opService.image().distancetransform(opService.convert().float32(voxelization))
# TODO: double-check that this returns the EDT for interior of vasculature by saving and inspecting
# uiService.show(volume_edt)

sum_volume_edt = opService.stats().sum(volume_edt)

print( 'Distance to surface is: ' + str(sum_volume_edt) )

# Skeletonizing

skeleton = opService.morphology().thinMorphological(opService.convert().bit(input))

skeleton_edt = opService.math().multiply(skeleton, volume_edt)

sum_skeleton_edt = opService.stats().sum(skeleton_edt)

print( 'Skeleton edt is: ' + str(sum_skeleton_edt) )

