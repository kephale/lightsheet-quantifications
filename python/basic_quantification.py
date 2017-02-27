
from net.imglib2 import Point
from net.imglib2.algorithm.region.hypersphere import HyperSphere
from ij import IJ
from mpicbg.imglib.interpolation.nearestneighbor import NearestNeighborInterpolatorFactory;
import random
import time
from net.imglib2.util import Intervals
from ij3d import Image3DUniverse
from ij3d import ImageJ3DViewer
from org.scijava.vecmath import Color3f
import customnode
from ij3d.gui import InteractiveMeshDecimation
from ij3d.gui import InteractiveMeshVoxelization
from ij import WindowManager
import os.path
from process3d import Flood_Fill
from ij.plugin import ImageCalculator
import time
from ij.measure import Measurements
from ij.process import StackStatistics

hist_min = 0
hist_max = 25
hist_step = 1
hist_nBins = int( ( hist_max - hist_min ) / hist_step )
#edt_stats = edt.getStatistics( Measurements.MEDIAN, hist_nBins, hist_min, hist_max )
hist_binLabels = [ ( hist_min + hist_step * el ) for el in range( hist_nBins ) ]

all_edt_histograms = {}
all_skeleton_edt_histograms = {}
for k in hist_binLabels:
    all_edt_histograms[k] = []
    all_skeleton_edt_histograms[k] = []
all_edt_histograms['name'] = [  ]
all_skeleton_edt_histograms['name'] = [  ]

def process_imgdir( directory, close_after ) :
    print( 'Processing: ' + directory )

    original_image_filename = directory + '/uniform.tif'
    probabilities_image_filename = directory + '/probabilities.tif'
    cube_rois_filename = directory + '/CubeROIs.csv'

    cubes = []
    f = open(cube_rois_filename)
    lines = f.readlines()
    for line in lines[1:]:
        #els = [ int(el)-1 for el in line.strip().split(',') ]
        els = [ int(el) for el in line.strip().split(',') ]
        cubes += [ { 'idx': els[0], 'class': els[1], 'x1':min(els[2],els[5]), 'y1':min(els[3],els[6]), 'z1':min(els[4],els[7]), 'x2':max(els[5],els[2]), 'y2':max(els[6],els[3]), 'z2':max(els[4],els[7]) } ]
    f.close()
    print( 'Read ' + str(len(cubes)) + ' cubes from file' )

    #for cube in cubes:
    #for cube_idx in range( len( cubes ) ):
    for cube_idx in [0]:
        cube = cubes[cube_idx]
        print( cube )

        all_edt_histograms['name'] += [ directory.split( '/' )[-1] + '_' + str(cube_idx) ]
        all_skeleton_edt_histograms['name'] += [ directory.split( '/' )[-1] + '_' + str(cube_idx) ]

        cube_basename = directory + '/cube_' + str(cube_idx) + '_'

        # Check if we have a parameters file, and evaluate it if so, otherwise make one with defaults
        param_filename = cube_basename + 'parameters.py'
        # Make a default to start with
        no_params_set = False
        if not os.path.isfile( param_filename ):
            f = open( param_filename, 'w' )
            f.write( 'mesh_thresh=200\n' )
            f.write( 'interior_coords=[ (1,1,1) ]\n' )
            f.close()
            no_params_set = True

        if os.path.isfile( param_filename ):
            f = open( param_filename )
            lines = f.readlines()
            for line in lines:
                if 'mesh_thresh' in line:
                    mesh_thresh = eval(line.split('=')[1])
                if 'interior_coords' in line:
                    interior_coords = eval(line.split('=')[1])

        load_uniform = True
        if load_uniform:
            #open uniform
            uniform = IJ.openImage( original_image_filename )

            #crop uniform
            uniform.setRoi( cube['x1'], cube['y1'], cube['x2'], cube['y2'] )
            IJ.run( uniform, "Crop", "")
            IJ.run( uniform, "Make Substack...", "slices=" + str(cube['z1']) + '-' + str(cube['z2']) )
            uniform.close()
            uniform = IJ.getImage()
            IJ.saveAsTiff( uniform, cube_basename + 'raw.tif' )

            #save max_projection of uniform
            IJ.run(uniform, "Z Project...", "projection=[Max Intensity]")
            uniform_mp = IJ.getImage()
            IJ.run(uniform_mp, "Enhance Contrast", "saturated=0.35")
            IJ.saveAsTiff( uniform_mp, cube_basename + 'raw_maxproj.tif' )

        #open probabilites
        probabilities = IJ.openImage( probabilities_image_filename )

        #crop probabilites
        probabilities.setRoi( cube['x1'], cube['y1'], cube['x2'], cube['y2'] )
        IJ.run( probabilities, "Crop", "")
        IJ.run( probabilities, "Make Substack...", "slices=" + str(cube['z1']) + '-' + str(cube['z2']) )
        time.sleep( 5 )
        probabilities.close()
        probabilities = IJ.getImage()
        #IJ.saveAsTiff( probabilities, cube_basename + 'probabilities.tif' )

        time.sleep( 5 )

        #save max_projection of probabilites
        IJ.run(probabilities, "Z Project...", "projection=[Max Intensity]")
        probabilities_mp = IJ.getImage()
        IJ.run(probabilities_mp, "Enhance Contrast", "saturated=0.35")
        IJ.saveAsTiff( probabilities_mp, cube_basename + 'probabilities_maxproj.tif' )

        time.sleep( 5 )

        # blur probabilities to help smooth the isosurfaces
        blur_radius = 2
        IJ.run( probabilities, "Gaussian Blur 3D...", 'x=' + str(blur_radius) + ' y=' + str(blur_radius) + ' z=' + str(blur_radius) )

        #threshold probabilities
        #IJ.run(probabilities, "Enhance Contrast", "saturated=0.35")
        IJ.run( probabilities, "8-bit", "")
        #IJ.saveAsTiff( probabilities, cube_basename + 'probabilities.tif' )
        #for k in range( probabilities.getNSlices() ):
        #    probabilities.setSlice(k+1)
        #    probabilities.getProcessor().setThreshold( mesh_thresh - 1, mesh_thresh + 1, 0 )
        #IJ.run( probabilities, 'Convert to Mask', 'method=Default background=Dark black')
        #time.sleep( 5 )
        IJ.saveAsTiff( probabilities, cube_basename + 'probabilities.tif' )

        time.sleep( 5 )

        #mesh binary_probabilities
        #IJ.run( "3D Viewer", "")


        threedType = 2
        threedName = 'cube_' + str(cube_idx)

        IJ.runPlugIn( "ij3d.ImageJ3DViewer", probabilities.getTitle() )
        univ = Image3DUniverse.universes.get(0)
        univ.addContent( probabilities, Color3f(1,1,1), threedName, mesh_thresh, [True, False, False], 2, threedType )
        ImageJ3DViewer.select( threedName )
        ImageJ3DViewer.exportContent( 'wavefront', cube_basename + 'mesh.obj' )

        #smooth mesh
        c = univ.getSelected()
        n = c.getContent()
        ctm = n.getMesh()
        fim = customnode.FullInfoMesh( ctm.getMesh() )
        ec = customnode.EdgeContraction( fim, False )

        initial_num_verts = ec.getVertexCount()
        num_to_remove = int( initial_num_verts * 0.1 )

        #v = InteractiveMeshDecimation.simplify( ec, num_to_remove )

        enable_smoothing = False
        if enable_smoothing:
            part = num_to_remove / 10
            last = num_to_remove % 10
            ret = 0
            for i in range(10):
                IJ.showProgress(i + 1, 10)
                ret = ec.removeNext(part)
            if (last != 0):
                ret = ec.removeNext(last)
            IJ.showProgress(1)

            ctm.setMesh( fim.getMesh() )
            ImageJ3DViewer.exportContent( 'wavefront', cube_basename + 'smooth_mesh.obj' )

        # 3d viewer screenshot
        screenshot_3d = univ.takeSnapshot()
        IJ.saveAsTiff( screenshot_3d, cube_basename + 'mesh_screenshot.tif' )

        #voxelize mesh
        voxelizer = InteractiveMeshVoxelization()
        voxelizer.voxelize( ctm, probabilities.getWidth(), probabilities.getHeight(), probabilities.getStackSize() )

        #voxelization = WindowManager.getImage( ctm.getName() + '_voxelization' )
        voxelization = WindowManager.getImage( 'null_voxelization' )
        voxelization.setCalibration( probabilities.getCalibration().copy() )
        #manually 3D fill vasculature, fill with thresholdable-color
        #save selected vasculature

        #call('process3d.Flood_Fill.fill', x,y,z)
        if not no_params_set:
            fill_color = 100
            Flood_Fill.fill( voxelization, interior_coords[0][0], interior_coords[0][1], interior_coords[0][2], fill_color )

            # Threshold to extract
            #IJ.setAutoThreshold( voxelization, 'Default dark' )
            #Prefs.blackBackground = true
            #IJ.run( voxelization, 'Convert to Mask', 'method=Default background=Dark black')

            for k in range( voxelization.getNSlices() ):
                voxelization.setSlice(k+1)
                voxelization.getProcessor().setThreshold( fill_color - 1, fill_color + 1, 0 )
            IJ.run( voxelization, 'Convert to Mask', 'method=Default background=Dark black')
            IJ.saveAsTiff( voxelization, cube_basename + 'voxelization.tif' )

            # Calculate and record volume HERE

            IJ.run(voxelization, "Z Project...", "projection=[Max Intensity]")
            voxelization_mp = IJ.getImage()
            IJ.run(voxelization_mp, "Enhance Contrast", "saturated=0.35")
            IJ.saveAsTiff( voxelization_mp, cube_basename + 'voxelization_maxproj.tif' )

            # Calculate EDT
            IJ.run( voxelization, "Exact Euclidean Distance Transform (3D)", "" )
            edt = WindowManager.getImage( 'EDT' )
            IJ.run( edt, 'Fire', '' )

            # Get the histogram data in an array

            hist_nBins = int( ( hist_max - hist_min ) / hist_step )
            #edt_stats = edt.getStatistics( Measurements.MEDIAN, hist_nBins, hist_min, hist_max )
            edt_stats = StackStatistics( edt, hist_nBins, hist_min, hist_max )
            hist_data = edt_stats.getHistogram()
            hist_binLabels = [ ( hist_min + hist_step * el ) for el in range( hist_nBins ) ]

            max_radius = 20
            IJ.run( edt, 'Histogram', 'bins=' + str(hist_nBins) + ' x_min=' + str(hist_min) + ' x_max=' + str(hist_max) + ' y_max=Auto stack' )
            edt_histogram = WindowManager.getImage( 'Histogram of EDT' )
            IJ.saveAsTiff( edt_histogram, cube_basename + 'edt_histogram.tif' )

            x_unit_scale = uniform.getCalibration().getX(1)
            y_unit_scale = uniform.getCalibration().getY(1)
            z_unit_scale = uniform.getCalibration().getZ(1)

            f = open( cube_basename + 'edt_histogram.csv', 'w' )
            for k in range( len( hist_data ) ):
                f.write( str(hist_binLabels[k]) + '\t' + str(hist_data[k]*x_unit_scale*y_unit_scale*z_unit_scale) + '\n' )
                #f.write( str(hist_binLabels[k]) + '\t' + str(hist_data[k]) + '\n' )
                all_edt_histograms[hist_binLabels[k]] += [ hist_data[k] ]

            f.close()

            # Handling skeletons
            IJ.run( voxelization, "Skeletonize (2D/3D)", "")
            skeleton = voxelization # For simplicity later, but note that voxelization has been mutated
            IJ.run(skeleton, "32-bit", "")
            IJ.run(skeleton, "Calculator Plus", 'i1=' + str(skeleton.getTitle()) + ' i2=' + str(skeleton.getTitle()) + ' operation=[Scale: i2 = i1 x k1 + k2] k1=0.003921568627 k2=0' )

            IJ.run(skeleton, "Z Project...", "projection=[Max Intensity]")
            skeleton_mp = IJ.getImage()
            IJ.run(skeleton_mp, "Enhance Contrast", "saturated=0.35")
            IJ.saveAsTiff( skeleton_mp, cube_basename + 'skeleton_maxproj.tif' )

            ic = ImageCalculator()
            skeleton_edt = ic.run("Multiply 32-bit stack", edt, skeleton)
            IJ.run( skeleton_edt, 'Fire', '' )


            # Get the histogram data in an array
            hist_nBins = int( ( hist_max - hist_min ) / hist_step )
            #skeleton_edt_stats = edt.getStatistics( Measurements.MEDIAN, hist_nBins, hist_min, hist_max )
            skeleton_edt_stats = StackStatistics( edt, hist_nBins, hist_min, hist_max )
            hist_data = skeleton_edt_stats.getHistogram()
            hist_binLabels = [ ( hist_min + hist_step * el ) for el in range( hist_nBins ) ]

            IJ.run( skeleton_edt, 'Histogram', 'bins=' + str(hist_nBins) + ' x_min=' + str(hist_min) + ' x_max=' + str(hist_max) + ' y_max=Auto stack' )
            skeleton_edt_histogram = WindowManager.getImage( 'Histogram of EDT' )
            IJ.saveAsTiff( skeleton_edt_histogram, cube_basename + 'skeleton_edt_histogram.tif' )

            unit_scale = uniform.getCalibration().getX(1)

            f = open( cube_basename + 'skeleton_edt_histogram.csv', 'w' )
            for k in range( len( hist_data ) ):
                f.write( str(hist_binLabels[k]) + '\t' + str(hist_data[k] * unit_scale) + '\n' )
                all_skeleton_edt_histograms[hist_binLabels[k]] += [ hist_data[k] * unit_scale ]
            f.close()

        else:
            IJ.saveAsTiff( voxelization, cube_basename + 'voxelization.tif' )

        # Free up resources
        if close_after:
            ImageJ3DViewer.close()
            IJ.run( 'Close All', '' )
            IJ.freeMemory()

parent_directory = '/Volumes/Storage/LightsheetQuant/'

# dirlist = [ parent_directory + 'Near_optic_nerve_20x', \
#             parent_directory + 'Near_optic_nerve_20xb', \
#             parent_directory + 'Near_optic_nerve_20xc2', \
#             parent_directory + 'Near_optic_nerve_20xd', \
#             parent_directory + 'Near_optic_nerve_20xe', \
#             parent_directory + 'Near_optic_nerve_20xf', \
#             parent_directory + 'Vascular_front_20x', \
#             parent_directory + 'Vascular_front_20xb', \
#             parent_directory + 'Vascular_front_20xc', \
#             parent_directory + 'Vascular_front_20xc2', \
#             parent_directory + 'Vascular_front_20xd', \
#             parent_directory + 'Vascular_front_20xe' ]

# For crops only
dirlist = [ parent_directory + 'Near_optic_nerve_20x', \
            parent_directory + 'Near_optic_nerve_20xc2', \
            parent_directory + 'Near_optic_nerve_20xd', \
            parent_directory + 'Near_optic_nerve_20xe' ]

#process_imgdir( dirlist[1], False )
#process_imgdir( dirlist[10], False )

for directory in dirlist:
    process_imgdir( directory, True )

skeleton_f = open( parent_directory + 'skeleton_edt_histograms.csv', 'w' )
edt_f = open( parent_directory + 'edt_histograms.csv', 'w' )
skeleton_f.write( 'Name' + '\t' + '\t'.join( all_skeleton_edt_histograms['name'] ) + '\n' )
edt_f.write( 'Name' + '\t' + '\t'.join( all_edt_histograms['name'] ) + '\n' )
for k in hist_binLabels:
    skeleton_f.write( str(k) + '\t' + '\t'.join( [ str(el) for el in all_skeleton_edt_histograms[k] ] ) + '\n' )
    edt_f.write( str(k) + '\t' + '\t'.join( [ str(el) for el in all_edt_histograms[k] ] ) + '\n' )
skeleton_f.close()
edt_f.close()
