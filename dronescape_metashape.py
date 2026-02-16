#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to initialize a Metashape project with RGB and multispectral images in separate chunks.

Author: Juan C. Montes Herrera (University of Tasmania)

Assumes TERN directory structure:
    <plot>/YYYYMMDD/imagery/
        ├── rgb/level0_raw/
        └── multispec/level0_raw/
User provides:
    --imagery_dir: path to YYYYMMDD/imagery/
    --crs: EPSG code for target CRS (optional, defaults to 4326)
    --out: Metashape project location, not orthomosaics
    --enable_oblique: flag to enable oblique cameras (default: disabled)
Project will be named as "YYYYMMDD-plot.psx"

Example usage:
    # Basic usage with required arguments
    -imagery_dir /path/to/SITE-01/20230615/imagery/ -out /path/to/output/

    # Enable oblique cameras
    -imagery_dir /path/to/SITE-01/20230615/imagery/ -out /path/to/output/ -enable_oblique
"""

import argparse
import os
import sys
import json
import subprocess
import pandas as pd
from pathlib import Path

import Metashape

from functions.gpu_setup import setup_gpu
from functions.utils import find_filtered_images
from functions.camera_ops import id_multispectral_camera
from functions.camera_ops import enable_oblique_cameras
from functions.camera_ops import filter_multispec

####### SfM Processing Parameters #############
## Alignment
downscale_align = 2 # Use 4 for tests
keypoint_limit = 50000 # Default 40K
tiepoint_limit = 0 #Default 4K

## Depth maps
downscale_depthmaps = 2 # 2 High

## Model
model_surface = Metashape.Arbitrary
model_source_data= Metashape.PointCloudData
model_face_count=Metashape.MediumFaceCount

## Reflectance
sun_sensor = False # Only recommended for cloudy conditions.

################################################

# RGB camera offset
p1_cam_offset = (0.087, 0 , 0)

def main():
    # Clear any existing document to ensure clean state for batch processing
    doc = Metashape.app.document
    doc.clear()
    
    # Set up GPU acceleration
    setup_gpu()
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Initialize Metashape project with RGB and multispectral images.")
    parser.add_argument('-imagery_dir', required=True, help='Path to YYYYMMDD/imagery/ directory')
    parser.add_argument('-crs', default="4326", help='EPSG code default to 4326(WGS84), need tests for other codes')
    parser.add_argument('-out', required=True, help='Directory to save the Metashape project')
    parser.add_argument('-enable_oblique', action='store_true', help='Enable oblique cameras (default: disabled)')
    args = parser.parse_args()

    # Extract YYYYMMDD and plot from input path
    imagery_dir = Path(args.imagery_dir).resolve()
    if imagery_dir.name != "imagery":
        print("The -imagery_dir must point to the 'imagery' directory (e.g., <plot>/YYYYMMDD/imagery/)")
    yyyymmdd = imagery_dir.parent.name
    plot = imagery_dir.parent.parent.name
    project_name = f"{yyyymmdd}-{plot}.psx"

    # Set up paths for RGB and multispectral imagery
    rgb_dir = imagery_dir / "rgb" / "level0_raw"
    multispec_dir = imagery_dir / "multispec" / "level0_raw"

    if not rgb_dir.is_dir():
        print(f"RGB directory not found: {rgb_dir}")
    if not multispec_dir.is_dir():
        print(f"Multispec directory not found: {multispec_dir}")

    # Find all RGB images (jpg files)
    rgb_images = find_filtered_images(rgb_dir, extensions=('.jpg', '.jpeg'))
    
    # Find all multispectral images (tif files), excluding Panchro images (ending with _6.tif)
    multispec_images = find_filtered_images(multispec_dir, extensions=('.tif', '.tiff'), exclude_patterns=('_6.tif',))

    if not rgb_images:
        print(f"No RGB images found in {rgb_dir}")
    if not multispec_images:
        print(f"No multispectral images found in {multispec_dir}")

    print(f"Found {len(rgb_images)} RGB images")
    print(f"Found {len(multispec_images)} multispectral images (excluding Panchro band)")
    
    # Initialize Metashape project and create output directory
    doc = Metashape.app.document
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    project_path = out_dir / project_name
    doc.save(str(project_path))

    # Remove default empty chunk if it exists
    if len(doc.chunks) == 1 and doc.chunks[0].label == "Chunk 1" and len(doc.chunks[0].cameras) == 0:
        doc.remove(doc.chunks[0])

    # Create separate chunks for RGB and multispectral images
    rgb_chunk = doc.addChunk()
    rgb_chunk.label = "rgb_images"
    
    multispec_chunk = doc.addChunk()
    multispec_chunk.label = "multispec_images"

    # Add RGB images to rgb chunk
    print(f"Adding {len(rgb_images)} RGB images to the project...")
    rgb_chunk.addPhotos(rgb_images, load_reference=True, load_xmp_calibration=True, 
                        load_xmp_orientation=True, load_xmp_accuracy=True, load_xmp_antenna=True)
    
    # Detect and handle oblique cameras in RGB chunk
    print("Processing RGB cameras...")
    oblique_count, oblique_cameras = enable_oblique_cameras(
        rgb_chunk,
        args.enable_oblique  # If enable_oblique is True, we enable oblique cameras
    )

    doc.save()

    # Add multispectral images to multispectral chunk
    print(f"Adding {len(multispec_images)} multispectral images to the project...")
    multispec_chunk.addPhotos(multispec_images, layout=Metashape.MultiplaneLayout, load_reference=True, 
                              load_xmp_calibration=True, load_xmp_orientation=True, load_xmp_accuracy=True, 
                              load_xmp_antenna=True)
    
    # Disable multispectral imagery reference locations
    for camera in multispec_chunk.cameras:
        camera.reference.location_enabled = False

    if len(rgb_chunk.cameras) == 0:
        print("RGB chunk is empty after adding images.")
    if len(multispec_chunk.cameras) == 0:
        print("Multispectral chunk is empty after adding images.")

    # Detect multispectral camera band label and indices
    id_multispectral_camera(multispec_chunk)
    
    # # Locate reflectance panels
    multispec_chunk.locateReflectancePanels()
    print("Reflectance panel detection complete.")

    # Filter multispectral images based on RGB capture time window using the first oblique camera as end time
    # NOTE: A 17.78 leap second offset is applied to the multispectral images to align with the RGB images.
    print("Filtering multispectral images outside RGB capture time window...")
    disabled_count = filter_multispec(rgb_dir, multispec_chunk, rgb_chunk, oblique_cameras)

    print(f"Disabled {disabled_count} multispectral images outside RGB capture window")

    doc.save()

    # Merge RGB and multispectral chunks
    print("Merging RGB and multispectral chunks...")
    # mergeChunks creates a new chunk, store the result directly
    merged_chunk = doc.mergeChunks([rgb_chunk, multispec_chunk])
    doc.chunks[2].label = "merged_chunk"
    
    doc.save()
    # Remove the original chunks
    print("Removing original RGB and multispectral chunks...")
    
    doc.remove([rgb_chunk, multispec_chunk])
    doc.save()
    
    print('Completed project setup and camera synchronization.')

    merged_chunk = doc.chunk

    # Apply camera offset to the first camera
    merged_chunk.sensors[0].antenna.location_ref = Metashape.Vector(p1_cam_offset)
    print(f"Camera offset applied to the first camera: {merged_chunk.sensors[0].antenna.location_ref}")

    # Set CRS
    crs_code = args.crs
    target_crs = Metashape.CoordinateSystem(f"EPSG::{crs_code}")
    merged_chunk.crs = target_crs

    print("Aligning images...")
    # Match photos with specified settings
    merged_chunk.matchPhotos(
        downscale=downscale_align,
        generic_preselection=True,  # Enable generic preselection
        reference_preselection=True,  # Enable reference preselection
        reference_preselection_mode=Metashape.ReferencePreselectionSource,  # Source mode
        keypoint_limit=keypoint_limit,  # Key points limit
        tiepoint_limit=tiepoint_limit,  # Tie points limit
        filter_stationary_points=True,  # Exclude stationary points
        keep_keypoints=False,
        guided_matching=False,
    )
    
    # Align cameras
    merged_chunk.alignCameras(
        min_image=2, #Min image count per tie point - Default 2
        adaptive_fitting=True, # Default? Refines distortion locally
        subdivide_task=True,
    )

    print("Image alignment complete!")

    # Filter tie points by reproj. error and uncertainty
    f = Metashape.TiePoints.Filter()
    f.init(merged_chunk, Metashape.TiePoints.Filter.ReprojectionError)
    # How accurately a point can be located in the image
    f.selectPoints(0.7) # Agisoft best-practice 0.5 / Permissive 1.0
    f.removePoints(threshold=0.7)

    f.init(merged_chunk, Metashape.TiePoints.Filter.ReconstructionUncertainty)
    # Estimated 3-D positional uncertainty for each tie point
    f.selectPoints(15) # Agisoft best-practice 10 / Permissive 20-30
    f.removePoints(threshold=15)

    print("Best-practice gradual filter complete!")
    
    # Optimize cameras with specified parameters
    merged_chunk.optimizeCameras(
        fit_f=True,  # Fit focal length
        fit_cx=True, fit_cy=True, # Fit principal point x, y
        fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=False, # Fit radial distortion k1, k2, k3 
        fit_p1=True, fit_p2=True, # Fit tangential distortion p1, p2
        fit_b1=True, fit_b2=True, # Fit affinity b1, b2
        fit_corrections=True,  # Fit additional corrections
        )
    
    print("Optimization complete!")

    print("Building depth maps...")
    merged_chunk.buildDepthMaps(
        downscale=downscale_depthmaps,  # Quality (4=Low, 2=Medium, 1=High, 0=Ultra)
        filter_mode=Metashape.NoFiltering,  # Options: NoFiltering, MildFiltering, ModerateFiltering, AggressiveFiltering
        reuse_depth=False,
        max_neighbors=20, # Default 16
        subdivide_task=True,
        workitem_size_cameras=30, #Default 20
        max_workgroup_size=150 #Default 100
        )
    
    merged_chunk.buildPointCloud(
        source_data=Metashape.DepthMapsData,
        point_colors=False,
        point_confidence=False,
        keep_depth=False,
        max_neighbors=20, # Default 16
        uniform_sampling=False,
        # points_spacing=0.005,
        subdivide_task=True,
        workitem_size_cameras=30, # Default 20
        max_workgroup_size=150 # Default 100
        )

    # Build the model
    merged_chunk.buildModel(
        surface_type=model_surface,
        source_data=model_source_data, 
        face_count=model_face_count,
        # face_count=Metashape.CustomFaceCount,
        # face_count_custom=200000,
        interpolation=Metashape.DisabledInterpolation, # Default Enabled
        trimming_radius=5, #Default 10
        build_texture=False,
        vertex_colors=False,
        vertex_confidence=False,
        keep_depth=False,
        subdivide_task=True,
        workitem_size_cameras=30, # Default 20
        max_workgroup_size=150 # Default 100
    )
    
    print("Model building complete!")

    merged_duplicate = merged_chunk.copy()
    merged_duplicate.label = "merged_duplicate"
    doc.chunks.append(merged_duplicate)
    doc.save()

    # Remove multispec (.TIF) cameras from RGB chunk (merged_chunk)
    tif_cams = [cam for cam in merged_chunk.cameras if cam.photo.path.lower().endswith(".tif")]
    if tif_cams:
        merged_chunk.remove(tif_cams)

    doc.save()

    # Remove RGB (.JPG) cameras from Multispec chunk (merged_duplicate)
    jpg_cams = [cam for cam in merged_duplicate.cameras if cam.photo.path.lower().endswith(".jpg")]
    if jpg_cams:
        merged_duplicate.remove(jpg_cams)

    doc.save()

    print("Project setup complete!")
    print("###########################")

     # Calibrate reflectance of Multispec chunk 
    merged_duplicate.calibrateReflectance(use_reflectance_panels=True, use_sun_sensor=sun_sensor)

    # Raster transform multispectral images
    print("Updating Raster Transform for relative reflectance")
    raster_transform_formula = []
    num_bands = len(merged_duplicate.sensors)
    print(f"Number of bands: {num_bands}")
    for band in range(1, num_bands):
        raster_transform_formula.append("B" + str(band) + "/32768")

    merged_duplicate.raster_transform.formula = raster_transform_formula
    merged_duplicate.raster_transform.calibrateRange()
    merged_duplicate.raster_transform.enabled = True
    doc.save()
    print(f"Applied raster transform formulas: {raster_transform_formula}")
    
    # Build RGB orthomosaic
    merged_chunk.buildOrthomosaic(
        surface_data=Metashape.DataSource.ModelData,
        blending_mode=Metashape.MosaicBlending,
        refine_seamlines=True,
        fill_holes=True
    )
    doc.save() 
    
    # Build Multispectral orthomosaic for duplicate chunk
    merged_duplicate.buildOrthomosaic(
        surface_data=Metashape.DataSource.ModelData,
        # blending_mode=Metashape.MosaicBlending,
        blending_mode=Metashape.AverageBlending,
        refine_seamlines=True,
        fill_holes=True,
        # resolution=0.05
    )
    doc.save()
    print("RGB ortho:", merged_chunk.orthomosaic)
    print("MS ortho:", merged_duplicate.orthomosaic)

    rgb_res = round(merged_chunk.orthomosaic.resolution, 2)
    ms_res = round(merged_duplicate.orthomosaic.resolution, 2)
    print("RGB ortho resolution:", rgb_res)
    print("MS ortho resolution:", ms_res)


    # Set up output paths for RGB and multispectral orthomosaics using string paths
    rgb_out_dir = os.path.join(str(imagery_dir), "rgb", "level1_proc")
    multispec_out_dir = os.path.join(str(imagery_dir), "multispec", "level1_proc")
    
    # Create output directories if they don't exist
    os.makedirs(rgb_out_dir, exist_ok=True)
    os.makedirs(multispec_out_dir, exist_ok=True)
    
    print("RGB out directory:", rgb_out_dir)
    print("MS out directory:", multispec_out_dir)

    # Simple string paths for output files
    rgb_ortho_file = os.path.join(rgb_out_dir, f"{yyyymmdd}_{plot}_rgb_ortho.tif")
    ms_ortho_file = os.path.join(multispec_out_dir, f"{yyyymmdd}_{plot}_multispec_ortho.tif")
    print("RGB ortho file:", rgb_ortho_file)
    print("MS ortho file:", ms_ortho_file)

    # Orthomosaic file export compression
    compression = Metashape.ImageCompression()
    compression.tiff_compression = Metashape.ImageCompression.TiffCompressionNone  # disable LZW
    compression.jpeg_quality = 95  # set JPEG quality
    compression.tiff_big = True
    compression.tiff_tiled = True
    compression.tiff_overviews = True

    # Verify orthomosaic exists before exporting and export
    if merged_chunk.orthomosaic:
        try:
            print("Exporting RGB orthomosaic...")
            # Simple direct export
            merged_chunk.exportRaster(
                path=rgb_ortho_file,
                image_format=Metashape.ImageFormatTIFF,
                save_alpha=True,
                source_data=Metashape.OrthomosaicData,
                image_compression=compression,
                save_kml=False,
                save_world=False,
                global_profile=True,
                image_description="RGB orthomosaic",
                white_background=False,
                north_up=True,
                clip_to_boundary=False
            )
            print("Exported RGB orthomosaic:", rgb_ortho_file)
        except Exception as e:
            print(f"Error exporting RGB orthomosaic: {e}")
    else:
        print("RGB orthomosaic not available. Check previous processing steps.")

    # Verify duplicate orthomosaic exists before exporting
    if merged_duplicate.orthomosaic:
        try:
            print("Exporting multispectral orthomosaic...")
            
            # Simple direct export
            merged_duplicate.exportRaster(
                path=ms_ortho_file,
                image_format=Metashape.ImageFormatTIFF,
                raster_transform=Metashape.RasterTransformValue,
                save_alpha=False,
                source_data=Metashape.OrthomosaicData,
                image_compression=compression,
                save_kml=False,
                save_world=False,
                global_profile=True,
                image_description="Multispectral orthomosaic",
                white_background=False,
                north_up=True,
                clip_to_boundary=False
            )
            print("Exported multispectral orthomosaic:", ms_ortho_file)
        except Exception as e:
            print(f"Error exporting multispectral orthomosaic: {e}")
    else:
        print("Multispectral orthomosaic not available. Check previous processing steps.")

    # Remove orthophotos from orthomosaics
    print("Removing orthophotos from orthomosaics...")
    merged_chunk.orthomosaic.removeOrthophotos()
    merged_duplicate.orthomosaic.removeOrthophotos()
    doc.save()
    print("Orthophotos removed from orthomosaics.")

    print("Processing complete!")

    # Export PDF reports for both chunks
    print("Generating PDF reports...")
    
    # Create metadata directory if it doesn't exist
    metadata_dir = imagery_dir.parent / "metadata"
    os.makedirs(str(metadata_dir), exist_ok=True)
    print(f"Using metadata directory: {metadata_dir}")
    
    # RGB chunk report
    try:
        # Get RGB chunk (the first chunk, which contains only JPG cameras)
        rgb_chunk = doc.chunks[0]
        rgb_report_path = str(metadata_dir / f"{yyyymmdd}_{plot}_rgb_report.pdf")
        rgb_chunk.exportReport(
            path=rgb_report_path,
            title="RGB Project Report",
            description="Automated Metashape Processing - RGB Dataset",
            font_size=12,
            page_numbers=True,
            include_system_info=True,
            user_settings=[("Processing:", "DATEHERE")]
        )
        print(f"RGB PDF report exported to {rgb_report_path}")
    except Exception as e:
        print(f"Error exporting RGB PDF report: {e}")
    
    # Multispectral chunk report
    try:
        # Get multispectral chunk (the second chunk, which contains only TIF cameras)
        ms_chunk = doc.chunks[1]
        ms_report_path = str(metadata_dir / f"{yyyymmdd}_{plot}_multispec_report.pdf")
        ms_chunk.exportReport(
            path=ms_report_path,
            title="Multispectral Project Report",
            description="Automated Metashape Processing - Multispectral Dataset",
            font_size=12,
            page_numbers=True,
            include_system_info=True,
            user_settings=[("Survey Date", yyyymmdd), ("Processing:", "Juan C. Montes-Herrera")]
        )
        print(f"Multispectral PDF report exported to {ms_report_path}")
    except Exception as e:
        print(f"Error exporting multispectral PDF report: {e}")

if __name__ == "__main__":
    main()