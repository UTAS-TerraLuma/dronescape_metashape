#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to initialize a Metashape project with RGB and multispectral images in separate chunks.
Assumes TERN directory structure:
    <plot>/YYYYMMDD/imagery/
        ├── rgb/level0_raw/
        └── multispec/level0_raw/
User provides:
    --imagery_dir: path to YYYYMMDD/imagery/
    --crs: EPSG code for target CRS (optional, defaults to 4326)
    --out: output directory for Metashape project
    --enable_oblique: flag to enable oblique cameras (default: disabled)
Project will be named as "YYYYMMDD-plot.psx"

Example usage:
    # Basic usage with required arguments
    -imagery_dir /path/to/SITE-01/20230615/imagery/ -out /path/to/output/

    # With custom CRS
    -imagery_dir /path/to/SITE-01/20230615/imagery/ -out /path/to/output/ -crs 3577

    # Enable oblique cameras
    -imagery_dir /path/to/SITE-01/20230615/imagery/ -out /path/to/output/ -enable_oblique

    # Run diagnostics only without creating project
    -imagery_dir /path/to/SITE-01/20230615/imagery/ -out /path/to/output/ -diagnose_only
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
# from functions.processing import DICT_SMOOTH_STRENGTH

# def resume_proc():
    # # Calibrate reflectance 
    # multispec_chunk.calibrateReflectance(use_reflectance_panels=True, use_sun_sensor=use_sun_sensor)

    # # Raster transform multispectral images
    # print("Updating Raster Transform for relative reflectance")
    # raster_transform_formula = []
    # num_bands = len(multispec_chunk.sensors)
    # for band in range(1, num_bands+1):
    #     raster_transform_formula.append("B" + str(band) + "/32768")

    # chunk.raster_transform.formula = raster_transform_formula
    # chunk.raster_transform.calibrateRange()
    # chunk.raster_transform.enabled = True
    # doc.save()
    # print(f"Applied raster transform formulas: {raster_transform_formula}")

def main():
    # Set up GPU acceleration
    setup_gpu()
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Initialize Metashape project with RGB and multispectral images.")
    parser.add_argument('-imagery_dir', required=True, help='Path to YYYYMMDD/imagery/ directory')
    parser.add_argument('-crs', default="4326", help='EPSG code for target projected CRS (default: 4326, WGS84)')
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

    print("Aligning images...")
    # Match photos with specified settings
    merged_chunk.matchPhotos(
        # downscale=1,  # High accuracy
        downscale=8, # Lowest accuracy
        generic_preselection=True,  # Enable generic preselection
        reference_preselection=True,  # Enable reference preselection
        reference_preselection_mode=Metashape.ReferencePreselectionSource,  # Source mode
        keypoint_limit=50000,  # Key points limit
        tiepoint_limit=5000,  # Tie points limit
        filter_stationary_points=True,  # Exclude stationary points
        guided_matching=False,  # Disable guided image matching,
    )
    
    # Align cameras
    merged_chunk.alignCameras()
    
    # Optimize cameras with specified parameters
    merged_chunk.optimizeCameras(
        fit_f=True,  # Fit focal length
        fit_cx=True,  # Fit principal point x, y
        fit_cy=True,
        fit_k1=True,  # Fit radial distortion k1, k2, k3
        fit_k2=True,
        fit_k3=True,
        fit_k4=False, # False
        fit_p1=True,  # Fit tangential distortion p1, p2
        fit_p2=True,
        fit_b1=True,  # Fit affinity b1, b2
        fit_b2=True,
        fit_corrections=True,  # Fit additional corrections
    )
    
    print("Image alignment complete!")

    # Build the model
    merged_chunk.buildModel(
        surface_type=Metashape.HeightField,
        source_data=Metashape.TiePointsData,
        face_count=Metashape.MediumFaceCount,
        interpolation=Metashape.EnabledInterpolation,
        build_texture=False,
        vertex_colors=False
    )
    
    # Smooth model based on specified strength
    # print(f"Smoothing model with {smooth_strength} strength...")
    # smooth_val = DICT_SMOOTH_STRENGTH[smooth_strength]
    # merged_chunk.smoothModel(smooth_val, fix_borders=True)
    
    print("Model building complete!")

    merged_duplicate = merged_chunk.copy()
    merged_duplicate.label = "merged_duplicate"
    doc.chunks.append(merged_duplicate)
    doc.save()

    # # Set CRS for both chunks
    # crs_code = args.crs
    # target_crs = Metashape.CoordinateSystem(f"EPSG::{crs_code}")
    # rgb_chunk.crs = target_crs
    # multispec_chunk.crs = target_crs

    # Remove TIF cameras from merged_chunk
    tif_cams = [cam for cam in merged_chunk.cameras if cam.photo.path.lower().endswith(".tif")]
    if tif_cams:
        merged_chunk.remove(tif_cams)

    doc.save()

    # Remove JPG cameras from merged_duplicate
    jpg_cams = [cam for cam in merged_duplicate.cameras if cam.photo.path.lower().endswith(".jpg")]
    if jpg_cams:
        merged_duplicate.remove(jpg_cams)

    doc.save()

    print("Project setup complete!")
    print("###########################")
    print("###########################")

     # Calibrate reflectance 
    merged_duplicate.calibrateReflectance(use_reflectance_panels=True, use_sun_sensor=False)

    # Raster transform multispectral images
    print("Updating Raster Transform for relative reflectance")
    raster_transform_formula = []
    num_bands = len(merged_duplicate.sensors)
    for band in range(1, num_bands+1):
        raster_transform_formula.append("B" + str(band) + "/32768")

    merged_duplicate.raster_transform.formula = raster_transform_formula
    merged_duplicate.raster_transform.calibrateRange()
    merged_duplicate.raster_transform.enabled = True
    doc.save()
    print(f"Applied raster transform formulas: {raster_transform_formula}")

    print("Building RGB orthomosaic...")
    try:
        # Check if model exists (should exist from earlier in the script)
        if not merged_chunk.model:
            print("WARNING: Model not found for RGB chunk - this should not happen")
        
        # Build orthomosaic using model data as requested
        merged_chunk.buildOrthomosaic(
            surface_data=Metashape.ModelData,  # Use model as surface
            refine_seamlines=True, 
            blending_mode=Metashape.MosaicBlending
        )
        doc.save()
        
        # Verify orthomosaic was created successfully
        if merged_chunk.orthomosaic and merged_chunk.orthomosaic.key:
            print("RGB orthomosaic built successfully. Details:")
            print(f"Projection: {merged_chunk.orthomosaic.projection}")
            print(f"Resolution: {merged_chunk.orthomosaic.resolution}")
            print(f"Size (pixels): {merged_chunk.orthomosaic.width} x {merged_chunk.orthomosaic.height}")
        else:
            print("WARNING: RGB orthomosaic build may have failed - object exists but seems invalid")
    except Exception as e:
        print(f"ERROR building RGB orthomosaic: {str(e)}")

    print("Building multispectral orthomosaic...")
    try:
        # Check if model exists (should exist from earlier in the script)
        if not merged_duplicate.model:
            print("WARNING: Model not found for multispectral chunk - this should not happen")
            
        # Build orthomosaic using model data as requested
        merged_duplicate.buildOrthomosaic(
            surface_data=Metashape.ModelData,  # Use model as surface
            refine_seamlines=True, 
            blending_mode=Metashape.MosaicBlending
        )
        doc.save()
        
        # Verify orthomosaic was created successfully
        if merged_duplicate.orthomosaic and merged_duplicate.orthomosaic.key:
            print("Multispectral orthomosaic built successfully. Details:")
            print(f"Projection: {merged_duplicate.orthomosaic.projection}")
            print(f"Resolution: {merged_duplicate.orthomosaic.resolution}")
            print(f"Size (pixels): {merged_duplicate.orthomosaic.width} x {merged_duplicate.orthomosaic.height}")
        else:
            print("WARNING: Multispectral orthomosaic build may have failed - object exists but seems invalid")
    except Exception as e:
        print(f"ERROR building multispectral orthomosaic: {str(e)}")

    # Set up output paths for RGB and multispectral orthomosaics
    rgb_out = imagery_dir / "rgb" / "level1_proc"
    multispec_out = imagery_dir / "multispec" / "level1_proc"
    
    # Ensure output directories exist
    rgb_out.mkdir(parents=True, exist_ok=True)
    multispec_out.mkdir(parents=True, exist_ok=True)

    rgb_res = round(merged_chunk.orthomosaic.resolution, 2)
    ms_res = round(merged_duplicate.orthomosaic.resolution, 2)
    print("RGB ortho resolution:", rgb_res)
    print("MS ortho resolution:", ms_res)

    # file naming format: <projname>_multispec_ortho_<res_in_m>.tif
    rgb_ortho_file = rgb_out / (yyyymmdd + "_" + plot + "_rgb_ortho.tif")
    ms_ortho_file = multispec_out / (yyyymmdd + "_" + plot + "_multispec_ortho.tif")
    print("RGB ortho file:", rgb_ortho_file)
    print("MS ortho file:", ms_ortho_file)

    compression = Metashape.ImageCompression()
    compression.tiff_compression = Metashape.ImageCompression.TiffCompressionNone  # disable LZW
    # compression.jpeg_quality = 95  # set JPEG quality
    compression.tiff_big = True
    compression.tiff_tiled = True
    compression.tiff_overviews = True

    # Make the RGB chunk the active chunk before export
    doc.chunk = merged_chunk
    
    print("RGB chunk activated, exporting...")
    try:
        # Get the active chunk directly like in your working script
        active_chunk = Metashape.app.document.chunk
        
        active_chunk.exportRaster(
            path=str(rgb_ortho_file), 
            resolution_x=rgb_res, 
            resolution_y=rgb_res,
            image_format=Metashape.ImageFormatTIFF,
            save_alpha=False, 
            source_data=Metashape.OrthomosaicData, 
            image_compression=compression,
            save_world=False, 
            save_kml=False, 
            image_description="RGB orthomosaic",
            white_background=False, 
            north_up=True
        )
        print("Exported RGB orthomosaic: " + str(rgb_ortho_file))
    except Exception as e:
        print(f"ERROR exporting RGB orthomosaic: {str(e)}")

    # Make the multispectral chunk the active chunk before export
    doc.chunk = merged_duplicate
    
    print("Multispectral chunk activated, exporting...")
    try:
        # Get the active chunk directly like in your working script
        active_chunk = Metashape.app.document.chunk
        
        active_chunk.exportRaster(
            path=str(ms_ortho_file), 
            resolution_x=ms_res, 
            resolution_y=ms_res,
            image_format=Metashape.ImageFormatTIFF,
            raster_transform=Metashape.RasterTransformValue,
            save_alpha=False, 
            source_data=Metashape.OrthomosaicData, 
            image_compression=compression,
            save_world=False, 
            save_kml=False, 
            image_description="Multispectral orthomosaic",
            white_background=False, 
            north_up=True
        )
        print("Exported multispec orthomosaic: " + str(ms_ortho_file))
    except Exception as e:
        print(f"ERROR exporting multispectral orthomosaic: {str(e)}")

    print("Processing complete!")
    doc.save()

if __name__ == "__main__":
    main()