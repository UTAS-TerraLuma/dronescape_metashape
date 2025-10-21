#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to initialize a Metashape project with RGB and multispectral images from Mavic 3M.

Author: Juan C. Montes Herrera (University of Tasmania)

Assumes directory structure:
    <plot>/YYYYMMDD/m3m/level0_raw/ (contains both .jpg and .tif files)

User provides:
    --imagery_dir: path to m3m directory (e.g., <plot>/YYYYMMDD/m3m/)
    --crs: EPSG code for target CRS (optional, defaults to 4326)
    --out: Metashape project location, not orthomosaics
    --enable_oblique: flag to enable oblique cameras (default: disabled)
Project will be named as "YYYYMMDD-plot-m3m.psx"

Example usage:
    # Basic usage with required arguments
    -imagery_dir /path/to/SITE-01/20230615/m3m/ -out /path/to/output/

    # With custom CRS
    -imagery_dir /path/to/SITE-01/20230615/m3m/ -out /path/to/output/ -crs 3577

    # Enable oblique cameras
    -imagery_dir /path/to/SITE-01/20230615/m3m/ -out /path/to/output/ -enable_oblique
"""

import argparse
import os
import sys
from pathlib import Path

import Metashape

from functions.gpu_setup import setup_gpu
from functions.utils import find_filtered_images
from functions.camera_ops import enable_oblique_cameras

## SfM Processing Parameters
downscale_align = 4 # Use 4 for tests

keypoint_limit = int(50000)
tiepoint_limit = int(5000)

def main():
    # Set up GPU acceleration
    setup_gpu()
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Initialize Metashape project with RGB and multispectral images from Mavic 3M.")
    parser.add_argument('-imagery_dir', required=True, help='Path to m3m directory (<plot>/YYYYMMDD/m3m/)')
    parser.add_argument('-crs', default="4326", help='EPSG code for target projected CRS (default: 4326, WGS84)')
    parser.add_argument('-out', required=True, help='Directory to save the Metashape project')
    parser.add_argument('-enable_oblique', action='store_true', help='Enable oblique cameras (default: disabled)')
    args = parser.parse_args()

    # Extract YYYYMMDD and plot from input path
    m3m_dir = Path(args.imagery_dir).resolve()
    yyyymmdd = m3m_dir.parent.name
    plot = m3m_dir.parent.parent.name
    print(f"Processing plot: {plot}, date: {yyyymmdd}")
    project_name = f"{yyyymmdd}-{plot}-m3m.psx"

    # Set up path for Mavic 3M imagery (contains both JPG and TIF files)
    level0_raw_dir = m3m_dir / "level0_raw"

    if not level0_raw_dir.is_dir():
        print(f"level0_raw directory not found: {level0_raw_dir}")
        sys.exit(1)

    # Find RGB images (jpg files) and multispectral images (tif files)
    rgb_images = find_filtered_images(level0_raw_dir, extensions=('.jpg', '.jpeg'))
    multispec_images = find_filtered_images(level0_raw_dir, extensions=('.tif', '.tiff'), exclude_patterns=())

    if not rgb_images:
        print(f"No RGB images found in {level0_raw_dir}")
    if not multispec_images:
        print(f"No multispectral images found in {level0_raw_dir}")
        
    print(f"Found {len(rgb_images)} RGB images")
    print(f"Found {len(multispec_images)} multispectral images")
    
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
        
    print("\nPhotos loaded successfully into respective chunks.")

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
        guided_matching=False,  # Disable guided image matching,
    )
    
    # Align cameras
    merged_chunk.alignCameras()
    
    # Optimize cameras with specified parameters
    merged_chunk.optimizeCameras(
        fit_f=True,  # Fit focal length
        fit_cx=True, fit_cy=True, # Fit principal point x, y
        fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=False, # Fit radial distortion k1, k2, k3 
        fit_p1=True, fit_p2=True, # Fit tangential distortion p1, p2
        fit_b1=True, fit_b2=True, # Fit affinity b1, b2
        fit_corrections=True,  # Fit additional corrections
        )
    
    print("Image alignment complete!")

    # print("Building depth maps...") # TODO: Add to research tests
    # merged_chunk.buildDepthMaps(
    #     downscale=4,  # Medium quality (use 1 or 0 for higher detail)
    #     filter_mode=Metashape.MildFiltering  # Options: NoFiltering, MildFiltering, ModerateFiltering, AggressiveFiltering
    #     )

    # Build dense point cloud
    # merged_chunk.buildPointCloud(
    #     source_data=Metashape.DepthMapsData,  # Extract dense cloud from depth maps
    #     point_colors=False,                    # Store colors for points
    #     point_confidence=False,               # Skip confidence values (set True if needed)
    #     keep_depth=True,                      # Keep depth maps for future use
    #     max_neighbors=100,                    # Limit neighbor images
    #     uniform_sampling=True,                # More even density distribution
    #     points_spacing=0.1                    # Desired point spacing in meters
    # )

    # Build the model
    merged_chunk.buildModel(
        surface_type=Metashape.Arbitrary,
        source_data=Metashape.TiePointsData, # TODO: Add to research tests
        # source_data=Metashape.DepthMapsData, # TODO: Add to research tests
        face_count=Metashape.MediumFaceCount, # TODO: Add to research tests
        interpolation=Metashape.EnabledInterpolation, # TODO: Add to research tests
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

    # Build orthomosaic with explicit projection
    merged_chunk.buildOrthomosaic(
        surface_data=Metashape.DataSource.ModelData,
        blending_mode=Metashape.MosaicBlending,
        refine_seamlines=True,
    )
    doc.save() # check blend mode and parameters for Ortho
    
    # Build orthomosaic with explicit projection for duplicate chunk
    merged_duplicate.buildOrthomosaic(
        surface_data=Metashape.DataSource.ModelData,
        blending_mode=Metashape.MosaicBlending,
        refine_seamlines=True
    )
    doc.save()
    print("RGB ortho:", merged_chunk.orthomosaic)
    print("MS ortho:", merged_duplicate.orthomosaic)

    rgb_res = round(merged_chunk.orthomosaic.resolution, 2)
    ms_res = round(merged_duplicate.orthomosaic.resolution, 2)
    print("RGB ortho resolution:", rgb_res)
    print("MS ortho resolution:", ms_res)


    # Set up output path for both RGB and multispectral orthomosaics
    out_dir = os.path.join(str(m3m_dir), "level1_proc")
    
    # Create output directory if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)
    
    print("Output directory:", out_dir)

    # Simple string paths for output files
    rgb_ortho_file = os.path.join(out_dir, f"{yyyymmdd}_{plot}_m3m_rgb_ortho.tif")
    ms_ortho_file = os.path.join(out_dir, f"{yyyymmdd}_{plot}_m3m_multispec_ortho.tif")
    print("RGB ortho file:", rgb_ortho_file)
    print("MS ortho file:", ms_ortho_file)

    compression = Metashape.ImageCompression()
    compression.tiff_compression = Metashape.ImageCompression.TiffCompressionNone  # disable LZW
    compression.jpeg_quality = 95  # set JPEG quality
    compression.tiff_big = True
    compression.tiff_tiled = True
    compression.tiff_overviews = True

    # Verify orthomosaic exists before exporting
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

    # Verify     duplicate orthomosaic exists before exporting
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
    metadata_dir = m3m_dir.parent / "metadata"
    os.makedirs(str(metadata_dir), exist_ok=True)
    print(f"Using metadata directory: {metadata_dir}")
    
    # RGB chunk report
    try:
        # Get RGB chunk (the first chunk, which contains only JPG cameras)
        rgb_chunk = doc.chunks[0]
        rgb_report_path = str(metadata_dir / f"{yyyymmdd}_{plot}_m3m_rgb_report.pdf")
        rgb_chunk.exportReport(
            path=rgb_report_path,
            title="RGB Project Report",
            description="Automated Metashape Processing - RGB Dataset (Mavic 3 Multispectral)",
            font_size=12,
            page_numbers=True,
            include_system_info=True,
            user_settings=[("Processing:", "Alice Robbins")]
        )
        print(f"RGB PDF report exported to {rgb_report_path}")
    except Exception as e:
        print(f"Error exporting RGB PDF report: {e}")
    
    # Multispectral chunk report
    try:
        # Get multispectral chunk (the second chunk, which contains only TIF cameras)
        ms_chunk = doc.chunks[1]
        ms_report_path = str(metadata_dir / f"{yyyymmdd}_{plot}_m3m_multispec_report.pdf")
        ms_chunk.exportReport(
            path=ms_report_path,
            title="Multispectral Project Report",
            description="Automated Metashape Processing - Multispectral Dataset (Mavic 3 Multispectral)",
            font_size=12,
            page_numbers=True,
            include_system_info=True,
            user_settings=[("Survey Date", yyyymmdd), ("Processing:", "Alice Robbins")]
        )
        print(f"Multispectral PDF report exported to {ms_report_path}")
    except Exception as e:
        print(f"Error exporting multispectral PDF report: {e}")

if __name__ == "__main__":
    main()