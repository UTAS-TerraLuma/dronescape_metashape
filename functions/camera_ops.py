import datetime
import json
import subprocess
from pathlib import Path

import Metashape
import pandas as pd

def id_multispectral_camera(chunk):
    """
    Detects the multispectral camera bands and prints their layer indices:
    """
    print("Detecting multispectral camera.")
    
    # Get all multispectral sensors
    sensor_index = {}
    multispec_sensors = [s for s in chunk.sensors]
    for sensor in multispec_sensors:
        sensor_index[sensor.label] = sensor.layer_index
        print(sensor.label)
        print(sensor.layer_index)
    # print(sensor_index)

    # If no multispectral sensors found, return
    if not multispec_sensors:
        print("No multispectral sensors found.")
        return


def enable_oblique_cameras(chunk, enable=False):
    """
    Detects oblique cameras in the chunk and enables/disables them based on input.
    Targets cameras with absolute pitch angle < 46 degrees (more vertical/nadir cameras).
    
    Args:
        chunk: Metashape chunk containing cameras to process
        enable: Whether to enable oblique cameras (default: False - disables them)
    
    Returns:
        tuple: (oblique_count, oblique_cameras) - number of oblique cameras detected and list of oblique camera objects
    """
    oblique_count = 0
    oblique_cameras = []
    angle_threshold = 46.0  # Absolute angle threshold
    
    print(f"Detecting oblique cameras (absolute pitch threshold: {angle_threshold}°, enable: {enable})...")
    
    for camera in chunk.cameras:
        if not camera.enabled:
            continue
            
        if not hasattr(camera.photo, 'meta'):
            print(f"Warning: No metadata for camera {camera.label}")
            continue
            
        try:
            # Try to get gimbal pitch from EXIF data - different drones may use different tag names
            pitch = None
            
            # Check common EXIF tag names for gimbal pitch
            pitch_tags = [
                "DJI/GimbalPitchDegree", 
                "Gimbal/Pitch", 
                "EXIF:GimbalPitchDegree"
            ]
            
            for tag in pitch_tags:
                if tag in camera.photo.meta:
                    try:
                        pitch = float(camera.photo.meta[tag])
                        # print(f"Found pitch in tag: {tag} = {pitch:.1f}°")
                        break
                    except (ValueError, TypeError):
                        continue
            
            # If no pitch found in EXIF tags, try to calculate from rotation matrix if available
            if pitch is None and camera.transform:
                # Extract rotation matrix from camera transform
                rotation_matrix = camera.transform.rotation()
                yaw, pitch, roll = Metashape.Utils.mat2ypr(rotation_matrix)
                # print(f"Calculated pitch from transform matrix: {pitch:.1f}°")
            
            # If we still don't have pitch data, skip this camera
            if pitch is None:
                print(f"Warning: Could not determine pitch for camera {camera.label}")
                continue
                
            # Check if camera is oblique based on absolute pitch value
            if abs(pitch) < angle_threshold:
                oblique_cameras.append(camera)
                oblique_count += 1
                
                # Set camera enabled state based on enable parameter
                camera.enabled = enable
                status = "ENABLED" if enable else "DISABLED"
                # print(f"Oblique camera: {camera.label} (pitch: {pitch:.1f}°) - {status}")
                    
        except Exception as e:
            print(f"Warning: Could not process camera {camera.label}: {str(e)}")

    # print(f"Found {oblique_count} cameras with absolute pitch < {angle_threshold}°")
    action = "Enabled" if enable else "Disabled"
    print(f"{action} {oblique_count} cameras")
    
    return oblique_count, oblique_cameras


def filter_multispec(rgb_dir, ms_chunk, rgb_chunk, oblique_cameras=None):
    """
    Filter multispectral cameras to the RGB capture window using file timestamps.
    Returns the number of disabled multispectral cameras.
    """
    print("Running simple filtering of multispectral images...")

    def parse_exif_timestamp(value):
        if not value:
            return None
        try:
            return pd.to_datetime(str(value).replace(':', '-', 2))
        except Exception:
            return None

    def first_available_timestamp(metadata, fields):
        for field in fields:
            if metadata.get(field):
                timestamp = parse_exif_timestamp(metadata.get(field))
                if timestamp is not None:
                    return timestamp
        return None

    print(f"Extracting RGB timestamps from {rgb_dir}")
    rgb_result = subprocess.run(
        ['exiftool', '-j', '-r', '-ext', 'jpg', str(rgb_dir)],
        capture_output=True,
        text=True,
    )
    if rgb_result.returncode != 0 or not rgb_result.stdout.strip():
        print(f"Error extracting RGB timestamps: {rgb_result.stderr}")
        return 0

    rgb_json = json.loads(rgb_result.stdout)
    rgb_records = [
        {
            'filename': Path(m.get('SourceFile', '')).name,
            'path': m.get('SourceFile', ''),
            'timestamp': parse_exif_timestamp(m.get('UTCAtExposure')),
        }
        for m in rgb_json
        if m.get('UTCAtExposure')
    ]
    rgb_records = [r for r in rgb_records if r['timestamp'] is not None]
    if not rgb_records:
        print("No RGB timestamps found!")
        return 0

    rgb_df = pd.DataFrame(rgb_records)

    time_offset = pd.Timedelta(seconds=17.78)
    print("Applying hardcoded time offset: RGB timestamps are 17.78 seconds ahead of multispectral timestamps")

    first_rgb_time = rgb_df['timestamp'].min() - time_offset

    if not oblique_cameras:
        print("Detecting oblique cameras to find oblique phase start time...")
        _, oblique_cameras = enable_oblique_cameras(rgb_chunk, enable=False)
    else:
        print(f"Using {len(oblique_cameras)} pre-identified oblique cameras for timestamp filtering")

    oblique_start_time = None
    if oblique_cameras:
        oblique_filenames = {Path(camera.photo.path).name for camera in oblique_cameras}
        oblique_df = rgb_df[rgb_df['filename'].isin(oblique_filenames)]
        if not oblique_df.empty:
            oblique_start_time = oblique_df['timestamp'].min() - time_offset
            print(f"Using first (earliest) oblique RGB camera time as oblique phase start time: {oblique_start_time}")

    if oblique_start_time is None:
        oblique_start_time = rgb_df['timestamp'].max() - time_offset
        print("No oblique cameras found or could not determine timestamp. Using last RGB image time instead.")

    if oblique_start_time <= first_rgb_time:
        print("Warning: oblique start time is not after the first RGB time. Adjusting to use last RGB time instead.")
        oblique_start_time = rgb_df['timestamp'].max() - time_offset

    print(f"First RGB image time (adjusted): {first_rgb_time}")
    print(f"Oblique phase start time for multispectral filtering: {oblique_start_time}")

    camera_map = {Path(cam.photo.path).name: cam for cam in ms_chunk.cameras}
    for cam in camera_map.values():
        cam.enabled = True

    disabled_count = 0
    enabled_count = 0
    ms_min_timestamp = None
    ms_max_timestamp = None
    out_before_count = 0
    out_after_count = 0

    ms_paths = [cam.photo.path for cam in ms_chunk.cameras]
    batch_size = 20

    total_batches = (len(ms_paths) + batch_size - 1) // batch_size
    for batch_index in range(0, len(ms_paths), batch_size):
        batch_paths = ms_paths[batch_index: batch_index + batch_size]
        print(f"Processing batch {batch_index // batch_size + 1}/{total_batches}")

        batch_result = subprocess.run(
            ['exiftool', '-j'] + batch_paths,
            capture_output=True,
            text=True,
        )
        if batch_result.returncode != 0 or not batch_result.stdout.strip():
            print(f"Error in batch {batch_index // batch_size + 1}: {batch_result.stderr}")
            continue

        try:
            for metadata in json.loads(batch_result.stdout):
                filename = Path(metadata.get('SourceFile', '')).name
                timestamp = first_available_timestamp(
                    metadata,
                    (
                        'UTCAtExposure',
                        'SubSecModifyDate',
                    ),
                )
                if timestamp is None:
                    continue

                ms_min_timestamp = timestamp if ms_min_timestamp is None or timestamp < ms_min_timestamp else ms_min_timestamp
                ms_max_timestamp = timestamp if ms_max_timestamp is None or timestamp > ms_max_timestamp else ms_max_timestamp

                camera = camera_map.get(filename)
                if camera is None:
                    continue

                is_before_window = timestamp < first_rgb_time
                is_after_or_at_oblique = timestamp >= oblique_start_time

                if is_before_window:
                    out_before_count += 1
                    camera.enabled = False
                    disabled_count += 1
                elif is_after_or_at_oblique:
                    out_after_count += 1
                    camera.enabled = False
                    disabled_count += 1
                else:
                    camera.enabled = True
                    enabled_count += 1
        except Exception as exc:
            print(f"Error processing batch: {exc}")

    if ms_min_timestamp is not None and ms_max_timestamp is not None:
        print(f"Multispectral timestamps range from {ms_min_timestamp} to {ms_max_timestamp}")
    print(f"Disabled {disabled_count} multispectral images outside RGB time range")
    print(f" - Before first RGB time: {out_before_count}")
    print(f" - At/after oblique start time: {out_after_count}")
    print(f"Kept {enabled_count} multispectral images within RGB time range")
    return disabled_count
