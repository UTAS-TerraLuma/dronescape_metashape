"""
TERN Dronescape Metashape Processing Package
"""

# Expose key functions from modules
from .camera_ops import detect_multispectral_camera, log_camera_timestamps, enable_oblique_cameras, remove_images_outside_rgb_times
from .utils import find_images, extract_timestamp_from_exif, analyze_time_difference, get_common_altitude_range
from .time_sync import detect_time_offset, synchronize_cameras, filter_cameras_by_time_alignment
from .gpu_setup import setup_gpu
from .processing import detect_reflectance_panels 