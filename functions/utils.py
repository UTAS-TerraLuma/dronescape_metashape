import os
import datetime
import exifread

def find_filtered_images(folder, extensions=(), exclude_patterns=()):
    """
    Helper function to recursively find image files with given extensions,
    excluding files that match any patterns in exclude_patterns.
    
    Args:
        folder: Path to search for images
        extensions: Tuple of allowed file extensions
        exclude_patterns: Tuple of patterns to exclude (file endings)
        
    Returns:
        List of image file paths that match the criteria
    """
    image_list = []
    for root, _, files in os.walk(folder):
        for fname in files:
            # Check if file has allowed extension
            if extensions and not fname.lower().endswith(extensions):
                continue
                
            # Check if file matches any exclude pattern
            should_exclude = False
            for pattern in exclude_patterns:
                if fname.endswith(pattern):
                    should_exclude = True
                    break
                    
            if not should_exclude:
                image_list.append(os.path.join(root, fname))
    return image_list


# Constants
DICT_SMOOTH_STRENGTH = {
    'low': 50,      # For low-lying vegetation (grasslands, shrublands)
    'medium': 100,  # For mixed vegetation
    'high': 200     # For forested sites
}