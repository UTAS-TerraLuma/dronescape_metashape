import Metashape

chunk = Metashape.app.document.chunk

# Set up TIFF compression
compression = Metashape.ImageCompression()
compression.tiff_compression = Metashape.ImageCompression.TiffCompressionNone  # LZW compression
compression.tiff_big = True        # Enable BigTIFF
compression.tiff_tiled = True      # Enable tiled TIFF
compression.tiff_overviews = True  # Generate overviews (pyramids)
compression.jpeg_quality = 95

# Export orthomosaic with compression and metadata options
chunk.exportRaster(
    path=r"C:\Users\jcmontes\Desktop/newband_script2.tif",
    image_format=Metashape.ImageFormatTIFF,
    raster_transform=Metashape.RasterTransformValue,
    save_alpha=False,
    source_data=Metashape.OrthomosaicData,
    image_compression=compression,
    save_world=False,      # Save accompanying world file
    save_kml=False,        # Save KML file for geolocation
    image_description="Orthomosaic exported from Metashape with full TIFF metadata",
    white_background=False,
    north_up=True
)

# # Export camera metadata (positions, orientation, etc.)
# chunk.exportCameras(
#     path=str(cameras_file),
#     format=Metashape.CamerasFormatXML,
#     save_points=True,
#     save_markers=True,
#     use_labels=True
# )

print(f"Exported orthomosaic")
# print(f"Exported camera metadata: {cameras_file}")

