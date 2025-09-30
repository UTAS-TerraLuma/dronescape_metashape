## dronescape_metashape: RGB and Multispectral imagery processing in Agisoft Metashape Pro

These scripts were designed for Agisoft Metashape Professional (V.1.1) to generate spatially co-aligned RGB (DJI Zenmuse-P) and 10-band multispectral (MicaSense Red Edge-P Dual) orthomosaics.

The steps outlined here support the creation of analysis-ready data products for TERN's ecosystem surveillance program. Procedures for deploying uncrewed aerial vehicle (UAV) and collecting simultaneous imagery are detailed in TERN's UAV operational guidelines. For more information please refer to the [Drone Data Collection Protocol](https://www.tern.org.au/field-survey-apps-and-protocols/). Georeferencing is achieved using the DJI D-RTK 2 base station, which integrates Global Navigation Satellite System Real-Time Kinematic (GNSS RTK) to provide centimeter-level positioning accuracy for the UAV. This setup ensures precise georeferencing without the need for ground control points.

Detailed instructions on project settings, image selection and alignment, model creation, reflectance conversion, and file export are provided. A key feature of this workflow is the alignment of imagery from both cameras into a unified model used for the orthomosaics surface. Only tie points are used for model creation, reducing computational power and time requirements. The resulting unified model is used to create two separate products: an RGB orthomosaic and a 10-band multispectral orthomosaic. For a summary of the processing steps and information on Agisoft Metashape setup, please refer to the [Drone RGB and Multispectral Processing Protocol](https://www.tern.org.au/field-survey-apps-and-protocols/).

**Funding**: This project was funded by TERN Surveillance  
**Authors**: J.C. Montes-Herrera, Poornima Sivanandam, Alice Robbins, Arko Lucieer, Darren Turner.
School of Geography, Planning and Spatial Sciences, University of Tasmania  
**Acknowledgements**: Ben Sparrow (TERN Surveillance)
