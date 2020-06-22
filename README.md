# tea
----------------------------------------------------------------------------------------------------------------------------------------
                                                  Target Embedding Algorithm
----------------------------------------------------------------------------------------------------------------------------------------
The Target Embedding Algorithm (T.E.A.) is a method to automatically insert target letters based on differences in contrast between the letter and the image. This allows the experimenter to adjust the perceived difficulty of a target location (e.g. the median contrast change will yield a 'medium difficulty' location). It also allows for targets of varying sizes to keep the same target location due to the TEA summing each contrast difference map, and operating on that.

Please see Clayden, Fisher, Nuthmann (in review) or Clayden (2019): https://era.ed.ac.uk/handle/1842/35932 for an overview.

**How to use:**

The algorithm currently supports up to 4 target sizes, with a number of variables for you to adjust in the TargetEmbeddingAlgorithm() function

These are:

* imagePath - path where your images are stored

* exportPath - path you wish to export final images to

* img_info - path containing a list of image names from (e.g. .xls file)

* columnName - name of the column to read image names from

* image width and height (e.g. 800x600) a) desiredImageWidth b) desiredImageHeight

* monitor values: x_monitor (visible monitor width (cm)) y_monitor (visible monitor height (cm))
    
* distance between monitor and subject (cm): MonitorSubjectDistance

* radius of the central circular exclusion mask in degrees: radius_deg

* size of a region buffer to place around the target to define the region size: region_buffer
    
* thickness of the letter
*These are example values that you can change*
horizontalWidth = [1 2 2 4];
verticalWidth = [1 3 3 5];
    
* actual width and height
*These are example values that you can change*
widthVector = [7 13 19 33];
heightVector = [9 16 23 40];
    
* the index of each target size, with the length being the count
letterMatrix = [1 2 3 4];
    
* valueType - supports: Min, Lower Quartile, Median, Upper Quartile, Max
    
* placeWithTypeThreshold - Change this to find suitable locations based on a certain target size. Largest size is recommended and is defined by the number within letterMatrix. E.g. 1, 2, 3, etc

After you have entered the desired values, the algorithm will then:

* Prepare the central exclusion mask
* Prepare the images based on the import path by converting them all to the desired width/height and to 
  greyscale. 
* Generate the contrast difference maps (CDMs)
* Generate the summed difference map based on the CDMs for each target size
* Obtain location coordinates, embed target and write final images to export path
