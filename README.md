# HIRONEX: Historical road network extractor
HIRONEX is a python tool for automatic, fully unsupervised extraction of historical road networks from historical maps. 
HIRONEX reads a shapefile containing the modern road network, and a scanned & georeferenced historical map from year T, covering the same area, in the same spatial reference system.
Then, HIRONEX samples the color information from the historical map, and calculates a probability for each modern road segment of overlapping with a road symbol in the historical map. Then, 1d cluster analysis is applied to extract road segments that are likely to have existed in the year T of the historical map

The HIRONEX output is not perfect, as HIRONEX is work in progress and is currently being improved.
