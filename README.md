# SV_bedgr
Script takes N vcf files, splits them into 3*N BedGraph files with duplication, deletion and inversion. Then it unite them in 3 files and calculates frequency of SV in each interval. Finally, it merges these files. 
BedGraph files after splitting are stored in separate folder (bedgr_dir). Other files are stored in another folder (output).
In script three arguments are required: vcfs_dir - path to folder containing vcf files, bedgr_dir - path to folder containing BedGraph files after splitting, bedtools - path to bedtools exec, output - path to output folder.
Optional arguments: quality - intervals with quality above X will be used, taken from vcf files, default = 0, filter_thresh - intervals with frequency above X will be used, default = 0.01, threads - number of parallel threads, default = 50.
Example: $ python3 rotation.py ~/vcf_dir ~/bedgr_dir /path to bedtools ~/output_dir --quality=100 --filter_thresh=0.1
