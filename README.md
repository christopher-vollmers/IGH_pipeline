# IGH_pipeline

Takes MiSeq 2x300 paired end reads with 14nt UMIs and merges them based on their UMIs. 
Old code but functional. Can be adapted to other amplicons but pResto might be a better bet if you need a flexible read merger. 
Uploaded mostly for archival reasons. 

Needs FLASH read merger in path 


How to run:

python3 IGH_pipeline.py file_containing_path_to_folders_containing_R1_and_R2_files.
One directory per line. 
