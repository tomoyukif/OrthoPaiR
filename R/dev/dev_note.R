2025-06-13
Pairs counted as false positive in OrthoPaiR results are likely to be reasonable based on locations and sequence similarity.
These "false positive" are rather detections of orthologous genes that were largely disrupted in one genome.

The next tasks are:
  Modify eval_busco_list() to give penalty on pairing of MtoM pairs if all combinations of MtoM pairs were identified.
  Build a new function to summarize differences between the list of orthologous pairs detected by the tools.
  Build a new function assess the validity of differently detected pairs and also false positive and false negative detections.
    
2025-06-16
I have modified eval_busco_list() to omit busco pairs that shows translocations
(inter-chromosomal pairing and 1 Mb < location difference of paired genes on genomes for non-1to1 pairs).
I found that still I need verify the pairs found in OrthoPaiR. 
Thus, I do the followings.
  Build a new function to summarize differences between the list of orthologous pairs detected by the tools.
  Build a new function assess the validity of differently detected pairs and also false positive and false negative detections.
  
2025-06-17
I developed a function to summarize differences between the list of orthologous pairs detected by the tools.
Next task is to build a new function assess the validity of differently detected pairs and also false positive and false negative detections. How to assess the validity?
    
2025-07-08
I found that I should revise the rbh() function that contains incorrect data handling to filter RBBHs.
I also started to simplify the data processing in the rbh() function by removing unnecessary steps.
The next step is to revise RBBH filterin process especially when genes have multiple hits with identical scores.