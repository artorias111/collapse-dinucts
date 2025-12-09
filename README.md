["AC", "AG", "AT", 
                    "CA", "CG", "CT", 
                    "GA", "GC", "GT", 
                    "TA", "TC", "TG"]


#### What are the flags?
- `-reads`: Path to reads (as a `fastq` or `fastq.gz`)
- `threshold`: If set to `0.0`, collapses dinucleotides. If set to any number `n > 0`, returns reads that consist of dinucleotides only below `n%` of the total read length. Default set to `0.0`
