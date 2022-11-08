# Readme for NPUTE

## REQUIREMENTS

 - Python (2.5-2.7, NOT version 3.0 --> update: it now also works with Python >= 3)
 - Numpy (numpy.scipy.org)

## COMMAND-LINE OPTIONS


| Option | Identifier | Parameters |
| - |  - |  - |
| Mode | `-m` | 0 for impute, 1 for test |
| Single Window | `-w` | A positive integer for the window size |
| Window File | `-W` | The name of file containing a list of window sizes separated by newline characters. |
| Window Range | `-r` | StartRange:EndRange |
| Input File | `-i` | The name of a comma-separated file containing the input data. ï¿½Rows are SNPs, columns are haplotypes. |
| Output File | `-o` | The name of the file to output to. ï¿½DATA WILL BE OVERWRITTEN. |

## TUTORIAL

We recommend that you use NPUTE by first testing a large range of windows and using that with the highest estimated accuracy to do the actual imputation.ï¿½ In this tutorial, we will guide you through the process of imputing the included sample data (sample_data.csv).

Once you have installed Python and Numpy, you may go to the Start Menu > Run*, type in cmd, and press enter.ï¿½ This should bring up the command-line interface.ï¿½ Using DOS commands (if in Windows), navigate to the directory containing sample_data.csv.ï¿½ Entering the following command will run an imputation test on the data for a range of windows between 5 and 30 and output the results in the file out.csv.

    python NPUTE.py -m 1 -r 5:30 -i sample_data.csv -o out.csv

After 5-15 minutes, the process will be complete and the output will be stored in out.csv.ï¿½ When you open the file in Excel, you'll see that the highest estimated accuracy is for a window size of 12 (97.11% correct).

Now to do the actual imputation at that window size, use the following command:

    python NPUTE.py -m 0 -w 12 -i sample_data.csv -o imputed_data_w12.csv

The unknowns will be imputed and all of the data will be written to imputed_data_w12.csv with imputed values in lower case.

_* If using Mac OS or Linux, use a terminal window (ex. xterm) instead._

## CONTACT

If you have any questions or comments, please feel free to conact us at ajr@unc.edu or mcmillan@cs.unc.edu.

## REFERENCE

> Roberts, Adam , Leonard McMillan, Wei Wang, Joel Parker, Ivan Rusyn, and David Threadgill.<br />
> "Inferring missing genotypes in large SNP panels using fast nearest-neighbor searches over sliding windows."<br />
> Bioinformatics (ISMB/ECCB 2007 Conference), 23(13), pp. i401-l407, July, 2007. [doi:10.1093/bioinformatics/btm220](https://doi.org/10.1093/bioinformatics/btm220)

© 2007, 2010, 2012, Adam Roberts, Leonard McMillan, and and UNC Computational Genetics
