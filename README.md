readme for github gg cloning project

This is a two-tool program that 1) searches if the provided 4nt sequence (i.e. type iis re that leaves a 4nt overhang sequence) is compatible with one another in an in-silico cloning design and 2) creates a set of compatible 4nt overhang sequence dependent on a sequence of 2AA sequences provided. The latter is useful for generating compatible 4nt sequences for junction sequences that one desires to dictate the AA used in the junction sites.

Below is a common list of type iis re that leaves a 4nt overhang sequence:
BsaI
BbsI
BsmbI

This program is taylored for type iis re that leaves a 4nt overhang sequence. There are plans on furthering the capability of this program by encompassing type iis re that leave 3nt overhang sequences and overhang sequences >4nt in length.

Download the code, or copy the code into your favorite text editor. Please check the directory written in this code - this is the only feature that must be changed to successfully run the code. There is a directory for your CSV, JSON, and python file. Please set var csv_dir and json_dir to your CSV and JSON directories, respectively.

CSV: please format your CSV as the following. The column headers is not case sensitive, however col 1 must be names 'NT sequence' (for NT_checker.py script) and col 2 must be named 'AA sequence' (for AA_to_NT_generator.py script). In your CSV file, each row corrensponds to individual sequences (whether NT or AA) that is to be checked. Please list all NT sequences in column one. The case is not sensitive, but make sure to write a 4nt sequence per row. At least two rows are required, or an exception will be thrown.
For the AA column, please write your AA sequence as two single-letter AA code. Two and only two AA sequence must be written per row. The cells are not case sensitive. More than two rows are required to run the program, or an execption will be thrown. (please reference figure below)

JSON: please download the JSON file 'gg_ntdict.json.' Do not change the name of this file or an exception will be thrown. You are welcome to change the name of your file, but please make sure to change every call to the file in the program to your new name. Store this file appropriately in your JSON directory.

Running the program: please run the program in your terminal. Change directories to your script location. Run AA to NT junction sequence generator script by typing 'python AA_to_NT_generator.py' in your terminal. Run the NT junction site checker script by typing 'python NT_checker.py' in your terminal (of course in a directory that has the script being called).

NT_checker.py: This script will look at the content in col1 (where only NT sequences are listed). A detailed result list will print on your terminal. If the NT sequence fails, reasons for failure will be delineated on the terminal. 

AA_to_NT_generator.py: This script will look at the content in col2 (where only AA sequences are listed). A load bar and status of your script will be printed on your terminal. If script is completed successfully, you will find a summary .csv file written to your CSV directory called '"input file name"_results.csv.' Please look at this file for reference when working with your successful NT sequences generated from the input AA sequences. (please reference figure below).
The beauty of this program is its memoization. Whether the AA sequences were successful or a failure, the JSON file will append these results to a larger and larger dictionary for future use. It will never check the same two AA sequence set twice. The more this program is used, the more elaborate your dictionary becomes. 

For example, lets compare two 3' to 5' sequences:
1) AGCG
2) TGCG


Rule 1: check nucleotide for nucleotide directly comparing NT matches in the 3' to 5' direction. Note: we are not comparing complementary similarity. No more than two matches are allowed. Two offset comparisons are employed, as well. If any of these conditions show greater than two matches, the sequence (in its entirity) will be deemed void.

A) 3'- AGCG -5'<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; |  |  |  |<br>
     3'- TGCG -5'<br>
     3 of 4 NT match - no good [FAIL]

B) 3'- AGCG -5'
             |  |  |
     3'-    TGCG -5'
     0 of 4 NT match - good

C) 3'-    AGCG -5'
             |  |  |
     3'- TGCG    -5'
     0 of 4 NT match - good


Rule 2: Similar to rule 1, compare NT to NT comparison as is, however, convert the second sequence in comparison to a reverse complement. Aside from this caveat, all conditions are identical to Rule 1.

A) 3'- AGCG -5'
          |  |  |  |
     3'- CGCA -5'
     2 of 4 NT match - good

B) 3'- AGCG -5'
             |  |  |
     3'-    CGCA -5'
     0 of 4 NT match - good

C) 3'-    AGCG -5'
             |  |  |
     3'- CGCA    -5'
     0 of 4 NT match - good

For this example, Rule 2 did not fail, but Rule 1 did. If this was checked in the NT_checker.py script, a warning would be printed on the terminal stating Rule 1 was violated. 
If this set of sequence were to arise in the AA_to_NT_generator.py, the generated NT sequence will not be considered, and a new set of AA compatible NT sequences will be checked.

Rule 3: Palindromic sequences will be rejected. If checked in the NT_checker.py script, a warning will be printed onto the terminal. For AA_to_NT_generator.py, the generated palindromic NT sequence will not be considered and a new set of AA compatible NT sequences will be checked.

Palindromic sequence:

3'- AGCT -5', which in reverse complement, is:
3'- AGCT -5'

This is the summary for now.
