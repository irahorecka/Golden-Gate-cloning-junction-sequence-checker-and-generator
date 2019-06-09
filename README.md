**README for Golden Gate Cloning Junction Sequence Checker and Generator**

Please read the following article for a precis of [Golden Gate cloning](https://en.wikipedia.org/wiki/Golden_Gate_Cloning).

This is a two-tool program that:
1) Searches if the provided 4nt sequence (i.e. type iis re that leaves a 4nt overhang sequence) is compatible with one another in an in-silico cloning design.
2) Creates a set of compatible 4nt overhang sequence dependent on a sequence of 2AA sequences provided. The latter is useful for generating compatible 4nt sequences for junction sequences that one desires to dictate the AA used in the junction sites.

Below is a common list of type iis re that leaves a 4nt overhang sequence:
- [BsaI](https://www.neb.com/products/r0535-bsai#Product%20Information)
- [BbsI](https://www.neb.com/products/r0539-bbsi#Product%20Information)
- [BsmbI](https://www.neb.com/products/r0580-bsmbi#Product%20Information)

This program is taylored for [Type IIS restriction enzyme](https://www.thermofisher.com/us/en/home/life-science/cloning/restriction-enzyme-digestion-and-ligation/restriction-enzyme-cloning/anza-restriction-enzyme-system/type-iis-restriction-enzymes.html) that leaves a 4nt overhang sequence. There are plans on furthering the capability of this program by encompassing Type IIS restriction enzymes that leave a 3nt overhang sequences and overhang sequences >4nt in length.

Download or copy the code into your favorite text editor. Please check the directories assigned to variables csv_dir and json_dir written in this program; these is the only variables that must be modified to successfully run the code. There is a directory for your CSV, JSON, and python file. Please set var csv_dir and json_dir to your CSV and JSON directories, respectively.

* *CSV:* * Please format your CSV as the following. <br>
- The column headers is not case sensitive, however col 1 must be names 'NT sequence' (for NT_checker.py script) and col 2 must be named 'AA sequence' (for AA_to_NT_generator.py script). 
- In your CSV file, each row corrensponds to individual sequences (whether NT or AA) that is to be checked. 
    - Please list all NT sequences in column one. The case is not sensitive, but make sure to write a 4nt sequence per row. At least two rows are required, or an exception will be thrown.
    - For the AA column, please write your AA sequence as two single-letter AA code. Two and only two AA sequence must be written per row. The cells are not case sensitive. More than two rows are required to run the program, or an execption will be thrown. (please reference figure below)

* *JSON:* * Please download the JSON file 'gg_ntdict.json.' Do not change the name of this file or an exception will be thrown. However, if you are inclined to change the name of this file, you are welcome to do so. Please make sure to change every instance to the file in the program to your new name. Store this file appropriately in your JSON directory.

* *Running the program:* * Please run the program in your terminal. Change directories to your script(s)' location. Run AA_to_NT_generator.py script by typing 'python AA_to_NT_generator.py' in your terminal. Run NT_checker.py script by typing 'python NT_checker.py' in your terminal.

* *NT_checker.py:* * This script will look at content in col1 (where only NT sequences are listed). A detailed result will print onto your terminal. If the NT sequence fails, reasons for failure will be outlined on the terminal. 

* *AA_to_NT_generator.py:* * This script will look at content in col2 (where only AA sequences are listed). A load bar and script status will be printed onto your terminal. If the script completed successfully, you will find a summary .csv file written to your CSV directory called '"input file name"_results.csv.' Please look at this file for reference when working with your successful NT sequences generated from the input AA sequences (please reference figure below).

The beauty of the AA_to_NT_generator.py script is its memoization. Whether the AA sequence checks were successful or a failure, the JSON file will append these results to an expanding dictionary. When initianting the script, the script will first check the JSON file for the input AA sequences, i.e. it will never check the same two AA sequence set twice. The more the script is run with different AA combinations, the more elaborate your dictionary becomes.

* *Criteria for Selection* *

For example, lets compare two 3' to 5' sequences:
1) AGCG
2) TGCG


Rule 1: check nucleotide for nucleotide directly comparing NT matches in the 3' to 5' direction. Note: we are not comparing complementary similarity. No more than two matches are allowed. Two offset comparisons are employed, as well. If any of these conditions show greater than two matches, the sequence (in its entirity) will be deemed void.

A) 3'- AGCG -5'<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; |  |  |  |<br>
    &nbsp;&nbsp;&nbsp;&nbsp; 3'- TGCG -5'<br>
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3 of 4 NT match - no good [FAIL]<br>

B) 3'- AGCG -5'<br>
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;           |  |  |<br>
      &nbsp;&nbsp;3'-    TGCG -5'<br>
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0 of 4 NT match - good<br>

C) 3'-    AGCG -5'<br>
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;          |  |  |<br>
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3'- TGCG    -5'<br>
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0 of 4 NT match - good<br>


Rule 2: Similar to rule 1, compare NT to NT comparison as is, however, convert the second sequence in comparison to a reverse complement. Aside from this caveat, all conditions are identical to Rule 1.

A) 3'- AGCG -5'<br>
       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;   |  |  |  |<br>
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3'- CGCA -5'<br>
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2 of 4 NT match - good<br>

B) 3'- AGCG -5'<br>
        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      |  |  |<br>
      &nbsp;&nbsp;3'-    CGCA -5'<br>
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0 of 4 NT match - good<br>

C) 3'-    AGCG -5'<br>
         &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;     |  |  |<br>
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3'- CGCA    -5'<br>
       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0 of 4 NT match - good<br>

For this example, Rule 2 did not fail, but Rule 1 did. If this was checked in the NT_checker.py script, a warning would be printed on the terminal stating Rule 1 was violated. 
If this set of sequence were to arise in the AA_to_NT_generator.py, the generated NT sequence will not be considered, and a new set of AA compatible NT sequences will be checked.

Rule 3: Palindromic sequences will be rejected. If checked in the NT_checker.py script, a warning will be printed onto the terminal. For AA_to_NT_generator.py, the generated palindromic NT sequence will not be considered and a new set of AA compatible NT sequences will be checked.

* *Palindromic sequence:* *

3'- AGCT -5', which in reverse complement, is: <br>
3'- AGCT -5'

This is the summary for now.
