#Uploaded to GitHub 2019-05-30
#A script to check junction site compatibility from a list of junction site sequences
#start with what you know from your OSX computer
import os.path
os.chdir('/Users/irahorecka/Desktop/Harddrive_Desktop/Python/GG Cloning/Junction Site Scripts/CSV')
import csv
import pandas as pd

class CSV:
    def __init__(self, file):
        self.file = file
        self.list = list()
    def csv_list(self):
        if self.file[-4:] != '.csv':
            self.file = self.file + '.csv'
        self.csv_file = pd.read_csv(self.file, sep = ',')
        for index, row in self.csv_file.iterrows():
            if str(row['NT sequence']) == 'nan':
                continue
            else:
                self.list.append(str(row['NT sequence']))
        return self.list

class Sequence:
    def __init__(self, string):
        self.string = string.lower()
        self.table = {
            'a':'t',
            'g':'c',
            'c':'g',
            't':'a'
        }
    #complementary sequencer
    def complement(self):
        comp = ''
        for i in self.string:
            comp += self.table.get(i, '0')
        return comp
    #reverse complementary sequencer - return reverse of complement sequence
    def r_complement(self):
        comp = ''
        for i in self.string:
            comp += self.table.get(i, '0')
        return comp[::-1]

class Junction:
    def __init__(self, string, string2):
        self.string = string.lower() #or .upper()
        self.sequence = Sequence(self.string)
        self.string2 = string2.lower()
        self.sequence2 = Sequence(self.string2)
    #function with string input (4nt sequences only)
    #checks if sequence is palindromic - returns 0 if it is, 1 if it isn't
    def palindromic(self):
        counter = 0
        pal_string = self.sequence.r_complement()
        for i in range(4):
            if self.string[i] == pal_string[i]:
                counter += 1
            if counter == 4:
                return 0
            else:
                return 1
    #test rule1: that no two junction sites can have complement sequencec be complementary
    #test rule2: that no two junction sites can have reverse sequence be complementary
    def seq_check(self):
        counter1 = counter1_1 = counter1_2 = 0
        counter2 = counter2_1 = counter2_2 = 0
        for i in range(4):
            #rule1
            if self.string[i] == self.string2[i]:
                counter1 += 1
            try:
                if self.string[i + 1] == self.string2[i]:
                    counter1_1 += 1
            except IndexError:
                pass
            try:
                if self.string[i] == self.string2[i + 1]:
                    counter1_2 += 1
            except IndexError:
                pass
            #rule2
            if self.string[i] == self.sequence2.r_complement()[i]:
                counter2 += 1
            try:
                if self.string[i + 1] == self.sequence2.r_complement()[i]:
                    counter2_1 += 1
            except IndexError:
                pass
            try:
                if self.string[i] == self.sequence2.r_complement()[i + 1]:
                    counter2_2 += 1
            except IndexError:
                pass
        '''print(counter1, counter1_1, counter1_2)
        print(counter2, counter2_1, counter2_2)
        print('------------------')'''
        if (counter1 > 2 or counter1_1 > 2 or counter1_2 > 2):
            return 0.1
        elif (counter2 > 2 or counter2_1 > 2 or counter2_2 > 2):
            return 0.2
        else:
            return 0

#the almighty sequence checker - go with just a function for now
def seq_looker(junc_list):
    gc_set = set()
    rule1 = 'The following combination(s) violates rule 1:\n'
    sub_rule1 = ''
    rule2 = 'The following combination(s) violates rule 2:\n'
    sub_rule2 = ''
    for i in junc_list:
        gc_count = 0
        for n in i:
            if n.upper() == 'G':
                gc_count += 1
            elif n.upper() == 'C':
                gc_count += 1
        if gc_count > 2:
            gc_set.add(i.upper())
        for j in junc_list:
            junction = Junction(i,j)
            pass_fail = junction.seq_check()
            if junction.palindromic() == 0:
                return 'Palindrome found: ' + i.upper()
            if i == j:
                pass
            elif pass_fail == 0.1:
                sub_rule1 += f'{i.upper()} {j.upper()}\n'
            elif pass_fail == 0.2:
                sub_rule2 += f'{i.upper()} {j.upper()}\n'
            else:
                pass
    sub_rule1 = sub_rule1[:-(int(len(sub_rule1)/2))]
    sub_rule2 = sub_rule2[:-(int(len(sub_rule2)/2))]
    if len(sub_rule1) == 0:
        rule1 = rule1 + 'None'
    if len(sub_rule2) == 0:
        rule2 = rule2 + 'None'
    if len(sub_rule1) + len(sub_rule2) + len(gc_set) == 0:
        print(f'\nInput junction sequences: \n{junc_list}\n')
        print('\nAll rules satisfied.\n')
    else:
        print(f'\nInput junction sequences: \n{junc_list}\n')
        print('\n' + (rule1 + sub_rule1) + '\n' + (rule2 + sub_rule2) + '\n')
        print('>2 GC content: ', list(gc_set), '\n')


input_file = input('Please input .csv file name with GG cloning junction sequences: ')
filename = CSV(input_file)
test_list = filename.csv_list()
seq_looker(test_list)
