#Uploaded to GitHub 2019-05-30
#updated 2019-06-08
import os.path
csv_dir = '/Users/irahorecka/Desktop/Harddrive_Desktop/Python/GG Cloning/Junction Site Scripts/CSV'
json_dir = '/Users/irahorecka/Desktop/Harddrive_Desktop/Python/GG Cloning/Junction Site Scripts/JSON'
os.chdir(csv_dir)
import collections
import csv
import json
#look to replacing pandas with numpy - load time for import pandas is long.
import pandas as pd
import random
import time


rr_table = {
    'I': ['ATA','ATC','ATT'],
    'M': ['ATG'],
    'T': ['ACA','ACC','ACG','ACT'],
    'N': ['AAC','AAT'],
    'K': ['AAA','AAG'],
    'S': ['AGC','AGT'],
    'R': ['AGA','AGG'],
    'L': ['CTA','CTC','CTG','CTT'],
    'P': ['CCA','CCC','CCG','CCT'],
    'H': ['CAC','CAT'],
    'Q': ['CAA','CAG'],
    'R': ['CGA','CGC','CGG','CGT'],
    'V': ['GTA','GTC','GTG','GTT'],
    'A': ['GCA','GCC','GCG','GCT'],
    'D': ['GAC','GAT'],
    'E': ['GAA','GAG'],
    'G': ['GGA','GGC','GGG','GGT'],
    'S': ['TCA','TCC','TCG','TCT'],
    'F': ['TTC','TTT'],
    'L': ['TTA','TTG'],
    'Y': ['TAC','TAT'],
    '_': ['TAA','TAG','TGA'],
    'C': ['TGC','TGT'],
    'W': ['TGG']
}

# FUNCTIONS

#progress bar for micro_search_loop method in Search class
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 0, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()


#the almighty sequence checker - go with just a function for now
def seq_looker(junc_list):
    double_count = 0
    for i in junc_list:
        for j in junc_list:
            if i == j:
                double_count += 1
                continue
            junction = Junction(i,j)
            pass_fail, nt_fail = junction.seq_check()
            #palindromic looks at j seq to kill iteration early if found.
            if junction.palindromic() == 1:
                return 1
            elif pass_fail == 0.1:
                return 1
            elif pass_fail == 0.2:
                return 1
            else:
                pass
    if double_count > len(junc_list):
        return 1
    else:
        return 0



#CLASS

class CSV:
    def __init__(self, file):
        self.list = list()
        if file[-4:] != '.csv':
            self.file = file + '.csv'
        else:
            self.file = file

    def csv_read(self):
        self.csv_file = pd.read_csv(self.file, sep = ',')
        for index, row in self.csv_file.iterrows():
            if str(row['AA sequence']) == 'nan':
                continue
            else:
                self.list.append(str(row['AA sequence']))
        return self.list

    def csv_template(self, aa_list, scramble_list, results_list):
        res1 = (f"NH3,{','.join([i.upper() for i in aa_list])},COOH")
        res2 = (f"5',{','.join(scramble_list)},3'")
        res3 = (f"3',{','.join([Sequence(i).complement().upper() for i in scramble_list])},5'")
        csv_format = ['Amino Acid,6nt Sequence,4nt Junction Sequence\n']
        csv_format.extend(results_list)
        csv_format.append(',')
        csv_format.append(',')
        csv_format.append(res1)
        csv_format.append(',')
        csv_format.append(res2)
        csv_format.append(res3)
        return csv_format

    def csv_writer(self, results):
        new_file = self.file[:-4] + '_results.csv'
        with open(new_file,'w',newline = '') as csv_file:
            writer = csv.writer(csv_file, delimiter = ',')
            writer.writerows([i.split(',') for i in results])
        csv_file.close()


class Codon:
    def __init__(self, aa_list, dictionary, GC_bypass):
        self.aa_list = aa_list
        self.dictionary = dictionary
        self.GC_bypass = GC_bypass
        self.nt_list = list()
        self.nt_dict = dict()
        self.val_dict = dict()

    def codon_table(self):
        for aa in self.aa_list:
            for a in aa:
                try:
                    self.nt_list.append(self.dictionary[a.upper()])
                except KeyError:
                    pass
            chain_seq = [i + j for i in self.nt_list[0] for j in self.nt_list[1]]
            self.nt_dict[''.join(aa).upper()] = chain_seq
            self.nt_list = list()
        for key, value in self.nt_dict.items():
            for string in value:
                string1 = string[0:4]
                string2 = string[1:5]
                string3 = string[2:6]
                if self.GC_bypass == 0:
                    #add a GC check inactivator if called due to impossible combinations
                    self.nt_list.append(string1) if string1.count("G") + string1.count("C") < 3 and string1 != Sequence(string1).complement()[::-1] else None
                    self.nt_list.append(string2) if string2.count("G") + string2.count("C") < 3 and string2 != Sequence(string2).complement()[::-1] else None
                    self.nt_list.append(string3) if string3.count("G") + string3.count("C") < 3 and string3 != Sequence(string3).complement()[::-1] else None
                elif self.GC_bypass == 1:
                    self.nt_list.append(string1)
                    self.nt_list.append(string2)
                    self.nt_list.append(string3)
                if len(self.nt_list) != 0:
                    self.val_dict[string] = list(set(self.nt_list))
                self.nt_list = list()
            self.nt_dict[key] = self.val_dict
            self.val_dict = dict()
        return self.nt_dict


class Junction:
    def __init__(self, string, string2):
        self.string = string.upper()
        self.sequence = Sequence(self.string)
        self.string2 = string2.upper()
        self.sequence2 = Sequence(self.string2)
    #function with string input (4nt sequences only)
    #checks if sequence is palindromic - returns 1 if it is, 0 if it isn't
    def palindromic(self):
        pal_string2 = self.sequence2.r_complement()
        if self.string2 == pal_string2:
            return 1
        else:
            return 0
    #test rule1: that no two junction sites can have complement sequencec be complementary
    #test rule2: that no two junction sites can have reverse sequence be complementary
    def seq_check(self):
        counter1 = counter1_1 = counter1_2 = 0
        counter2 = counter2_1 = counter2_2 = 0
        for i in range(4): #maybe len(self.string)?
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
        if (counter1 > 2 or counter1_1 > 2 or counter1_2 > 2):
            return 0.1, max(counter1, counter1_1, counter1_2)
        elif (counter2 > 2 or counter2_1 > 2 or counter2_2 > 2):
            return 0.2, max(counter2, counter2_1, counter2_2)
        else:
            return 0, None


class Sequence:
    def __init__(self, string):
        self.string = string.upper()
        self.table = {
            'A':'T',
            'G':'C',
            'C':'G',
            'T':'A'
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


class Compute:
    def __init__(self, aa_list, product):
        self.aa_list = aa_list
        self.product = product

    def nt_string(self):
        string_list = list()
        for i in self.aa_list:
            codon_list = list()
            for key, value in self.product[i.upper()].items():
                codon_list.extend(value)
            string_list.append(set(codon_list))
        return string_list

    def scrambler(self, nt_list):
        self.nt_list = nt_list
        scramble_list = list()
        for i in range(len(self.aa_list)):
            scramble_list.append(random.sample(self.nt_list[i], 1)[0])
        return scramble_list

def seq_setter(aa_list):
    counts = collections.Counter(aa_list)
    sorted_list = sorted(aa_list, key=lambda x: -counts[x])
    return sorted_list

def sort_grouped(sort_list): 
    sort_set = sorted(set(sort_list))
    count_list = list()
    grouped_list = list()
    for i in sort_set:
        count_list.append(sort_list.count(i))
    index = 0
    for j in count_list:
        grouped_list.append([sort_set[index]] * j)
        index += 1
    return grouped_list


class Search:
    def __init__(self, codon_dict, search_no):
        self.codon_dict = codon_dict
        self.search_no = search_no

    def file_to_list(self):
        self.input_file = input('Please input .csv file name with GG cloning junction sequences: ')
        self.filename = CSV(self.input_file)
        self.aa_list = self.filename.csv_read()
        self.orig_list = self.aa_list

    def search_prime(self, GC_check):
        self.GC_check = GC_check
        self.product = Codon(self.aa_list, self.codon_dict, self.GC_check).codon_table()
        self.compute = Compute(self.aa_list, self.product)
        self.codon_tot = self.compute.nt_string()
    
    def search_loop(self):
        self.verif_set = set()
        self.redo_count = list()
        start_t = time.time()
        self.search_prime(0)
        self.stop_check = 1
        self.micro_stop_check = 1
        self.scramble_list = self.compute.scrambler(self.codon_tot)
        self.verif_set.add(tuple(self.scramble_list))
        test = seq_looker(self.scramble_list)
        tuple_count = 0

        while test != 0:
            if self.stop_check == 1:
                printProgressBar(self.stop_check, self.search_no, prefix = '  Progress:', suffix = 'Complete', length = 50)
            if self.stop_check == 50:
                self.search_prime(1)
            if self.stop_check > self.search_no:
                end_t = time.time()
                self.delta_t = end_t - start_t
                self.scramble_list = ["Impossible"]
                self.hit_dict = {''.join(self.aa_list) : self.scramble_list}
                break
            if self.stop_check == 10000:
                result = self.micro_search_loop()
                if result == "impossible":
                    self.aa_list = self.orig_list
                    self.scramble_list = ["Impossible"]
                    break
                else:
                    self.aa_list = self.orig_list
                    self.search_prime(1)
                    self.scramble_list = self.compute.scrambler(self.codon_tot)
                    self.stop_check += 1
                    continue
            self.scramble_list = self.compute.scrambler(self.codon_tot)
            if tuple(self.scramble_list) not in self.verif_set:
                tuple_count += 1
                self.verif_set.add(tuple(self.scramble_list))
                test = seq_looker(self.scramble_list)
            else:
                pass
            self.stop_check += 1
            printProgressBar(self.stop_check, self.search_no, prefix = '  Progress:', suffix = 'Complete', length = 50)

        if self.stop_check < self.search_no:
            end_t = time.time()
            self.delta_t = end_t - start_t
            self.hit_dict = {''.join(self.aa_list) : self.scramble_list}

    def micro_search_loop(self):
        temp_aalist = seq_setter(self.aa_list)
        self.grouped_list = sort_grouped(temp_aalist)
        for i in self.grouped_list:
            self.hit_dict = {}
            self.micro_verif_set = set()
            micro_tuple_count = 0
            count = i.count(i[0])
            if count == 1:
                continue
            self.aa_list = list()
            for j in i:
                if len(self.aa_list) == 0:
                    self.aa_list.append(j)
                    continue
                self.aa_list.append(j)
                self.search_prime(1)
                self.scramble_list = self.compute.scrambler(self.codon_tot)
                micro_test = seq_looker(self.scramble_list)
                os.chdir(json_dir)
                json_key = ''.join(self.aa_list).lower()
                with open('gg_ntdict.json') as json_dict:
                    master_dict = json.load(json_dict)
                    if json_key in master_dict:
                        if 'impossible' in [i.lower() for i in master_dict[json_key]]:
                            return 'impossible'
                        else:
                            self.scramble_list = master_dict[json_key]
                            print('\nEstablished sequence: miniloop')
                            print(f'\n{self.scramble_list}\n\n')
                            continue
                    else:
                        pass

                statcount_max = len(self.codon_tot[0])**len(self.codon_tot)
                while micro_test != 0:
                    if self.micro_stop_check == 1:
                        printProgressBar(self.micro_stop_check, self.search_no + 1, prefix = '  Progress:', suffix = 'Complete', length = 50)
                    if self.micro_stop_check == self.search_no:
                        return 'impossible'
                    self.scramble_list = self.compute.scrambler(self.codon_tot)
                    if tuple(self.scramble_list) not in self.micro_verif_set:
                        micro_tuple_count += 1
                        self.micro_verif_set.add(tuple(self.scramble_list))
                        if micro_tuple_count > statcount_max * 0.99:
                            return 'impossible'
                        micro_test = seq_looker(self.scramble_list)
                    else:
                        pass
                    
                    self.micro_stop_check += 1
                    printProgressBar(self.micro_stop_check, self.search_no + 1, prefix = '  Progress:', suffix = 'Complete', length = 50)
                
                if micro_test == 0:
                    self.hit_dict = {''.join(self.aa_list) : self.scramble_list}
                    with open('gg_ntdict.json', 'w') as json_dict:
                        for key, value in self.hit_dict.items():
                            master_dict[key] = value
                            json.dump(master_dict, json_dict)
        return "ok"

    def execute_script(self):
        os.chdir(json_dir)
        json_key = ''.join(self.aa_list).lower()
        results_list = list()
        with open('gg_ntdict.json') as json_dict:
            master_dict = json.load(json_dict)
            self.search_prime(1)
            if json_key in master_dict:
                self.scramble_list = master_dict[json_key]
                self.hit_dict = {}
                if self.scramble_list[0].lower() == 'impossible':
                    pass
                else:
                    for i in range(len(self.aa_list)):
                        for key, value in self.product[self.aa_list[i].upper()].items(): 
                            if self.scramble_list[i].upper() in key:
                                results_list.append(f'{self.aa_list[i].upper()},{key},{self.scramble_list[i]}')
                                break
                    os.chdir(csv_dir)
                    new_results_list = self.filename.csv_template(self.aa_list, self.scramble_list, results_list)
                    self.filename.csv_writer(new_results_list)
                    print('\nEstablished NT sequence:')
                    print(f'\n{self.scramble_list}\n\n')
                    print('Written to file\n')
            else:
                self.search_loop()
                if self.scramble_list[0].lower() == 'impossible':
                    with open('gg_ntdict.json', 'w') as json_dict:
                        for key, value in self.hit_dict.items():
                            master_dict[key] = value
                            json.dump(master_dict, json_dict)
                    self.hit_dict = {}
                    pass
                else:
                    for i in range(len(self.aa_list)):
                        for key, value in self.product[self.aa_list[i].upper()].items(): 
                            if self.scramble_list[i].upper() in key:
                                results_list.append(f'{self.aa_list[i].upper()},{key},{self.scramble_list[i]}')
                                break
                    new_results_list = self.filename.csv_template(self.aa_list, self.scramble_list, results_list)
                    self.filename.csv_writer(new_results_list)
        os.chdir(json_dir)
        if len(self.hit_dict) != 0:
            print(f'Operation time: {"%.2f" % self.delta_t} sec.')
            print(f'Number of sequence checks: {self.stop_check + self.micro_stop_check}\n')
            print('Success - check .csv file\n\n')
            with open('gg_ntdict.json', 'w') as json_dict:
                for key, value in self.hit_dict.items():
                    master_dict[key] = value
                    json.dump(master_dict, json_dict)
        else:
            try:
                print(f'Operation time: {"%.2f" % self.delta_t} sec.')
            except AttributeError:
                pass
            if self.scramble_list[0].lower() == 'impossible':
                print(f'\nImpossible AA sequence: \n\n{[i.upper() for i in self.aa_list]}')
                print('\n\nNothing Written to File\n')


# Execute script
#run_script(AA to codon table, # of sequence checks)
seq_check = 500000
run_script = Search(rr_table, seq_check)
run_script.file_to_list()
run_script.execute_script()
