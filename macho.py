#!/usr/bin/env python3
# AUTHORS of Multiple Alignments Column Hits Observer (MACHO)
# Daniil Nagornyi
# Vsevolod Maslenikov
# Vitalii Gagarochkin
#v1.4

import os
import sys
import argparse

def get_seqs(path):

    seqs = dict()

    with path as input:
        for line in input:
            if line[0] == ">":
                name = line.strip()
                seqs[name] = ""
                continue
            seqs[name] += line.strip()
    return seqs

def sort_seqs_id(seq_dict):
    names = sorted(seq_dict.keys())
    new_dict = dict()
    for i in range(len(names)):
        name = names[i]
        new_dict[i] = seq_dict[name]
    return new_dict

def find_overlap(seq_1, seq_2):
    counter_1 = 0
    counter_2 = len(seq_1)
    if seq_1 in seq_2:
        counter_1 = seq_2.find(seq_1)
        counter_2 = 0
    elif seq_2 in seq_1:
        counter_1 = 0
        counter_2 = seq_1.find(seq_2)
    else:
        max_overlap = min(len(seq_1), len(seq_2)) - 1
        for o in range(max_overlap, max_overlap//5, -1):
            if seq_1.endswith(seq_2[:o]):
                counter_1 = 0
                counter_2 = len(seq_1) - o
                break
            elif seq_2.endswith(seq_1[:o]):
                counter_1 = len(seq_2) - o
                counter_2 = 0
                break
            else:
                continue
    return counter_1, counter_2

def numerate_seq(seq_dict_1, seq_dict_2):
    for name in seq_dict_1.keys():
        seq_1 = seq_dict_1[name]
        seq_2 = seq_dict_2[name]
        seq_1_ungap = seq_1.replace('-','')
        seq_2_ungap = seq_2.replace('-','')
        if seq_1_ungap == seq_2_ungap:
            counter_1 = 0
            counter_2 = 0
        else:
            counter_1, counter_2 = find_overlap(seq_1_ungap, seq_2_ungap)
        seq_1 = list([x for x in seq_1])
        seq_2 = list([x for x in seq_2])
        for id in range(len(seq_1)):
            letter =  seq_1[id]
            if letter != '-':
                seq_1[id] = counter_1
                counter_1 += 1
            else:
                continue
        seq_dict_1[name] = seq_1
        for id in range(len(seq_2)):
            letter =  seq_2[id]
            if letter != '-':
                seq_2[id] = counter_2
                counter_2 += 1
            else:
                continue
        seq_dict_2[name] = seq_2
    return seq_dict_1, seq_dict_2

def numerate_seq_smart(seq_dict_1, seq_dict_2):
    for name in seq_dict_1.keys():
        seq_1 = seq_dict_1[name]
        seq_2 = seq_dict_2[name]
        seq_1_ungap = seq_1.replace('-','')
        seq_2_ungap = seq_2.replace('-','')
        if seq_1_ungap == seq_2_ungap:
            counter_1 = 0
            counter_2 = 0
        else:
            counter_1, counter_2 = find_overlap(seq_1_ungap, seq_2_ungap)
        seq_1 = list([x for x in seq_1])
        seq_2 = list([x for x in seq_2])
        for id in range(len(seq_1)):
            letter =  seq_1[id]
            if letter != '-':
                seq_1[id] = counter_1
                counter_1 += 1
            else:
                seq_1[id] += str(counter_1)
        seq_dict_1[name] = seq_1
        for id in range(len(seq_2)):
            letter =  seq_2[id]
            if letter != '-':
                seq_2[id] = counter_2
                counter_2 += 1
            else:
                seq_2[id] += str(counter_2)
        seq_dict_2[name] = seq_2
    return seq_dict_1, seq_dict_2

def get_columns(seq_dict):
    columns = list()
    length = len(list(seq_dict.values())[0] )
    for i in range(length):
        column = list()
        for name, value in seq_dict.items():
            column.append(value[i])
        columns.append(column)
    return columns

def process_output(matched_columns):

    areas_1 = list()
    area_1 = list()
    areas_2 = list()
    area_2 = list()
    matched_columns.append((0,0))

    for column_id in range(len(matched_columns)-1):
        area_1.append(matched_columns[column_id][0])
        area_2.append(matched_columns[column_id][1])
        if matched_columns[column_id][0] + 1 == matched_columns[column_id+1][0] and abs(matched_columns[column_id][1] - matched_columns[column_id+1][1]) == 1:
            continue
        else:
            areas_1.append(area_1)
            area_1 = list()
            areas_2.append(area_2)
            area_2 = list()
            
    return areas_1, areas_2



parser = argparse.ArgumentParser(description = 'MACHO: Multiple Alignments Column Hits Observer\nComparison of two multiple alignments of the same set of sequences', epilog = 'For more information, see the web documentation (on Russian): https://kodomo.fbb.msu.ru/~vitalii.g/term2/MACHO.html\nWe hope that you will enjoy using our tool', formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument('alignment_1', type = argparse.FileType('r'), help = 'The path to the file with the first alignment in FASTA format')
parser.add_argument('alignment_2', type = argparse.FileType('r'), help = 'The path to the file with the second alignment in FASTA format')
parser.add_argument('out', nargs = '?', type = argparse.FileType('w'), default = sys.stdout, help = 'The path to the file for recording the results in TSV format')
parser.add_argument('-hr', '--human-readable', action = 'store_true', help = 'Group matching columns into matching blocks in the output file')
parser.add_argument('-g', '--gaps-controller', action = 'store_true', help = 'Add the ability to specify the maximum number of gaps in a column')
parser.add_argument('-s', '--smart-mode', action = 'store_true', help = 'Consider the differences between indels')
args = parser.parse_args()


seq_path_first = args.alignment_1
seq_path_second = args.alignment_2
seq_path_out = args.out

seq_dict_first = get_seqs(seq_path_first)
seq_dict_second = get_seqs(seq_path_second)

max_gaps = len(seq_dict_first) - 1

if args.gaps_controller:
    print(f'Number of sequences in alignment: {len(seq_dict_first)}')
    max_gaps = int(input('Maximum number of gaps in a column: '))

seq_dict_first = sort_seqs_id(seq_dict_first)
seq_dict_second = sort_seqs_id(seq_dict_second)


if args.smart_mode:
    seq_dict_first, seq_dict_second = numerate_seq_smart(seq_dict_first, seq_dict_second)
else:
    seq_dict_first, seq_dict_second = numerate_seq(seq_dict_first, seq_dict_second)

columns_first = get_columns(seq_dict_first)
columns_second = get_columns(seq_dict_second)

#*@@@@@@@@@@@@@@@@@@@@@@
matched_columns = list()
#*@@@@@@@@@@@@@@@@@@@@@@

#* находим совпадения
for i in range(len(columns_first)):
    colum = columns_first[i]
    for j in range(len(columns_second)):
        if colum == columns_second[j] and ''.join(list(map(str, colum))).count('-') <= max_gaps:
            #* где i - позиция из 1 файла, а j - позиция из 2 файла
            match_position = tuple([i+1,j+1])
            matched_columns.append(match_position)


#* выводим совпадения в новый файл
areas_1, areas_2 = process_output(matched_columns)


print(f'First alignment length ({os.path.basename(seq_path_first.name)}): {len(list(seq_dict_first.values())[0])}')
print(f'Second alignment length ({os.path.basename(seq_path_second.name)}): {len(list(seq_dict_second.values())[0])}')
print(f'Percentage of matching columns for the first alignment ({os.path.basename(seq_path_first.name)}): {(len(matched_columns) - 1) / len(list(seq_dict_first.values())[0]) * 100:.02f} %')
print(f'Percentage of matching columns for the second alignment ({os.path.basename(seq_path_second.name)}): {(len(matched_columns) - 1) / len(list(seq_dict_second.values())[0]) * 100:.02f} %')


with seq_path_out as output:

    if args.human_readable:
        output.write(f"Block\tAlignment_1 ({os.path.basename(seq_path_first.name)})\tAlignment_2 ({os.path.basename(seq_path_second.name)})\n")
        for i in range(len(areas_1)):
            block = i + 1
            seq1_start = areas_1[i][0]
            seq1_end = areas_1[i][-1]
            seq2_start = areas_2[i][0]
            seq2_end = areas_2[i][-1]
            output.write(f"{block}\t{seq1_start}-{seq1_end}\t{seq2_start}-{seq2_end}\n")
    else:
        output.write(f"Alignment_1 ({os.path.basename(seq_path_first.name)})\tAlignment_2 ({os.path.basename(seq_path_second.name)})\n")
        for i in range(len(areas_1)):
            for j in range(len(areas_1[i])):
                output.write(f"{areas_1[i][j]}\t{areas_2[i][j]}\n")