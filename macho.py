#!/usr/bin/env python3
# AUTHORS of Multiple Alignment Column Hits Observer (MACHO)
# Daniil Nagornyi
# Vsevolod Maslenikov
# Vitalii Gagarochkin

import argparse

def get_seqs(path):

    seqs = dict()

    with open(path, "r") as input:
        for line in input:
            if line[0] == ">":
                name = line.strip()
                seqs[name] = ""
                continue
            seqs[name] += line.strip()
    return seqs

def sort_seqs_id(seq_dict):
    return {k: v for k, v in sorted(seq_dict.items(), key=lambda item: item[0])}

def numerate_seq(seq_dict):
    for name,value in seq_dict.items():
        seq = value
        seq = list([x for x in seq])
        counter = 0
        for id in range(len(seq)):
            letter =  seq[id]
            if letter != '-':
                seq[id] = counter
                counter += 1
            else:
                continue
        seq_dict[name] = seq
    return seq_dict
            
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
        if matched_columns[column_id][0] + 1 == matched_columns[column_id+1][0]:
            continue
        else:
            areas_1.append(area_1)
            area_1 = list()
            areas_2.append(area_2)
            area_2 = list()
            
    return areas_1, areas_2

parser = argparse.ArgumentParser(description = "Comparison of two multiple alignments of the same set of sequences")
parser.add_argument('alignment_1', type = str, help = 'The path to the file with the first alignment in FASTA format')
parser.add_argument('alignment_2', type = str, help = 'The path to the file with the second alignment in FASTA format')
parser.add_argument('out', type = str, help = 'The path to the file for recording the results in TSV format')
parser.add_argument('-hr', '--human-readable', action = 'store_true', help = 'Group matching columns into matching blocks in the output file')
parser.add_argument('-g', '--gaps-controller', action = 'store_true', help = 'Add the ability to specify the maximum number of gaps in a column')
args = parser.parse_args()

seq_path_first = args.alignment_1
seq_path_second = args.alignment_2
seq_path_out = args.out

seq_dict_first = get_seqs(seq_path_first)
seq_dict_second = get_seqs(seq_path_second)

max_gaps = len(seq_dict_first)

if args.gaps_controller:
    print(f'Number of sequences in alignment: {len(seq_dict_first)}')
    max_gaps = int(input('Maximum number of gaps in a column: '))

seq_dict_first = sort_seqs_id(seq_dict_first)
seq_dict_second = sort_seqs_id(seq_dict_second)

seq_dict_first = numerate_seq(seq_dict_first)
seq_dict_second = numerate_seq(seq_dict_second)

columns_first = get_columns(seq_dict_first)
columns_second = get_columns(seq_dict_second)

#*@@@@@@@@@@@@@@@@@@@@@@
matched_columns = list()
#*@@@@@@@@@@@@@@@@@@@@@@

#* находим совпадения
for i in range(len(columns_first)):
    colum = columns_first[i]
    for j in range(len(columns_second)):
        if colum == columns_second[j] and colum.count('-') <= max_gaps:
            #* где i - позиция из 1 файла, а j - позиция из 2 файла
            match_position = tuple([i+1,j+1])
            matched_columns.append(match_position)


#* выводим совпадения в новый файл
areas_1, areas_2 = process_output(matched_columns)

with open(seq_path_out, "w") as output:

    if args.human_readable:
        output.write("Block\tAlignment_1\tAlignment_2\n")
        for i in range(len(areas_1)):
            block = i + 1
            seq1_start = areas_1[i][0]
            seq1_end = areas_1[i][-1]
            seq2_start = areas_2[i][0]
            seq2_end = areas_2[i][-1]
            output.write(f"{block}\t{seq1_start}-{seq1_end}\t{seq2_start}-{seq2_end}\n")
    else:
        output.write("Alignment_1\tAlignment_2\n")
        for i in range(len(areas_1)):
            for j in range(len(areas_1[i])):
                output.write(f"{areas_1[i][j]}\t{areas_2[i][j]}\n")


print(f'First alignment length: {len(list(seq_dict_first.values())[0])}')
print(f'First alignment length: {len(list(seq_dict_second.values())[0])}')
print(f'Percentage of matching columns for the first alignment: {(len(matched_columns) - 1) / len(list(seq_dict_first.values())[0]) * 100:.02f} %')
print(f'Percentage of matching columns for the second alignment: {(len(matched_columns) - 1) / len(list(seq_dict_second.values())[0]) * 100:.02f} %')