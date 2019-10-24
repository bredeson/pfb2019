#!/usr/bin/env python3

import os
import sys

# It's good practice to write a command-line usage message.
# that we can re-use throughout the code. Here the `usage()`
# function can take an optional argument to print an error
# message to inform the user they did something wrong.
def usage(message=None):
    message = '' if message is None else 'ERROR: {}\n\n'.format(message)
    sys.stderr.write("""
Usage: {} <input.tsv> <column-number>

Notes: column-number is 1-based, so the first column is `1'


{}""".format(os.path.basename(sys.argv[0]), message))
    sys.exit(1)
    

def main(argv):
    if len(argv) == 0:
        usage()
        
    elif len(argv) != 2:
        usage("Unexpected number of arguments")

    try:
        data_file = open(argv[0], 'r')
        data_column = int(argv[1])
    except FileNotFoundError as error:
        usage(str(error))
    except ValueError:
        usage("Invalid column number, positive integer expected")

    frequencies = {}
    for line in data_file:
        if line.isspace() or \
           line.startswith('#'):
            continue

        fields = line.strip().split('\t')

        try:
            value = int(fields[data_column-1])
        except IndexError:
            usage("Invalid column number, too few fields in file")
        except ValueError:
            usage("Invalid column number, positive integer expected")

        if value not in frequencies:
            frequencies[value]  = 1
        else:
            frequencies[value] += 1

    pdensity = 0.0
    cdensity = 0.0
    sum_freq  = sum(frequencies.values())
    max_freq  = max(frequencies.values())
    max_value = max(frequencies.keys())
    max_value_length = len(str(max_value))
    print("{1:>{0}}|{2:>5}|{3:>5}|".format(
        max_value_length, 'x', 'PDF','CDF'))
    for value in range(max_value + 1):
        try:
            freq = frequencies[value]
        except KeyError:
            freq = 0

        pdensity = float(freq) / sum_freq
        cdensity += pdensity
        
        print("{1:>{0}}|{2:.3f}|{3:.3f}|{4:}".format(
            max_value_length,
            value,
            pdensity,
            cdensity,
            ']' * int((80 - max_value_length - 13) / max_freq * freq)))
    
            
if __name__ == '__main__':
    main(sys.argv[1:])  # remove program name from argument list
