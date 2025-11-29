# coding: utf-8
import json
import pickle
import sys

if len(sys.argv) != 3:
    print "Usage: python json_to_pickle.py <input_json_file> <output_pickle_file>"
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

try:
    with open(input_file, 'r') as f_in:
        data = json.load(f_in)

    with open(output_file, 'wb') as f_out:
        pickle.dump(data, f_out)

    print "Successfully converted %s to %s (preserving dictionary structure)" % (input_file, output_file)

except IOError as e:
    print "Error: Could not open or write to file: %s" % e
    sys.exit(1)
except ValueError as e:
    print "Error: Invalid JSON format in: %s" % e
    sys.exit(1)
except Exception as e:
    print "An error occurred: %s" % e
    sys.exit(1)