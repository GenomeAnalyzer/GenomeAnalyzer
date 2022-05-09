# coding: utf8
import DNA_mod
# import set
import array
import glob
import moduleDNA as m
import os
import getopt
import sys


short_options = "hos:"
long_options = ["help", "output", "sequences"]
full_cmd_arguments = sys.argv
argument_list = full_cmd_arguments[1:]

find = False
output = 0
help = False
try:
        arguments, values = getopt.getopt(
            argument_list, short_options, long_options)
except getopt.error as err:
        print(str(err))
        sys.exit(2)
for current_argument, current_value in arguments:
        if current_argument in ("-h", "--help"):
            print("usages:\n\t-h/--help  help\n\t-o/--output the program will output htmls \n\t-s/--sequences set how many sequences will be analyzed (without = max) ")
            help = True
        elif current_argument in ("-s", "--sequences"):
            try:
                if not current_value.isnumeric():
                    raise NameError('nan')
                else:
                    fin = current_value
                    find = True
            except NameError:
                print("Arg is not a number")
                raise
        elif current_argument in ("-o", "--output"):
            output = 1
if find == False:
        fin = 25698
if help == False:
    DNA_mod.py_launch(output,fin)
#    DNA_bin.launch(mpi4py.MPI.COMM_WORLD)

#    if __name__ == "__main__":
#        main()