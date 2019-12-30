#!/usr/bin/env python
import sys
path='./'
inputfile=sys.argv[1]

#Get input from input file "Heteromirpred/user_input.fasta"
openfile=open(path+inputfile,'r')

readfile=openfile.read()
openfile.close()

#Count inputs by '>'
number_input=readfile.count('>')
#print "number of input ="+str(number_input)
all_inputs=readfile.split('\n')
#print all_inputs
#print val1

#Get input from HeteroMirPred output "out1.txt"
read_output=open(path+'/'+inputfile+'.out','r')
output_file=read_output.read()
read_output.close()
output_file=output_file.split('\n')
#for i in output_file:
#	print i+'\n'
#print output_file
del output_file[0:5]
#print output_file
arrange_output=[]
for n in range (0,number_input):
	edit1=output_file[n].split()
	edit2=edit1[2].replace("'",'')
	edit3=edit2.split(":")
	arrange_output.append(edit3[1])
	if edit1[3]=='+':
		arrange_output.append(edit1[4])
	else:
		arrange_output.append(edit1[3])

input_list = ["Seq_number","Seq_name","Prediction","Probability"]
for i in range(1,number_input+1):
	i=i*2
	input_list.append(i/2)
	name_seq=all_inputs[i-2].replace('\r','')
	input_list.append(name_seq)
	input_list.append(arrange_output[i-2])
	input_list.append(arrange_output[i-1])

def print_table_in_file(data, col_len):
    counter = 1
    #print str(data)
    for element in data:
        print str(element)+'\t',
        if counter % col_len == 0:
            print ''
	
        counter += 1
    
    
print_table_in_file(input_list,4)
